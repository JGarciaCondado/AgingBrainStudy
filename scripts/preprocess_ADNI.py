import pandas as pd
import numpy as np

# Load baseline and longitudinal data
df_baseline = pd.read_csv('data/final/adni/adni_baseline_040725.csv', 
                          index_col='ID',
                          usecols=lambda x: x != 'Unnamed: 0')
df_long = pd.read_csv('data/final/adni/adni_long_040725.csv',
                      index_col='ID',
                      usecols=lambda x: x != 'Unnamed: 0')

# Remove baseline duplicated participants
df_baseline = df_baseline[~df_baseline.index.duplicated(keep='first')]

# Keep only participants with more than two time points
df_long = df_long.sort_values(by=['ID', 'Date'])
df_long['nTimePoints'] = df_long.groupby('ID')['zPACC'].transform(lambda x: x.notna().sum())
df_long = df_long[df_long['nTimePoints'] >= 2]
long_ids = df_long.index.unique()
df_baseline = df_baseline[df_baseline.index.isin(long_ids)]

# Extract followup time of subjects from longitudinal data
max_followup = df_long.groupby('ID')['MonthsFromBaseline_raw'].max().to_dict()
df_baseline['follow_up_time'] = df_baseline.index.map(max_followup)/12

# Process MRI features for BrainAge modelling by keeping only columns that start with MRI_
df_structural = df_baseline.filter(like='MRI_').copy()
df_structural.columns = df_structural.columns.str.replace('MRI_FS7_rnr_', '')

# Combine those with lh and rh cortical thickness
lh = df_structural.filter(like='lh_').columns
rh = df_structural.filter(like='rh_').columns
for l, r in zip(lh, rh):
    df_structural['bi_' + l.split('_')[1]] = df_structural[l] + df_structural[r]
    df_structural.drop(columns=[l, r], inplace=True)

# Combine those with Left and Right Hippocampus and Amygdala
left = df_structural.filter(like='Left').columns
right = df_structural.filter(like='Right').columns
for l, r in zip(left, right):
    df_structural['bi_' + l.split('_')[1]] = df_structural[l] + df_structural[r]
    df_structural.drop(columns=[l, r], inplace=True)

# Normalize by total intracranial volume with respect to CN
def normalize_mri(mri_col, icv_col, mask):
    from sklearn.linear_model import LinearRegression
    model = LinearRegression().fit(icv_col[mask].values.reshape(-1, 1), mri_col[mask])
    b = model.coef_[0]
    c = icv_col[mask].mean()
    adj_mri = mri_col - b * (icv_col - c)
    return adj_mri

mask_cn = df_baseline['Diagnosis'] == 'CN'
df_structural['bi_Hippocampus'] = normalize_mri(df_structural['bi_Hippocampus'], df_structural['ICV_vol'], mask_cn)
df_structural['bi_Amygdala'] = normalize_mri(df_structural['bi_Amygdala'], df_structural['ICV_vol'], mask_cn)

# Keep only MRI_Age and bi_ columns
df_structural = df_structural[['MRI_Age'] + [col for col in df_structural.columns if 'bi_' in col]]

# Rename MRI_Age to age and lower case column names
df_structural = df_structural.rename(columns={'MRI_Age': 'age'})
df_structural.columns = df_structural.columns.str.lower()

# Save as csv
df_structural.to_csv('data/final/adni/processed/structural_features.csv')

# Create clinical file with type of Apoe4
df_clinical = df_baseline[['e4_carrier', 'AMYLOID_STATUS', 'Diagnosis']].copy()
df_clinical['CN'] = ((df_clinical['e4_carrier'] == 0.0) & (df_clinical['AMYLOID_STATUS'] == 0.0) & (df_clinical['Diagnosis'] == 'CN')).astype(int)
df_clinical['e4+'] = ((df_clinical['e4_carrier'] == 1.0) & (df_clinical['AMYLOID_STATUS'] == 0.0) & (df_clinical['Diagnosis'] == 'CN')).astype(int)
df_clinical['e4+ab+'] = ((df_clinical['e4_carrier'] == 1.0) & (df_clinical['AMYLOID_STATUS'] == 1.0) & (df_clinical['Diagnosis'] == 'CN')).astype(int)
df_clinical['ab+'] = ((df_clinical['e4_carrier'] == 0.0) & (df_clinical['AMYLOID_STATUS'] == 1.0) & (df_clinical['Diagnosis'] == 'CN')).astype(int)
df_clinical['MCI'] = ((df_clinical['e4_carrier'] == 0.0) & (df_clinical['AMYLOID_STATUS'] == 0.0) & (df_clinical['Diagnosis'] == 'MCI')).astype(int)
df_clinical['MCIe4+'] = ((df_clinical['e4_carrier'] == 1.0) & (df_clinical['AMYLOID_STATUS'] == 0.0) & (df_clinical['Diagnosis'] == 'MCI')).astype(int)
df_clinical['MCIe4+ab+'] = ((df_clinical['e4_carrier'] == 1.0) & (df_clinical['AMYLOID_STATUS'] == 1.0) & (df_clinical['Diagnosis'] == 'MCI')).astype(int)
df_clinical['MCIab+'] = ((df_clinical['e4_carrier'] == 0.0) & (df_clinical['AMYLOID_STATUS'] == 1.0) & (df_clinical['Diagnosis'] == 'MCI')).astype(int)
df_clinical['CNunknown'] = ((df_clinical['e4_carrier'].isna() | df_clinical['AMYLOID_STATUS'].isna()) & (df_clinical['Diagnosis'] == 'CN')).astype(int)
df_clinical['MCIunknown'] = ((df_clinical['e4_carrier'].isna() | df_clinical['AMYLOID_STATUS'].isna()) & (df_clinical['Diagnosis'] == 'MCI')).astype(int)
df_clinical.drop(columns=['e4_carrier', 'AMYLOID_STATUS', 'Diagnosis'], inplace=True)
df_clinical.to_csv('data/final/adni/processed/clinical.csv')

# Process baseline so that across all cohorts the same type of naming and conventions

# Date of each point of measurment. They were already pulled to be closest in time to first known MRI
df_baseline['MRI_SessionDate'] = pd.to_datetime(df_baseline['MRI_SessionDate'], errors='coerce')
df_baseline['Clinical_Date'] = pd.to_datetime(df_baseline['Clinical_Date'], errors='coerce')
df_baseline['AB_SCANDATE'] = pd.to_datetime(df_baseline['AB_SCANDATE'], errors='coerce')
df_baseline['PLASMA_DATE'] = pd.to_datetime(df_baseline['PLASMA_DATE'], errors='coerce')
df_baseline['TAU_SCANDATE'] = pd.to_datetime(df_baseline['TAU_SCANDATE'], errors='coerce')


# Time differences
df_baseline['time_diff_pacc'] = (df_baseline['Clinical_Date'] - df_baseline['MRI_SessionDate']).dt.days/365.25
df_baseline['time_diff_ab'] = (df_baseline['AB_SCANDATE'] - df_baseline['MRI_SessionDate']).dt.days/365.25
df_baseline['time_diff_ptau'] = (df_baseline['PLASMA_DATE'] - df_baseline['MRI_SessionDate']).dt.days/365.25
df_baseline['time_diff_tau'] = (df_baseline['TAU_SCANDATE'] - df_baseline['MRI_SessionDate']).dt.days/365.25

# Calculate Tau composite by averaging over regions with _bh
# Starts with PVC_CTX and finishes with _SUVR_TAU
tau_features = [col for col in df_baseline.columns if 'PVC_CTX' in col and '_SUVR_TAU' in col]
df_baseline['tau_composite'] = df_baseline[tau_features].mean(axis=1)

# Make ptau217 values below 0 to NaN
df_baseline.loc[df_baseline['pT217_AB42_F'] < 0, 'pT217_AB42_F'] = np.nan

# Rename of columns of interest
# Ptau use pTau217 to AB42 ratio
df_baseline.rename(columns={
    'MRI_Age': 'mri_age',
    'Sex': 'sex',
    'Education': 'edu',
    'Diagnosis': 'diagnosis',
    'zPACC': 'PACC_mri',
    'AMYLOID_STATUS': 'ab_status',
    'SUMMARY_SUVR_AMYLOID': 'ab_composite',
    'pT217_AB42_F': 'ptau',
    'MRI_SessionDate': 'mri_date',
}, inplace=True)

# Rename index column to ID
df_baseline.index.name = 'ID'

# Convert to human redable codes and standardized codes 
df_baseline['sex'] = df_baseline['sex'].map({1: 'Female', 0: 'Male'})
df_baseline['e4_carrier'] = df_baseline['e4_carrier'].map({1: 'e4+', 0: 'e4-'})
df_baseline['ab_status'] = df_baseline['ab_status'].map({1: 'ab+', 0: 'ab-'})

# Columns of interest with all the other data
common_cols = ['mri_age', 'sex', 'e4_carrier', 'ab_status', 'edu', 'diagnosis', 'mri_date', 'PACC_mri', 'time_diff_pacc',
               'ab_composite', 'time_diff_ab', 'ptau', 'time_diff_ptau', 'tau_composite', 'time_diff_tau']
df_baseline = df_baseline[common_cols]

# Save as csv
df_baseline.to_csv('data/final/adni/processed/baseline.csv')