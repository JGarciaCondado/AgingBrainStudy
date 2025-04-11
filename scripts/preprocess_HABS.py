import pandas as pd
import numpy as np

# Load baseline and longitudinal data
df_baseline = pd.read_csv('data/final/habs/habs_baseline_040925.csv', 
                          index_col='SubjIDshort',
                          usecols=lambda x: x != 'Unnamed: 0')
df_long = pd.read_csv('data/final/habs/habs_long_040925.csv',
                      index_col='ID',
                      usecols=lambda x: x != 'Unnamed: 0')

# Rename index column to ID
df_baseline.index.name = 'ID'

# Keep only participants with more than two time points
df_long = df_long.sort_values(by=['ID', 'NP_SessionDate'])
df_long['nTimePoints'] = df_long.groupby('ID')['zPACC'].transform(lambda x: x.notna().sum())
df_long = df_long[df_long['nTimePoints'] >= 2]
long_ids = df_long.index.unique()
df_baseline = df_baseline[df_baseline.index.isin(long_ids)]

# Extract followup time of subjects from longitudinal data
max_followup = df_long.groupby('ID')['MonthsFromBaseline_raw'].max().to_dict()
df_baseline['follow_up_time'] = df_baseline.index.map(max_followup)/12

# Remove any participants with diagnosis of MCI because only interested in cognitively unimpaired
df_baseline = df_baseline[df_baseline['Diagnosis'] == 'CN']

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

# Normalize volumes by total intracranial volume
def normalize_mri(mri_col, icv_col):
    from sklearn.linear_model import LinearRegression
    b = LinearRegression().fit(icv_col.values.reshape(-1,1), mri_col).coef_[0]
    adj_mri = mri_col - b * (icv_col - icv_col.mean())
    return adj_mri
df_structural['bi_Hippocampus'] = normalize_mri(df_structural['bi_Hippocampus'], df_structural['ICV_vol'])
df_structural['bi_Amygdala'] = normalize_mri(df_structural['bi_Amygdala'], df_structural['ICV_vol'])

# Keep only MRI_Age and bi_ columns
df_structural = df_structural[['MRI_Age'] + [col for col in df_structural.columns if 'bi_' in col]]

# Rename MRI_Age to age and lower case column names
df_structural = df_structural.rename(columns={'MRI_Age': 'age'})
df_structural.columns = df_structural.columns.str.lower()

# Save as csv
df_structural.to_csv('data/final/habs/processed/structural_features.csv')

# Create clinical file with CN bein e4- and ab-
df_clinical = df_baseline[['e4_carrier', 'PIB_FS_DVR_Group']].copy()
df_clinical['cn'] = ((df_clinical['e4_carrier'] == 0) & (df_clinical['PIB_FS_DVR_Group'] == 'PIB-')).astype(int)
df_clinical['e4+'] = ((df_clinical['e4_carrier'] == 1) & (df_clinical['PIB_FS_DVR_Group'] == 'PIB-')).astype(int)
df_clinical['e4+ab+'] = ((df_clinical['e4_carrier'] == 1) & (df_clinical['PIB_FS_DVR_Group'] == 'PIB+')).astype(int)
df_clinical['ab+'] = ((df_clinical['e4_carrier'] == 0) & (df_clinical['PIB_FS_DVR_Group'] == 'PIB+')).astype(int)
df_clinical['unknown'] = ((df_clinical['e4_carrier'].isna() | df_clinical['PIB_FS_DVR_Group'].isna())).astype(int)
df_clinical.drop(columns=['e4_carrier', 'PIB_FS_DVR_Group'], inplace=True)
df_clinical.to_csv('data/final/habs/processed/clinical.csv')

# Process baseline so that across all cohorts the same type of naming and conventions

# Date of each point of measurment. They were already pulled to be closest in time to first known MRI
# Date of blood samples is same as the NP date
df_baseline['MRI_SessionDate'] = pd.to_datetime(df_baseline['MRI_SessionDate'], errors='coerce')
df_baseline['NP_SessionDate'] = pd.to_datetime(df_baseline['NP_SessionDate'], errors='coerce')
df_baseline['PIB_SessionDate'] = pd.to_datetime(df_baseline['PIB_SessionDate'], errors='coerce')
df_baseline['TAU_SessionDate'] = pd.to_datetime(df_baseline['TAU_SessionDate'], errors='coerce')

# Time differences
df_baseline['time_diff_pacc'] = (df_baseline['NP_SessionDate'] - df_baseline['MRI_SessionDate']).dt.days/365.25
df_baseline['time_diff_ab'] = (df_baseline['PIB_SessionDate'] - df_baseline['MRI_SessionDate']).dt.days/365.25
df_baseline['time_diff_ptau'] = (df_baseline['NP_SessionDate'] - df_baseline['MRI_SessionDate']).dt.days/365.25
df_baseline['time_diff_tau'] = (df_baseline['TAU_SessionDate'] - df_baseline['MRI_SessionDate']).dt.days/365.25

# Calculate Tau composite by averaging over regions with _bh
# Starts with TAU_HRC_FS_SUVR_PVC_ and finishes with _bh
tau_features = [col for col in df_baseline.columns if 'TAU_HRC_FS_SUVR_PVC_' in col and '_bh' in col]
df_baseline['tau_composite'] = df_baseline[tau_features].mean(axis=1)

# Rename of columns of interest
# Ptau use the p/np ratio and the C2N platform
df_baseline.rename(columns={
    'MRI_Age': 'mri_age',
    'Sex': 'sex',
    'Education': 'edu',
    'Diagnosis': 'diagnosis',
    'zPACC': 'PACC_mri',
    'PIB_FS_DVR_Group': 'ab_status',
    'PIB_FS_DVR_FLR': 'ab_composite',
    'p_tau217_ratio': 'ptau',
    'MRI_SessionDate': 'mri_date',
}, inplace=True)


# Convert to human redable codes and standardized codes 
df_baseline['sex'] = df_baseline['sex'].map({1: 'Female', 0: 'Male'})
df_baseline['e4_carrier'] = df_baseline['e4_carrier'].map({1: 'e4+', 0: 'e4-'})
df_baseline['ab_status'] = df_baseline['ab_status'].map({'PIB+': 'ab+', 'PIB-': 'ab-'})

# Create three new columns for type of analysis, secondary and exploratory
# Primary includes those who have MRI measures and longitudinal PACC which we already filtered for
# Secondary is those who have ab_composite AND ptau measures
df_baseline['secondary'] = (df_baseline['ab_composite'].notna()).astype(int)
# Exploratory is those who have ab_composite AND ptau AND tau_composite measures
df_baseline['exploratory'] = ((df_baseline['ab_composite'].notna()) & (df_baseline['ptau'].notna()) & (df_baseline['tau_composite'].notna())).astype(int)

# Columns of interest with all the other data
common_cols = ['mri_age', 'sex', 'e4_carrier', 'ab_status', 'edu', 'diagnosis', 'mri_date', 'PACC_mri', 'time_diff_pacc', 'follow_up_time',
               'ab_composite', 'time_diff_ab', 'ptau', 'time_diff_ptau', 'tau_composite', 'time_diff_tau', 'secondary', 'exploratory']
df_baseline = df_baseline[common_cols]

# Save as csv
df_baseline.to_csv('data/final/habs/processed/baseline.csv')