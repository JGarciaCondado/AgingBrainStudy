import pandas as pd
import numpy as np

# Load baseline and longitudinal data
df_baseline = pd.read_csv('data/final/a4/a4_baseline_041125.csv', 
                          index_col='ID',
                          usecols=lambda x: x != 'Unnamed: 0')
df_long = pd.read_csv('data/final/a4/a4_long_041125.csv',
                      index_col='ID',
                      usecols=lambda x: x != 'Unnamed: 0')

# Remove participant with outlier MRI B66388909_a4
df_baseline = df_baseline[df_baseline.index != 'B66388909_a4']

# Keep only participants with more than two time points
df_long = df_long.sort_values(by=['ID', 'MonthsFromBaseline_raw'])
df_long['nTimePoints'] = df_long.groupby('ID')['zPACC'].transform(lambda x: x.notna().sum())
df_long = df_long[df_long['nTimePoints'] >= 2]
long_ids = df_long.index.unique()
df_baseline = df_baseline[df_baseline.index.isin(long_ids)]

# Extract followup time of subjects from longitudinal data
max_followup = df_long.groupby('ID')['MonthsFromBaseline_raw'].max().to_dict()
df_baseline['follow_up_time'] = df_baseline.index.map(max_followup)/12

# Only keep subjects that have MRI and Amyloid measures
df_baseline = df_baseline[(df_baseline['summary_suvr_amyloid'].notna()) & (df_baseline['MRI_Age'].notna())]

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
df_structural.to_csv('data/final/a4/processed/structural_features.csv')

# Create clinical file with CN being e4- and ab-
df_clinical = df_baseline[['e4_carrier', 'Amyloid_group']].copy()
df_clinical['cn'] = ((df_clinical['e4_carrier'] == 0) & (df_clinical['Amyloid_group'] == 'Ab-')).astype(int)
df_clinical['not_cn'] = 1 - df_clinical['cn'].astype(int)
df_clinical.drop(columns=['e4_carrier', 'Amyloid_group'], inplace=True)
df_clinical.to_csv('data/final/a4/processed/clinical.csv')

# A4 has no dates so we have to use months from baseline for each
df_baseline['time_diff_pacc'] = (df_baseline['MonthsFromBaseline_raw'] - df_baseline['MonthsFromBaseline_MRI'])/12
df_baseline['time_diff_ab'] = (df_baseline['MonthsFromBaseline_amyloid'] - df_baseline['MonthsFromBaseline_MRI'])/12
df_baseline['time_diff_ptau'] = (df_baseline['MonthsFromBaseline_ptau217'] - df_baseline['MonthsFromBaseline_MRI'])/12
df_baseline['time_diff_tau'] = (df_baseline['MonthsFromBaseline_tau'] - df_baseline['MonthsFromBaseline_MRI'])/12

# Calculate Tau composite by averaging over regions with _bh
# Starts with PVC and finishes with _bh
tau_features = [col for col in df_baseline.columns if 'PVC_' in col and '_bh' in col]
df_baseline['tau_composite'] = df_baseline[tau_features].mean(axis=1)

# Rename of columns of interest
# Ptau use the p/np ratio and the C2N platform
df_baseline.rename(columns={
    'MRI_Age': 'mri_age',
    'Sex': 'sex',
    'Education': 'edu',
    'zPACC': 'PACC_mri',
    'Amyloid_group': 'ab_status',
    'Amyloid_Centiloid': 'ab_composite',
    'ptau217_read': 'ptau',
    'MonthsFromBaseline_MRI': 'time_baseline_to_mri',
}, inplace=True)

# Convert to human redable codes and standardized codes 
df_baseline['sex'] = df_baseline['sex'].map({1: 'Female', 0: 'Male'})
df_baseline['e4_carrier'] = df_baseline['e4_carrier'].map({1: 'e4+', 0: 'e4-'})
df_baseline['ab_status'] = df_baseline['ab_status'].map({'Ab+': 'ab+', 'Ab-': 'ab-'})

# All are technically CU at start of the study. Create new column with all CN
df_baseline['diagnosis'] = 'CN'

# Create cohort variable
conditions = [
    ((df_baseline['SUBSTUDY'] == 'SF') | (df_baseline['SUBSTUDY'] == 'LEARN')),
    ((df_baseline['SUBSTUDY'] == 'A4') & (df_baseline['TX'] == 'Placebo')),
    ((df_baseline['SUBSTUDY'] == 'A4') & (df_baseline['TX'] == 'Solanezumab'))
]

choices = ['LEARN/SF', 'A4 Placebo', 'A4 Treated']

# Apply the conditions
df_baseline['cohort'] = np.select(conditions, choices, default=None)

# Change ptau217 reads of <LLOQ and >ULOQ to NaN
df_baseline['ptau'] = df_baseline['ptau'].replace('<LLOQ', np.nan)
df_baseline['ptau'] = df_baseline['ptau'].replace('>ULOQ', np.nan)

# Create one new column for type of analysis: exploratory
# Primary includes those who have MRI measures, longitudinal PACC  and AB measures which we already filtered for
# Exploratory is those who also ptau AND tau_composite measures
df_baseline['exploratory'] = ((df_baseline['ptau'].notna()) & (df_baseline['tau_composite'].notna())).astype(int)

# Columns of interest with all the other data
common_cols = ['mri_age', 'sex', 'e4_carrier', 'ab_status', 'edu', 'diagnosis', 'cohort', 'TX', 'time_baseline_to_mri', 'PACC_mri', 'time_diff_pacc', 'follow_up_time',
               'ab_composite', 'time_diff_ab', 'ptau', 'time_diff_ptau', 'tau_composite', 'time_diff_tau', 'exploratory']
df_baseline = df_baseline[common_cols]

# Save as csv
df_baseline.to_csv('data/final/a4/processed/baseline.csv')