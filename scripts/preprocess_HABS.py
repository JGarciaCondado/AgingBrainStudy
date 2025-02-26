import pandas as pd
import numpy as np

# Load features but ignore first column
df_structural = pd.read_csv('data/HABS/raw/raw_features.csv', index_col='SubjIDshort')

# Remove MRI_FS6_ADNI_ from column names
df_structural.columns = df_structural.columns.str.replace('MRI_FS6_ADNI_', '')

# Combine those with lh and rh 
lh = df_structural.filter(like='lh_').columns
rh = df_structural.filter(like='rh_').columns
for l, r in zip(lh, rh):
    df_structural['bi_' + l.split('_')[1]] = df_structural[l] + df_structural[r]
    df_structural.drop(columns=[l, r], inplace=True)

# Combine those with Left and Right
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
df_structural.to_csv('data/HABS/processed/structural_features.csv')

# Load PET data
df_pet = pd.read_csv('data/HABS/raw/raw_pet_features.csv', index_col='SubjIDshort')
df_pet.columns = df_pet.columns.str.replace('TAU_FS_SUVR_PVC_', '')
df_pet = df_pet.rename(columns={'TAU_Age': 'age'})
df_pet.columns = df_pet.columns.str.lower()
df_pet = df_pet[['age'] + [col for col in df_pet.columns if '_bh' in col]]

# Remove outlier
outliers = ['B_RAUOBJ']

# Save as csv
df_pet.to_csv('data/HABS/processed/pet_features.csv')

# Load factors
df_factors = pd.read_csv('data/HABS/raw/raw_factors.csv', index_col='SubjIDshort')

# Create a FCTOTAL96 by adding all columns that start with NP_FCsrt_
df_factors['FCTOTAL96'] = df_factors.filter(like='NP_FCsrt_').sum(axis=1)
df_factors.drop(columns=df_factors.filter(like='NP_FCsrt_').columns, inplace=True)
df_factors = df_factors[df_factors['FCTOTAL96'] != 0]

# Create Tau Composite
tau_features = [col for col in df_pet.columns if '_bh' in col]
df_pet['tau_composite'] = df_pet[tau_features].mean(axis=1)

# Load blood biomarkers MSD
df_blood = pd.read_csv('data/HABS/raw/raw_blood_biomarkers.csv')
# Keep first sample
df_blood['StudyArc'] = df_blood['StudyArc'].str.replace('HAB_', '').str.replace('.0', '').astype(int)
df_blood = df_blood.sort_values('StudyArc').drop_duplicates(subset='SubjIDshort', keep='first')
# Format
df_blood.set_index('SubjIDshort', inplace=True)
df_blood = df_blood.filter(like='_conc')
df_blood.columns = df_blood.columns.str.replace('_conc', '')
df_blood.drop(columns=['ptau', 'tota'], inplace=True) # prefer blood biomarkers form C2N

# Load ptau data from C2N
df_ptau = pd.read_csv('data/HABS/raw/raw_pTau.csv', usecols=['SubjIDshort', 'StudyArc', 'p_tau217', 'p_tau217_ratio'])
df_ptau['StudyArc'] = df_ptau['StudyArc'].str.replace('HAB_', '').str.replace('.0', '').astype(int)
df_ptau = df_ptau.sort_values('StudyArc').drop_duplicates(subset='SubjIDshort', keep='first')
df_ptau.drop(columns='StudyArc', inplace=True)

# Merge pet and blood to factors
df_factors = pd.merge(df_factors, df_pet[tau_features + ['tau_composite']], on='SubjIDshort', how='outer')
df_factors = pd.merge(df_factors, df_blood, on='SubjIDshort', how='outer')
df_factors = pd.merge(df_factors, df_ptau, on='SubjIDshort', how='outer')

# Save as csv
df_factors.to_csv('data/HABS/processed/factors.csv')

# Create clinical file with CN where 1 is E4_status is e4-
df_clinical = df_factors[['SubjIDshort','E4_Status', 'PIB_FS_SUVR_Group']].copy()
df_clinical = df_clinical.dropna()
df_clinical['cn'] = ((df_clinical['E4_Status'] == 'e4-') & (df_clinical['PIB_FS_SUVR_Group'] == 'PIB-')).astype(int)
df_clinical['e4+'] = ((df_clinical['E4_Status'] == 'e4+') & (df_clinical['PIB_FS_SUVR_Group'] == 'PIB-')).astype(int)
df_clinical['e4+ab+'] = ((df_clinical['E4_Status'] == 'e4+') & (df_clinical['PIB_FS_SUVR_Group'] == 'PIB+')).astype(int)
df_clinical['ab+'] = ((df_clinical['E4_Status'] == 'e4-') & (df_clinical['PIB_FS_SUVR_Group'] == 'PIB+')).astype(int)
df_clinical.drop(columns=['E4_Status', 'PIB_FS_SUVR_Group'], inplace=True)
df_clinical.to_csv('data/HABS/processed/clinical.csv', index=False)

# Joint information on A4 and HABS

# Load A4 structural data 
df_a4_structural = pd.read_csv('data/A4/processed/structural_features.csv', index_col='BID')
df_a4_structural.columns = df_a4_structural.columns.str.lower()

# Stack dataframes
df_structural = pd.concat([df_structural, df_a4_structural], axis=0)
df_structural.index.name = 'SubjIDshort'
df_structural.to_csv('data/A4_HABS_joint/structural_features.csv')

# Create new dataframe with CN column with 1 for those in a4 and 0 for those in HABS
df_structural['cn'] = df_structural.index.isin(df_a4_structural.index).astype(int)
clinical = df_structural[['cn']].copy()
clinical['HABS'] = (~clinical['cn'].astype(bool)).astype(int)
df_structural.drop(columns='cn', inplace=True)
clinical.to_csv('data/A4_HABS_joint/A4_as_cn.csv')

# Load A4 PET data
df_a4_pet = pd.read_csv('data/A4/processed/pet_features.csv', index_col='BID')
df_a4_pet.columns = df_a4_pet.columns.str.lower()
df_a4_pet.columns = ['_'.join(col.split('_')[1:]) + '_bh' if 'bi_' in col else col for col in df_a4_pet.columns]

# Stack dataframes
df_pet.drop(columns='tau_composite', inplace=True)
df_pet = pd.concat([df_pet, df_a4_pet], axis=0)
df_pet.index.name = 'SubjIDshort'
df_pet.to_csv('data/A4_HABS_joint/pet_features.csv')