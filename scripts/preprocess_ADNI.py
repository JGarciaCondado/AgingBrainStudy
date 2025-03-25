import pandas as pd

# Load the dataset
df = pd.read_csv('data/ADNI/raw/adni_03242025.csv', index_col=0)

# Keep only baseline data
df_baseline = df[df["MonthsFromBaseline_PACC"] == 0.0].copy()

# Remove those with NaNs in MRI
df_baseline.dropna(subset=['MRI_FS7_rnr_ICV_vol'], inplace=True)


# Preprocess thickness features cotains thk
features_thickness = [col for col in df_baseline.columns if 'thk' in col]
features_thickness = [col.split(sep = '_')[4] for col in features_thickness if 'lh' in col]

# Merge left and right
for feature in features_thickness:
    df_baseline['bi_' + feature] = df_baseline['MRI_FS7_rnr_lh_' + feature + '_thk'] + df_baseline['MRI_FS7_rnr_rh_' + feature + '_thk']
    df_baseline.drop(columns=['MRI_FS7_rnr_lh_' + feature + '_thk', 'MRI_FS7_rnr_rh_' + feature + '_thk'], inplace=True)

# Find volumetric features
features_vol = ['Hippocampus', 'Amygdala']
for feature in features_vol:
    df_baseline['bi_' + feature.lower()] = df_baseline['MRI_FS7_rnr_Left_' + feature + '_vol'] + df_baseline['MRI_FS7_rnr_Right_' + feature + '_vol']
    df_baseline.drop(columns=['MRI_FS7_rnr_Left_' + feature + '_vol', 'MRI_FS7_rnr_Right_' + feature + '_vol'], inplace=True)

# Normalize by total intracranial volume
def normalize_mri(mri_col, icv_col, mask):
    from sklearn.linear_model import LinearRegression
    model = LinearRegression().fit(icv_col[mask].values.reshape(-1, 1), mri_col[mask])
    b = model.coef_[0]
    c = icv_col[mask].mean()
    adj_mri = mri_col - b * (icv_col - c)
    return adj_mri

mask_cn = df_baseline['Baseline_Diagnosis'] == 'CN'
df_baseline['bi_hippocampus'] = normalize_mri(df_baseline['bi_hippocampus'], df_baseline['MRI_FS7_rnr_ICV_vol'], mask_cn)
df_baseline['bi_amygdala'] = normalize_mri(df_baseline['bi_amygdala'], df_baseline['MRI_FS7_rnr_ICV_vol'], mask_cn)

# Create features file for use with AgeML
df_age = df_baseline[['ID', 'Age', 'bi_hippocampus', 'bi_amygdala'] + ['bi_' + feature for feature in features_thickness]]
df_age.to_csv('data/ADNI/processed/structural_features.csv', index=False)

# Create covariates file
covars = ['Age', 'Sex', 'Education', 'e4_carrier', 'AMYLOID_STATUS', 'Baseline_Diagnosis']
df_covars = df_baseline[['ID'] + covars]
df_covars.to_csv('data/ADNI/processed/covariates.csv', index=False)

# Create clinical file with type of Apoe4
df_clinical = df_covars[['ID', 'e4_carrier', 'AMYLOID_STATUS', 'Baseline_Diagnosis']].copy()
df_clinical['CN'] = ((df_clinical['e4_carrier'] == 0.0) & (df_clinical['AMYLOID_STATUS'] == 0.0) & (df_baseline['Baseline_Diagnosis'] == 'CN')).astype(int)
df_clinical['e4+'] = ((df_clinical['e4_carrier'] == 1.0) & (df_clinical['AMYLOID_STATUS'] == 0.0) & (df_baseline['Baseline_Diagnosis'] == 'CN')).astype(int)
df_clinical['e4+ab+'] = ((df_clinical['e4_carrier'] == 1.0) & (df_clinical['AMYLOID_STATUS'] == 1.0) & (df_baseline['Baseline_Diagnosis'] == 'CN')).astype(int)
df_clinical['ab+'] = ((df_clinical['e4_carrier'] == 0.0) & (df_clinical['AMYLOID_STATUS'] == 1.0) & (df_baseline['Baseline_Diagnosis'] == 'CN')).astype(int)
df_clinical['MCI'] = ((df_clinical['e4_carrier'] == 0.0) & (df_clinical['AMYLOID_STATUS'] == 0.0) & (df_baseline['Baseline_Diagnosis'] == 'MCI')).astype(int)
df_clinical['MCIe4+'] = ((df_clinical['e4_carrier'] == 1.0) & (df_clinical['AMYLOID_STATUS'] == 0.0) & (df_baseline['Baseline_Diagnosis'] == 'MCI')).astype(int)
df_clinical['MCIe4+ab+'] = ((df_clinical['e4_carrier'] == 1.0) & (df_clinical['AMYLOID_STATUS'] == 1.0) & (df_baseline['Baseline_Diagnosis'] == 'MCI')).astype(int)
df_clinical['MCIab+'] = ((df_clinical['e4_carrier'] == 0.0) & (df_clinical['AMYLOID_STATUS'] == 1.0) & (df_baseline['Baseline_Diagnosis'] == 'MCI')).astype(int)
df_clinical['CNunknown'] = ((df_clinical['e4_carrier'].isna() | df_clinical['AMYLOID_STATUS'].isna()) & (df_baseline['Baseline_Diagnosis'] == 'CN')).astype(int)
df_clinical['MCIunknown'] = ((df_clinical['e4_carrier'].isna() | df_clinical['AMYLOID_STATUS'].isna()) & (df_baseline['Baseline_Diagnosis'] == 'MCI')).astype(int)
df_clinical.drop(columns=['e4_carrier', 'AMYLOID_STATUS', 'Baseline_Diagnosis'], inplace=True)
df_clinical.to_csv('data/ADNI/processed/clinical.csv', index=False)

# Extract factors of interest
factors_cog = ['zPACC']
factors_blood = ['pT217_F', 'AB42_AB40_F', 'NfL_Q', 'GFAP_Q']
factors_pet_amyl = ['CENTILOIDS_AMYLOID']
factors_pet_tau = [col for col in df_baseline.columns if 'SUVR_TAU' in col]
df_factors = df_baseline[['ID'] + factors_cog + factors_blood + factors_pet_amyl + factors_pet_tau].copy()

# Calculate tau summary score
df_factors['tau_composite'] = df_factors[factors_pet_tau].mean(axis=1)

# Rename columns
df_factors.rename(columns={col: col.replace('_F', '').replace('_Q', '') for col in factors_blood}, inplace=True)
df_factors.rename(columns={col: col.replace('PVC_CTX_', '').replace('_SUVR_TAU', '').lower() + '_tau' for col in factors_pet_tau}, inplace=True)
df_factors.rename(columns={'SUMMARY_SUVR_AMYLOID': 'amyloid_composite'}, inplace=True)
df_factors.to_csv('data/ADNI/processed/factors.csv', index=False)

