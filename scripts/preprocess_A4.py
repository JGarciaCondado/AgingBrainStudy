import pandas as pd

# Load subject info
df_subinfo = pd.read_csv('data/A4/raw/SUBJINFO.csv', usecols=['BID', 'AGEYR', 'SEX', 'APOEGN'])
df_subinfo = df_subinfo.rename(columns={'AGEYR': 'AGE'})

# Load freeSurfer volume data
df_freesurfer = pd.read_csv('data/A4/raw/freesurfer/aseg_stats.txt', sep='\t')
df_freesurfer = df_freesurfer.rename(columns={'Measure:volume': 'BID'})

# Keep only specific volume columns and divide by total intracranial volume
vols = [ 'Left-Hippocampus', 'Left-Amygdala', 'Right-Hippocampus', 'Right-Amygdala', 'EstimatedTotalIntraCranialVol']
df_freesurfer = df_freesurfer[['BID'] + vols]
# Combine left and right volumes
df_freesurfer['bi_hippocampus'] = df_freesurfer['Left-Hippocampus'] + df_freesurfer['Right-Hippocampus']
df_freesurfer['bi_amygdala'] = df_freesurfer['Left-Amygdala'] + df_freesurfer['Right-Amygdala']
def normalize_mri(mri_col, icv_col):
    from sklearn.linear_model import LinearRegression
    b = LinearRegression().fit(icv_col.values.reshape(-1,1), mri_col).coef_[0]
    adj_mri = mri_col - b * (icv_col - icv_col.mean())
    return adj_mri
# Normalize volumes by total intracranial volume
df_freesurfer['bi_hippocampus'] = normalize_mri(df_freesurfer['bi_hippocampus'], df_freesurfer['EstimatedTotalIntraCranialVol'])
df_freesurfer['bi_amygdala'] = normalize_mri(df_freesurfer['bi_amygdala'], df_freesurfer['EstimatedTotalIntraCranialVol'])
df_freesurfer.drop(columns=vols, inplace=True)

# Load thickness right and left
df_thickness_lh = pd.read_csv('data/A4/raw/freesurfer/thickness_lh_stats.txt', sep='\t').rename(columns={'lh.aparc.thickness': 'BID'})
df_thickness_rh = pd.read_csv('data/A4/raw/freesurfer/thickness_rh_stats.txt', sep='\t').rename(columns={'rh.aparc.thickness': 'BID'})
df_thickness = pd.merge(df_thickness_lh, df_thickness_rh, on='BID', suffixes=('_lh', '_rh'))
# Keep only specific thickness columns
thickness = ['inferiortemporal', 'inferiorparietal', 'fusiform', 'middletemporal', 'entorhinal', 'parahippocampal']
# Contains some of this words
df_thickness = df_thickness[['BID'] + [col for col in df_thickness.columns if any([t in col for t in thickness])]]
# Combine left and right thickness
for t in thickness:
    df_thickness['bi_' + t] = df_thickness['lh_'+t+'_thickness'] + df_thickness['rh_'+t+'_thickness']
    df_thickness.drop(columns=['lh_'+t+'_thickness', 'rh_'+t+'_thickness'], inplace=True)
df_freesurfer = pd.merge(df_freesurfer, df_thickness, on='BID')

# Remove outliers
outliers = ['B83014615', 'B16568072', 'B11279364', 'B65897427', 'B19132670', 'B49144285', 'B12659422', 'B65278081', 'B87879018']
df_freesurfer = df_freesurfer[~df_freesurfer['BID'].isin(outliers)]

# Merge subject info and freeSurfer data to create features
df_features = pd.merge(df_subinfo[['BID', 'AGE']], df_freesurfer, on='BID')

# Save without index
df_features.to_csv('data/A4/processed/structural_features.csv', index=False)

# Save covariates
df_covars = df_subinfo[['BID', 'SEX', 'APOEGN']]

# Obtain amyloid status and ApoGen
df_status = pd.read_csv('data/A4/raw/pet_imaging/imaging_PET_VA.csv', usecols=['BID', 'overall_score'])
df_status.rename(columns={'overall_score': 'Amyloid'}, inplace=True)
df_covars = pd.merge(df_covars, df_status, on='BID').dropna()
df_covars['Amyloid'] = (df_status['Amyloid'] == 'positive').astype(int)
df_covars['e4'] = df_covars['APOEGN'].str.contains('E4').astype(int)
df_covars.drop(columns=['APOEGN'], inplace=True)
df_covars.to_csv('data/A4/processed/covariates.csv', index=False)

# Create clinical file with type of Apoe4
df_clinical = df_covars[['BID', 'e4']].copy()
df_clinical['CN'] = ((df_clinical['e4'] == 0)).astype(int)
df_clinical.to_csv('data/A4/processed/e4.csv', index=False)

# Obrain pTAU data make BID coulmn index and only extract ORRES column
df_ptau = pd.read_csv('data/A4/raw/plasma/biomarker_pTau217.csv', usecols=['BID', 'VISCODE', 'ORRES'])
df_ptau.rename(columns={'ORRES': 'pTau217'}, inplace=True)
# Keep only those with VISCODE 6
df_ptau = df_ptau[df_ptau['VISCODE'] == 6]
df_ptau.drop(columns=['VISCODE'], inplace=True)
# Convert ptau to float and remove nans
df_ptau['pTau217'] = pd.to_numeric(df_ptau['pTau217'], errors='coerce')
df_ptau.dropna(inplace=True)

# Obtain AB Test data 
df_ab = pd.read_csv('data/A4/raw/plasma/biomarker_AB_Test.csv', usecols=['BID', 'LBTESTCD', 'LBORRES'])
# Pivot to wide dataframe and remove duplicates
df_ab = df_ab.drop_duplicates(subset=['BID', 'LBTESTCD'], keep='first')
df_ab = df_ab.pivot(index='BID', columns='LBTESTCD', values='LBORRES')
# Convert to floats
df_ab = df_ab.apply(pd.to_numeric, errors='coerce')
# Only keep TP42/TP40 and FP42/FP40
df_ab = df_ab[['TP42/TP40', 'FP42/FP40']]

# Obtain Plasma Roche measures
df_roche = pd.read_csv('data/A4/raw/plasma/biomarker_Plasma_Roche_Results.csv', usecols=['BID', 'LBTESTCD', 'LABRESN'])

# Pivot to wide dataframe and remove duplicates
df_roche = df_roche.drop_duplicates(subset=['BID', 'LBTESTCD'], keep='first')
df_roche = df_roche.pivot(index='BID', columns='LBTESTCD', values='LABRESN')
df_roche = df_roche.apply(pd.to_numeric, errors='coerce')
# Calculate the ratio of AB42/AB40
df_roche['AB42/AB40'] = df_roche['AMYLB42'] / df_roche['AMYLB40']
df_roche.drop(columns=['AMYLB42', 'AMYLB40'], inplace=True)

# Load Cognitive data
df_cog = pd.read_csv('data/A4/raw/cognition/PACC.csv', usecols=['BID', 'VISCODE', 'PACC.raw', "FCTOTAL96.z","LDELTOTAL.z","DIGITTOTAL.z","MMSCORE.z"])
df_cog.rename(columns={'PACC.raw': 'PACC', 'FCTOTAL96.z': 'FCTOTAL96', 'LDELTOTAL.z': 'LDELTOTAL', 'DIGITTOTAL.z': 'DIGITTOTAL', 'MMSCORE.z': 'MMSCORE'}, inplace=True)
df_cog = df_cog[df_cog['VISCODE'] == 1]
df_cog.drop(columns=['VISCODE'], inplace=True)

# Load PET data Amyloid
df_pet_amyl = pd.read_csv('data/A4/raw/pet_imaging/imaging_SUVR_amyloid.csv', usecols=['BID', 'VISCODE', 'brain_region', 'suvr_cer'])
df_pet_amyl = df_pet_amyl[df_pet_amyl['VISCODE'] == 2]
df_pet_amyl = df_pet_amyl.drop_duplicates(subset=['BID', 'brain_region'], keep='first')
df_pet_amyl = df_pet_amyl.dropna(subset=['suvr_cer'])
df_pet_amyl = df_pet_amyl.drop_duplicates(subset=['BID', 'brain_region'], keep='first')
df_pet_amyl = df_pet_amyl.pivot(index='BID', columns='brain_region', values='suvr_cer')
df_pet_amyl = df_pet_amyl.apply(pd.to_numeric, errors='coerce')
df_pet_amyl.rename(columns={'Composite_Summary': 'amyloid_composite'}, inplace=True)

# Load PET Tau PETSurfer
df_pet_tau = pd.read_csv('data/A4/raw/pet_imaging/imaging_Tau_PET_PetSurfer.csv')
df_pet_tau = df_pet_tau[['BID'] + [col for col in df_pet_tau.columns if 'bi_' in col]]
# Only interested in specifci features
tau_features = ['bi_inferiortemporal', 'bi_inferiorparietal', 'bi_fusiform', 'bi_middletemporal', 'bi_entorhinal', 'bi_Amygdala', 'bi_parahippocampal']
df_pet_tau = df_pet_tau[['BID'] + tau_features]
# Create Tau composite
df_pet_tau['tau_composite'] = df_pet_tau[tau_features].mean(axis=1)

# Merge all into one csv
df_factors = pd.merge(df_ptau, df_ab, on='BID', how='outer')
df_factors = pd.merge(df_factors, df_roche, on='BID', how='outer')
df_factors = pd.merge(df_factors, df_cog, on='BID', how='outer')
df_factors = pd.merge(df_factors, df_pet_amyl, on='BID', how='outer')
df_factors = pd.merge(df_factors, df_pet_tau, on='BID', how='outer')
df_factors.to_csv('data/A4/processed/factors.csv', index=False)

# Add Age to create age model
df_age = df_subinfo[['BID', 'AGE']]
df_pet_tau = pd.merge(df_age, df_pet_tau, on='BID')
df_pet_tau.drop(columns='tau_composite', inplace=True)
df_pet_tau.to_csv('data/A4/processed/pet_features.csv', index=False)

# Define PRS type
prs_types = ['gm', 'wm', 'fc']
prs_threshold = 0.5 # 0.001, 0.05 0.1, 0.2, 0.3, 0.4, 0.5

# Create PRS data dataframe
for prs_type in prs_types:
    file_path = 'data/A4/raw/prs_scores/{}/prs2.pT{}.sscore'.format(prs_type, prs_threshold)
    df_prs = pd.read_csv(file_path, sep='\s+')
    df_prs = df_prs[['IID', 'SCORE1_AVG']]
    df_prs.rename(columns={'SCORE1_AVG': 'SCORE'}, inplace=True)

    # Convert SCORE to z-score
    df_prs['ZSCORE'] = (df_prs['SCORE'] - df_prs['SCORE'].mean()) / df_prs['SCORE'].std()

    # Rename columns
    df_prs.rename(columns={'SCORE': '{}_SCORE'.format(prs_type.upper()), 'ZSCORE': '{}_ZSCORE'.format(prs_type.upper())}, inplace=True)

    # Merge dataframes
    if prs_type == 'gm':
        df = df_prs
    else:
        df = pd.merge(df, df_prs, on='IID', how='outer')

# Dorp score columns
df.drop(columns=['GM_SCORE', 'WM_SCORE', 'FC_SCORE'], inplace=True)
df.rename(columns={'IID': 'BID'}, inplace=True)
df.to_csv(f'data/A4/processed/prs_{prs_threshold}.csv', index=False)