import pandas as pd

# Load subject info
df_subinfo = pd.read_csv('data/SUBJINFO.csv', usecols=['BID', 'AGEYR', 'SEX', 'APOEGN'])
df_subinfo = df_subinfo.rename(columns={'AGEYR': 'AGE'})

# Load freeSurfer data
df_freesurfer = pd.read_csv('data/freesurfer/aseg_stats.txt', sep='\t')
df_freesurfer = df_freesurfer.rename(columns={'Measure:volume': 'BID'})

# Keep only specific volume columns and divide by total intracranial volume
vols = ['Left-Lateral-Ventricle', 'Left-Inf-Lat-Vent',
       'Left-Cerebellum-White-Matter', 'Left-Cerebellum-Cortex',
       'Left-Thalamus', 'Left-Caudate', 'Left-Putamen', 'Left-Pallidum',
       '3rd-Ventricle', '4th-Ventricle', 'Brain-Stem', 'Left-Hippocampus',
       'Left-Amygdala', 'CSF', 'Left-Accumbens-area', 'Left-VentralDC',
       'Left-vessel', 'Left-choroid-plexus', 'Right-Lateral-Ventricle',
       'Right-Inf-Lat-Vent', 'Right-Cerebellum-White-Matter',
       'Right-Cerebellum-Cortex', 'Right-Thalamus', 'Right-Caudate',
       'Right-Putamen', 'Right-Pallidum', 'Right-Hippocampus',
       'Right-Amygdala', 'Right-Accumbens-area', 'Right-VentralDC',
       'Right-vessel', 'Right-choroid-plexus', 'CC_Posterior', 'CC_Mid_Posterior', 'CC_Central',
       'CC_Mid_Anterior', 'CC_Anterior', 'lhCortexVol', 'rhCortexVol', 'CortexVol',
       'lhCerebralWhiteMatterVol', 'rhCerebralWhiteMatterVol', 'CerebralWhiteMatterVol',
       'SubCortGrayVol', 'TotalGrayVol', 'SupraTentorialVol']
df_freesurfer = df_freesurfer[vols + ['BID']]

# Remove outliers
outliers = ['B83014615', 'B16568072', 'B11279364', 'B65897427', 'B19132670', 'B49144285', 'B12659422', 'B65278081', 'B87879018']
df_freesurfer = df_freesurfer[~df_freesurfer['BID'].isin(outliers)]


# Merge subject info and freeSurfer data to create features
df_features = pd.merge(df_subinfo[['BID', 'AGE']], df_freesurfer, on='BID')
id_features = set(df_features['BID'].to_list())

# Save without index
df_features.to_csv('data/ageml/features.csv', index=False)

# Save covariates
df_covariates = df_subinfo[['BID', 'SEX']]
df_covariates = df_covariates[df_covariates['BID'].isin(id_features)]
df_covariates.to_csv('data/ageml/covariates.csv', index=False)

# Create clinical files for gender, one column is Female = 1 and Male = 2
df_sex = df_subinfo[['BID', 'SEX']]
df_sex = df_sex[df_sex['BID'].isin(id_features)]
df_sex['Female'] = (df_sex['SEX'] == 1).astype(int)
df_sex['Male'] = (df_sex['SEX'] == 2).astype(int)
df_sex.drop(columns='SEX', inplace=True)
# Rename Male to Control
df_sex.rename(columns={'Female': 'CN'}, inplace=True)
df_sex.to_csv('data/ageml/clinical/sex.csv', index=False)

# Obtain amyloid status and ApoGen
df_status = pd.read_csv('data/pet_imaging/imaging_PET_VA.csv', usecols=['BID', 'overall_score'])
df_status.rename(columns={'overall_score': 'Amyloid'}, inplace=True)
df_status = pd.merge(df_status, df_subinfo[['BID', 'APOEGN']], on='BID').dropna()
df_status = df_status[df_status['BID'].isin(id_features)]
df_status['Amyloid'] = (df_status['Amyloid'] == 'positive').astype(int)
df_status['e4'] = df_status['APOEGN'].str.contains('E4').astype(int)

# Create clinical file with type of Amyloid and ApoE4
df_clinical = df_status.copy()
df_clinical['CN'] = ((df_clinical['Amyloid'] == 0) & (df_clinical['e4'] == 0)).astype(int)
df_clinical['A+e4-'] = ((df_clinical['Amyloid'] == 1) & (df_clinical['e4'] == 0)).astype(int)
df_clinical['A-e4+'] = ((df_clinical['Amyloid'] == 0) & (df_clinical['e4'] == 1)).astype(int)
df_clinical['A+e4+'] = ((df_clinical['Amyloid'] == 1) & (df_clinical['e4'] == 1)).astype(int)
df_clinical.drop(columns=['Amyloid', 'e4', 'APOEGN'], inplace=True)
df_clinical.to_csv('data/ageml/clinical/a_e4.csv', index=False)

# Create clinical file with type of Apoe4
df_clinical = df_status.copy()
df_clinical['CN'] = ((df_clinical['e4'] == 0)).astype(int)
df_clinical.drop(columns=['Amyloid', 'APOEGN'], inplace=True)
df_clinical.to_csv('data/ageml/clinical/e4.csv', index=False)

# Load plasma biomarkers (they are saved seperaetly because of NaN issues with AgeML)

# Obrain pTAU data make BID coulmn index and only extract ORRES column
df_ptau = pd.read_csv('data/plasma/biomarker_pTau217.csv', usecols=['BID', 'VISCODE', 'ORRES'])
df_ptau.rename(columns={'ORRES': 'pTau217'}, inplace=True)
# Keep only those with VISCODE 6
df_ptau = df_ptau[df_ptau['VISCODE'] == 6]
df_ptau.drop(columns=['VISCODE'], inplace=True)
# Convert ptau to float and remove nans
df_ptau['pTau217'] = pd.to_numeric(df_ptau['pTau217'], errors='coerce')
df_ptau.dropna(inplace=True)
# Save ptau
df_ptau.to_csv('data/ageml/factors/pTau217.csv', index=False)

# Obtain AB Test data 
df_ab = pd.read_csv('data/plasma/biomarker_AB_Test.csv', usecols=['BID', 'LBTESTCD', 'LBORRES'])
# Pivot to wide dataframe and remove duplicates
df_ab = df_ab.drop_duplicates(subset=['BID', 'LBTESTCD'], keep='first')
df_ab = df_ab.pivot(index='BID', columns='LBTESTCD', values='LBORRES')
# Convert to floats
df_ab = df_ab.apply(pd.to_numeric, errors='coerce')
# Save AB test data
df_ab.to_csv('data/ageml/factors/AB_test.csv')

# Obtain Plasma Roche measures
df_roche = pd.read_csv('data/plasma/biomarker_Plasma_Roche_Results.csv', usecols=['BID', 'LBTESTCD', 'LABRESN'])

# Pivot to wide dataframe and remove duplicates
df_roche = df_roche.drop_duplicates(subset=['BID', 'LBTESTCD'], keep='first')
df_roche = df_roche.pivot(index='BID', columns='LBTESTCD', values='LABRESN')
# Convert to floats
df_roche = df_roche.apply(pd.to_numeric, errors='coerce')
# Remove apoe4 column
df_roche.drop(columns=['APOE4'], inplace=True)
# Save Roche data
df_roche.to_csv('data/ageml/factors/roche.csv')

# Load Cognitive data
df_cog = pd.read_csv('data/cognition/PACC.csv', usecols=['BID', 'VISCODE', 'PACC.raw', "FCTOTAL96.z","LDELTOTAL.z","DIGITTOTAL.z","MMSCORE.z"])
df_cog.rename(columns={'PACC.raw': 'PACC', 'FCTOTAL96.z': 'FCTOTAL96', 'LDELTOTAL.z': 'LDELTOTAL', 'DIGITTOTAL.z': 'DIGITTOTAL', 'MMSCORE.z': 'MMSCORE'}, inplace=True)
df_cog = df_cog[df_cog['BID'].isin(id_features)]
df_cog = df_cog[df_cog['VISCODE'] == 1]
df_cog.drop(columns=['VISCODE'], inplace=True)
df_cog.to_csv('data/ageml/factors/cognition.csv', index=False)