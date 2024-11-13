import pandas as pd

# Load subject info
df_subinfo = pd.read_csv('data/SUBJINFO.csv', usecols=['BID', 'AGEYR'])
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
       'Right-vessel', 'Right-choroid-plexus', '5th-Ventricle',
       'Optic-Chiasm', 'CC_Posterior', 'CC_Mid_Posterior', 'CC_Central',
       'CC_Mid_Anterior', 'CC_Anterior', 'lhCortexVol', 'rhCortexVol', 'CortexVol',
       'lhCerebralWhiteMatterVol', 'rhCerebralWhiteMatterVol', 'CerebralWhiteMatterVol',
       'SubCortGrayVol', 'TotalGrayVol', 'SupraTentorialVol']
df_freesurfer = df_freesurfer[vols + ['BID']]

# Find outliers if there are more than 3std away from the mean in total grey volume
df_freesurfer['TotalGrayVol'] = df_freesurfer['TotalGrayVol'] / df_freesurfer['SupraTentorialVol']



# Merge subject info and freeSurfer data to create features
df_features = pd.merge(df_subinfo, df_freesurfer, on='BID')

# Save without index
df_features.to_csv('data/features.csv', index=False)

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
df_ptau.to_csv('data/plasma/pTau217.csv', index=False)

# Obtain AB Test data 
df_ab = pd.read_csv('data/plasma/biomarker_AB_Test.csv', usecols=['BID', 'LBTESTCD', 'LBORRES'])
# Pivot to wide dataframe and remove duplicates
df_ab = df_ab.drop_duplicates(subset=['BID', 'LBTESTCD'], keep='first')
df_ab = df_ab.pivot(index='BID', columns='LBTESTCD', values='LBORRES')
# Convert to floats
df_ab = df_ab.apply(pd.to_numeric, errors='coerce')
# Save AB test data
df_ab.to_csv('data/plasma/AB_test.csv')

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
df_roche.to_csv('data/plasma/roche.csv')