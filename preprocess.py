import pandas as pd

# Define PRS type
prs_types = ['gm', 'wm', 'fc']

# Create PRS data dataframe
for prs_type in prs_types:
    file_path = 'data/{}_prs.profile'.format(prs_type)
    df_prs = pd.read_csv(file_path, delim_whitespace=True)
    df_prs = df_prs[['IID', 'SCORE']]

    # Convert SCORE to z-score
    df_prs['ZSCORE'] = (df_prs['SCORE'] - df_prs['SCORE'].mean()) / df_prs['SCORE'].std()

    # Rename columns
    df_prs.rename(columns={'SCORE': '{}_SCORE'.format(prs_type.upper()), 'ZSCORE': '{}_ZSCORE'.format(prs_type.upper())}, inplace=True)

    # Merge dataframes
    if prs_type == 'gm':
        df = df_prs
    else:
        df = pd.merge(df, df_prs, on='IID', how='outer')


# Obrain pTAU data make BID coulmn index and only extract ORRES column
df_ptau = pd.read_csv('data/plasma/biomarker_pTau217.csv', usecols=['BID', 'VISCODE', 'ORRES'])
df_ptau.rename(columns={'ORRES': 'pTau217', 'BID': 'IID'}, inplace=True)
# Keep only those with VISCODE 6
df_ptau = df_ptau[df_ptau['VISCODE'] == 6]
df_ptau.drop(columns=['VISCODE'], inplace=True)
# Convert ptau to float and remove nans
df_ptau['pTau217'] = pd.to_numeric(df_ptau['pTau217'], errors='coerce')
df_ptau.dropna(inplace=True)

# Obtain AB Test data 
df_ab = pd.read_csv('data/plasma/biomarker_AB_Test.csv', usecols=['BID', 'LBTESTCD', 'LBORRES'])
df_ab.rename(columns={'BID': 'IID'}, inplace=True)
# Pivot to wide dataframe and remove duplicates
df_ab = df_ab.drop_duplicates(subset=['IID', 'LBTESTCD'], keep='first')
df_ab = df_ab.pivot(index='IID', columns='LBTESTCD', values='LBORRES')
# Convert to floats
df_ab = df_ab.apply(pd.to_numeric, errors='coerce')

# Obtain Plasma Roche measures
df_roche = pd.read_csv('data/plasma/biomarker_Plasma_Roche_Results.csv', usecols=['BID', 'LBTESTCD', 'LABRESN'])
df_roche.rename(columns={'BID': 'IID'}, inplace=True)

# Pivot to wide dataframe and remove duplicates
df_roche = df_roche.drop_duplicates(subset=['IID', 'LBTESTCD'], keep='first')
df_roche = df_roche.pivot(index='IID', columns='LBTESTCD', values='LABRESN')
# Convert to floats
df_roche = df_roche.apply(pd.to_numeric, errors='coerce')

# Obtain Cognition data
df_cog = pd.read_csv('data/cognition/PACC.csv', usecols=['BID', 'VISCODE', 'FCTOTAL96'])
df_cog.rename(columns={'BID': 'IID'}, inplace=True)
# Keep only those with VISCODE 001
df_cog = df_cog[df_cog['VISCODE'] == 1]
df_cog.drop(columns=['VISCODE'], inplace=True)

# Merge all dataframes
df = pd.merge(df, df_ptau, on='IID', how='outer')
df = pd.merge(df, df_cog, on='IID', how='outer')
df = pd.merge(df, df_ab, on='IID', how='outer')
df = pd.merge(df, df_roche, on='IID', how='outer')

# Save dataframe
df.to_csv('data/processed_data.csv', index=False)