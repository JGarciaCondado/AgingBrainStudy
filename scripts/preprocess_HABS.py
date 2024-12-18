import pandas as pd

# Load features but ignore first column
features = pd.read_csv('data/HABS/raw/raw_features.csv', index_col='SubjIDshort')

# Remove MRI_FS6_ADNI_ from column names
features.columns = features.columns.str.replace('MRI_FS6_ADNI_', '')

# Combine those with lh and rh 
lh = features.filter(like='lh_').columns
rh = features.filter(like='rh_').columns
for l, r in zip(lh, rh):
    features['bi_' + l.split('_')[1]] = features[l] + features[r]
    features.drop(columns=[l, r], inplace=True)

# Combine those with Left and Right
left = features.filter(like='Left').columns
right = features.filter(like='Right').columns
for l, r in zip(left, right):
    features['bi_' + l.split('_')[1]] = features[l] + features[r]
    features.drop(columns=[l, r], inplace=True)

def normalize_mri(mri_col, icv_col):
    from sklearn.linear_model import LinearRegression
    b = LinearRegression().fit(icv_col.values.reshape(-1,1), mri_col).coef_[0]
    adj_mri = mri_col - b * (icv_col - icv_col.mean())
    return adj_mri
# Normalize volumes by total intracranial volume
features['bi_Hippocampus'] = normalize_mri(features['bi_Hippocampus'], features['ICV_vol'])
features['bi_Amygdala'] = normalize_mri(features['bi_Amygdala'], features['ICV_vol'])

# Keep only MRI_Age and bi_ columns
features = features[['MRI_Age'] + [col for col in features.columns if 'bi_' in col]]

# Rename MRI_Age to age and lower case column names
features = features.rename(columns={'MRI_Age': 'age'})
features.columns = features.columns.str.lower()

# Save as csv
features.to_csv('data/HABS/processed/structural_features.csv')

# Load A4 data 
df_a4 = pd.read_csv('data/A4/structural_features.csv', index_col='BID')
df_a4.columns = df_a4.columns.str.lower()

# Stack dataframes
df = pd.concat([features, df_a4], axis=0)
df.index.name = 'SubjIDshort'
df.to_csv('data/HABS/processed/HABS_A4_structural_features.csv')

# Create new dataframe with CN column with 1 for those in a4 and 0 for those in HABS
df['cn'] = df.index.isin(df_a4.index).astype(int)
clinical = df[['cn']].copy()
clinical['HABS'] = (~clinical['cn'].astype(bool)).astype(int)
df.drop(columns='cn', inplace=True)
clinical.to_csv('data/HABS/processed/A4_as_cn.csv')