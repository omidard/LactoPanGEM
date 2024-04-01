import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from multiprocessing import Pool

# Load your existing data
inf2 = pd.read_csv('/home/omidard/niche/all_inf.csv')
inf2.drop(columns=['Unnamed: 0'], inplace=True)
inf2.set_index('id', inplace=True)

pam = pd.read_csv('/home/omidard/niche/all_reactions.csv')
pam.fillna(0, inplace=True)
pam.set_index('Unnamed: 0', inplace=True)
pamt = pam.T

# Remove '.json' from indexes
pamt.index = pamt.index.str.replace('.json', '', regex=True)

# Add isolation_source from inf2 to pamt
pamt = pamt.merge(inf2[['isolation_source']], left_index=True, right_index=True, how='left')

# Remove non-informative isolation sources
non_informative_sources = ['Not reported', 'missing', '-', 'unknown']
pamt = pamt[~pamt['isolation_source'].isin(non_informative_sources)]

# Keywords to search for in isolation sources
keywords = ['food', 'vagina', 'blood', 'kefir', 'yogurt', 'cheese', 'cream', 'gut', 'gastrointestinal']
synonyms = {'feces': 'fecal sample', 'stool': 'fecal sample'}

# Function to rename isolation sources based on the keyword
def rename_isolation_sources(source):
    for keyword, replacement in synonyms.items():
        if keyword in source.lower():
            return replacement
    for keyword in keywords:
        if keyword in source.lower():
            return keyword
    return source

# Apply the function to rename isolation sources
pamt['isolation_source'] = pamt['isolation_source'].apply(rename_isolation_sources)

# Drop reactions (columns) that are present in all GEMs
cols_to_drop = pamt.columns[pamt.sum(axis=0) == len(pamt)]
pamt.drop(columns=cols_to_drop, inplace=True)

# Only consider isolation sources with more than 10 isolates
isolation_counts = pamt['isolation_source'].value_counts()
valid_sources = isolation_counts[isolation_counts > 10].index
pamt = pamt[pamt['isolation_source'].isin(valid_sources)]

def calculate_odds_ratio(data):
    reaction, niche, pamt = data
    a = len(pamt[(pamt[reaction] > 0) & (pamt['isolation_source'] == niche)])
    b = len(pamt[(pamt[reaction] == 0) & (pamt['isolation_source'] == niche)])
    c = len(pamt[(pamt[reaction] > 0) & (pamt['isolation_source'] != niche)])
    d = len(pamt[(pamt[reaction] == 0) & (pamt['isolation_source'] != niche)])
    table = [[a, b], [c, d]]
    odds_ratio, p_value = fisher_exact(table)
    return (reaction, niche, odds_ratio, p_value)

# Get the list of unique niches/isolation sources
niches = pamt['isolation_source'].unique()

# Preparing the data for multiprocessing
data_for_processing = [(reaction, niche, pamt) for niche in niches for reaction in pamt.columns[:-1]]  # Exclude the 'isolation_source' column

# Use Pool to parallelize the calculations
with Pool(64) as pool:
    results = pool.map(calculate_odds_ratio, data_for_processing)

# Filter reactions with significant p-values (e.g., p < 0.05) and high odds ratios
some_threshold = 2  # Adjust this value as required
significant_reactions = [result for result in results if result[3] < 0.05 and result[2] > some_threshold]

# Convert results to DataFrame and save
df = pd.DataFrame(significant_reactions, columns=['Reaction', 'Niche', 'Odds Ratio', 'P-value'])
df.to_csv('/home/omidard/niche/significant_reactions.csv', index=False)
