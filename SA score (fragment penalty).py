from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdMolDescriptors, Descriptors
from collections import Counter
import pandas as pd
import zipfile
from rdkit import RDLogger
import numpy as np
# Downloading the SMILES of the molecules, based on which we will create fragment score data.
RDLogger.DisableLog('rdApp.*')

limit = 1000000
count = 0
smiles_list = []
with zipfile.ZipFile("/home/denovo/Downloads/pubchem_10m.txt.zip") as z: # You can change this to your 1 million SMILES data
    inner = z.namelist()[0]
    with z.open(inner, "r") as f:
        for raw in f:
            s = raw.decode().strip()
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                continue
            # using molecules which have molecular weights from 100 to 700
            mw = Descriptors.MolWt(mol)
            if 100 <= mw <= 700:
                smiles_list.append(s)
            count += 1
            if count >= limit:
                break
df = pd.DataFrame(smiles_list)

# Calculating the number of molecules in which each fragment appears.
def count_fragments(smiles_list, radius=2):
    fragment_counter = Counter()
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprint(mol, radius)
        frags = fp.GetNonzeroElements()
        fragment_counter.update(frags.keys())
    return fragment_counter
fragment_counts = count_fragments(smiles_list)
freq_df = pd.DataFrame(fragment_counts.items(), columns=['fragment', 'count'])
freq_df = freq_df.sort_values(by='count', ascending=False)

# Calculating fragment penalty
summ = freq_df['count'].sum()
freq_df['cumsum'] = freq_df['count'].cumsum()
threshold_idx = freq_df[freq_df['cumsum'] >= summ * 0.8].index[0]
mean = freq_df.loc[threshold_idx, 'count']
freq_df['fragment_penalty'] = np.log10(freq_df['count'] / threshold_idx)
freq_df = freq_df.drop(columns=['cumsum'])

# Saving our dataframe
freq_df.to_csv('freq_data.csv', index=False)
