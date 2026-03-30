from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from collections import Counter
import pandas as pd
import zipfile
from rdkit import RDLogger
import numpy as np

# Disabling RDKit logging
RDLogger.DisableLog('rdApp.*')

class FragmentPenaltyFile:
    def __init__(self, zip_path, limit=1000000):
        self.zip_path = zip_path
        self.limit = limit
        self.smiles_list = []
        self.freq_df = None

    def load_and_filter_smiles(self):
        """Extracts SMILES from zip and filters by Molecular Weight."""
        count = 0
        with zipfile.ZipFile(self.zip_path) as z:
            inner = z.namelist()[0]
            with z.open(inner, "r") as f:
                for raw in f:
                    s = raw.decode().strip()
                    mol = Chem.MolFromSmiles(s)
                    if mol is None:
                        continue

                    # Using molecules which have molecular weights from 100 to 700
                    mw = Descriptors.MolWt(mol)
                    if 100 <= mw <= 700:
                        self.smiles_list.append(s)

                    count += 1
                    if count >= self.limit:
                        break
        return self.smiles_list

    def count_fragments(self, radius=2):
        """Calculates the number of molecules in which each fragment appears."""
        fragment_counter = Counter()
        for smi in self.smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprint(mol, radius)
            frags = fp.GetNonzeroElements()
            fragment_counter.update(frags.keys())

        self.freq_df = pd.DataFrame(fragment_counter.items(), columns=['fragment', 'count'])
        self.freq_df = self.freq_df.sort_values(by='count', ascending=False)
        return self.freq_df

    def calculate_penalties(self):
        """Calculates fragment penalty based on frequency distribution."""
        if self.freq_df is None:
            raise ValueError("Fragment counts must be calculated before penalties.")

        summ = self.freq_df['count'].sum()
        # Adding temp cumsum for threshold calculation
        self.freq_df['cumsum'] = self.freq_df['count'].cumsum()

        threshold_idx = self.freq_df[self.freq_df['cumsum'] >= summ * 0.8].index[0]

        self.freq_df['fragment_penalty'] = np.log10(self.freq_df['count'] / threshold_idx)
        self.freq_df = self.freq_df.drop(columns=['cumsum'])
        return self.freq_df

    def save_results(self, filename='freq_data.csv'):
        """Saves the final dataframe to CSV."""
        if self.freq_df is not None:
            self.freq_df.to_csv(filename, index=False)
        else:
            print("No data to save.")
