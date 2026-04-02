from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from collections import Counter
import pandas as pd
import numpy as np
from rdkit import RDLogger
import zipfile
import gzip
import io

RDLogger.DisableLog('rdApp.*')

class FragmentPenalty:
    def __init__(self, zip_path, limit=1000000):
        self.zip_path = zip_path
        self.limit = limit
        self.freq_df = None

    def save_results(self, filename='freq_data.csv'):

        fragment_counter = Counter()
        processed_count = 0
        passed_filter_count = 0

        path = self.zip_path

        try:
            if path.endswith('.zip'):
                z = zipfile.ZipFile(path)
                # Zip files always return bytes, so we wrap it in a TextIOWrapper
                raw_handle = z.open(z.namelist()[0], "r")
                f = io.TextIOWrapper(raw_handle, encoding='utf-8', errors='ignore')
            elif path.endswith('.gz'):
                f = gzip.open(path, "rt", encoding='utf-8', errors='ignore')
            else:
                # Regular text/smi file
                f = open(path, "r", encoding='utf-8', errors='ignore')

            with f:
                for s in f:
                    s = s.strip()
                    if not s: continue

                    processed_count += 1
                    mol = Chem.MolFromSmiles(s)

                    if mol:
                        mw = Descriptors.MolWt(mol)
                        if 100 <= mw <= 700:
                            passed_filter_count += 1
                            fp = AllChem.GetMorganFingerprint(mol, 2)
                            fragment_counter.update(fp.GetNonzeroElements().keys())

                    if processed_count >= self.limit:
                        break

        except Exception as e:
            print(f"Failed to process file: {e}")

        if not fragment_counter:
            print(f"No data to save! Processed {processed_count} lines, but 0 molecules passed the MW filter.")
            return

        # 2. Build DataFrame
        df = pd.DataFrame(fragment_counter.items(), columns=['fragment', 'count'])
        df = df.sort_values(by='count', ascending=False).reset_index(drop=True)

        # 3. Calculate Penalties
        df['cumsum'] = df['count'].cumsum()
        total_count = df['count'].sum()
        threshold_idx = df[df['cumsum'] >= total_count * 0.8].index[0]
        df['fragment_penalty'] = np.log10(df['count'] / threshold_idx)
        df = df.drop(columns=['cumsum'])

        # 4. Save
        self.freq_df = df
        self.freq_df.to_csv(filename, index=False)
        print(f"Saved to {filename}")