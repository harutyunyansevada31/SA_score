# Molecular Synthetic Accessibility (SA) Scorer

This repository provides a high-performance Python tool for calculating the **Synthetic Accessibility (SA) Score** of chemical compounds. The score estimates how difficult a molecule is to synthesize based on fragment frequency (from a 1-million molecule reference set) and structural complexity (rings, stereo centers, macrocycles, and bridged systems).

## 🚀 Features

* **Fragment-Based Analysis:** Uses Morgan Fingerprints (Radius 2) to build a fragment frequency dictionary.
* **Structural Complexity Penalty:** Factors in:
* **Spiro and Bridgehead atoms** (using custom ring-topology logic).
* **Stereocenters** and **Macrocycles** (>8 atoms).
* **Size Penalty** based on heavy atom count.


* **Parallel Processing:** Optimized with `multiprocessing` to handle large datasets (e.g., 100k+ molecules) quickly.
* **Normalization:** Outputs a score between **1 (Easy to synthesize)** and **10 (Very difficult)**.

---

## 🛠 Prerequisites

Ensure you have the following libraries installed:

```bash
pip install rdkit pandas numpy tqdm

```

---

## 📂 Code Structure

The implementation consists of two main workflows:

### 1. Fragment Frequency Generation (`Pre-processing`)

The first script processes a large `.zip` file (e.g., `pubchem_10m.txt.zip`) to:

1. Filter molecules by Molecular Weight (100–700 Da).
2. Deconstruct molecules into Morgan fragments.
3. Calculate a **Fragment Penalty** based on the log-frequency of occurrence.
4. Export the reference data to `freq_data.csv`.

### 2. SA Score Calculation (`Main Pipeline`)

The second script loads the `freq_data.csv` and calculates scores for a target dataset:

* **`compute_sa_score(smiles)`**: The core function that aggregates fragment penalties and complexity penalties.
* **`is_bridged(smiles)`**: A deep-dive topological function to identify complex ring systems.
* **`process_molecules_in_parallel`**: Uses a pool of workers to maximize CPU utilization.

---

## 📈 Methodology

The SA Score is calculated as:

SA\_Score = fragment\_penalty - complexity\_penalty

| Component | Description |
| --- | --- |
| **Fragment Score** | Based on the rarity of fragments in the PubChem/Reference set. Rare fragments increase the score. |
| **Size Penalty** | $n^{1.005} - n$, where $n$ is the number of atoms. |
| **Stereo Penalty** | $\log_{10}(\text{chiral centers} + 1)$. |
| **Ring Penalty** | Specifically targets bridged, spiro, and macrocyclic systems. |

---

## 💻 Usage

1. **Generate Frequencies:**
Place your raw SMILES zip file in the directory and run the first block of code to generate `freq_data.csv`.
2. **Score Molecules:**
Update the path to your CSV (e.g., `chembl_1000000_random.csv`) in `main()` and run:
```bash
python sa_scorer.py

```


3. **View Results:**
The script will output processing statistics, including total time, average score, and the number of valid molecules processed.

---

## ⚠️ Notes

* **Default Penalty:** Fragments not found in the reference training set are assigned a "Default Penalty" based on the rarest known fragment.
* **Logging:** RDKit logs are disabled by default (`RDLogger.DisableLog`) to keep the terminal clean during large-scale processing.

Would you like me to help you write a script to visualize the distribution of these scores using Matplotlib or Seaborn?