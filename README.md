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

To include a **"How to Use"** section in your README that explains the technical setup on a local computer, use the following format.

This guide assumes the user has a Python environment ready.

---

## 💻 Installation & Usage Guide

Follow these steps to set up and run the Synthetic Accessibility (SA) Scorer on your local machine.

### 1. Environment Setup

First, ensure you have Python 3.8+ installed. It is recommended to use a virtual environment:

```bash
# Create and activate a virtual environment
python -m venv sa_env
source sa_env/bin/activate  # On Windows: sa_env\Scripts\activate

# Install required dependencies
pip install rdkit pandas numpy tqdm
```

### 2. Prepare the Reference Data (First Run Only)

The scorer requires a fragment frequency library to determine what "rare" chemistry looks like.

1. Download a large SMILES dataset (e.g., PubChem) and name it `pubchem_10m.txt.zip`.
2. Run the **Fragment Counting** portion of the script.
3. This will generate a file named `freq_data.csv`. **Do not delete this file**, as the scorer needs it for every calculation.

### 3. Running the Scorer

To score your own molecules, follow these steps:

1. **Prepare your input:** Create a CSV file (e.g., `molecules_to_score.csv`) with a column named `smiles`.
2. **Update Paths:** Open the script and update the file paths in the `main()` function:
```python
# Update these paths to match your computer's folders
freq_data_path = 'C:/Users/Name/Downloads/freq_data.csv' or '/home/Name/Downloads/freq_data.csv'
input_smiles_path = 'C:/Users/Name/Downloads/target_data.csv' or '/home/Name/Downloads/target_data.csv'
```


3. **Execute:**
```bash
python sa_scorer.py
```



### 4. How it Works (System Architecture)

The script uses a multi-stage pipeline to ensure efficiency on your local hardware:

1. **Data Loading:** Loads the `freq_data.csv` into memory as a lookup dictionary.
2. **Parallelization:** The script detects your CPU core count (e.g., 8 or 16 cores) and splits your SMILES list into "chunks."
3. **Processing:** Each core independently calculates the fragment and complexity penalties.
4. **Aggregation:** Results are collected and displayed as an average score and total processing time.

---

### 📂 File Requirements

* **`pubchem_10m.txt.zip`**: Dataset for creating freq_data.csv, Just in case, [Here is the dataset I used](https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/pubchem_10m.txt.zip) 
* **`freq_data.csv`**: The "knowledge base" created from the reference set.
* **`target_data.csv`**: Your molecules. Must contain a `smiles` column.

---

### 5. Licence

MIT licence