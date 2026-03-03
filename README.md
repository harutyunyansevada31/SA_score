# Custom Synthetic Accessibility (SA) Score Calculator
> A high-performance Python implementation for estimating molecular synthetic difficulty based on fragment frequency and structural complexity.

## 🔬 Overview
This project implements a Synthetic Accessibility (SA) scoring system inspired by Ertl and Schuffenhauer’s methodology but customized with a data-driven approach. It calculates a score between **1 (easy to synthesize)** and **10 (highly complex)** by balancing two major factors:

1.  **Fragment Score:** Based on the frequency of molecular fragments (Morgan Fingerprints) found in a training set of 1,000,000 PubChem molecules.
2.  **Complexity Penalty:** A mathematical penalty for "difficult" structural features such as bridgehead atoms, spiro centers, macrocycles, stereocenters, and high atom counts.

## 🛠️ Key Technical Features
* **Data-Driven Fragment Analysis:** Includes a preprocessing script to generate `freq_data.csv` from large-scale chemical databases (supports `.zip` and `.csv`).
* **Advanced Topology Detection:** Custom logic to identify bridged, fused, and spiro ring systems.
* **Parallel Processing:** Built with Python's `multiprocessing` and `tqdm` to handle datasets of 1,000,000+ molecules efficiently.
* **RDKit Integration:** Uses RDKit for high-fidelity molecular descriptor calculation and SMILES parsing.

## 📂 Repository Structure
* `SA score (main part).py`: The primary engine. Includes the complexity penalty calculator and the parallel processing pipeline.
* `fragment_generation.py`: (The first part of your code) Used to generate the fragment penalty reference file (`freq_data.csv`).
* `freq_data.csv`: The generated lookup table for fragment penalties (required for scoring).

## 🚀 Getting Started

### Prerequisites
* Python 3.8+
* RDKit
* Pandas
* Numpy
* TQDM

### Installation
```bash
pip install pandas numpy tqdm
# For RDKit (via Conda)
conda install -c rdkit rdkit