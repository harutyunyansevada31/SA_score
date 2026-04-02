# SA-Score: Synthetic Accessibility Scorer

**SA-Score** is a Python-based tool designed to estimate the synthetic accessibility of drug-like molecules. By combining fragment-based analysis with structural complexity penalties, it provides a heuristic score (typically between 1 and 10) to help prioritize molecules in virtual screening and generative design.

## 🧪 Features
* **Fragment-Based Scoring:** Analyzes the frequency of molecular fragments (based on Morgan fingerprints) to identify common vs. "rare" substructures.
* **Complexity Penalties:** Accounts for structural features that complicate synthesis, including:
    * Chiral centers
    * Ring complexity (bridgeheads, spiro atoms, and macrocycles)
    * Molecular size and symmetry
* **RDKit Integration:** Seamlessly processes SMILES strings and RDKit molecule objects.
* **Batch Processing:** Support for scoring large libraries of molecules efficiently.

## 🚀 Getting Started

### Prerequisites
* **Python 3.8+**
* **RDKit** (Core dependency for molecular manipulation)

### Installation
Clone the repository and install the necessary requirements:

```bash
git clone https://github.com/harutyunyansevada31/SA_score.git
cd SA_score
pip install -r requirements.txt
```

## 💻 Usage

`SA_score` is designed to be flexible. You can either use the pre-calculated fragment scores (derived from a dataset of 1 million molecules) or generate your own scores based on a custom library.

### 1. Standard Usage (Pre-calculated Data)
If you want to use the default fragment penalties and structural complexity logic, ensure you have downloaded the `freq_data` file. This file contains the pre-calculated fragment penalties derived from our baseline dataset.

```python
from my_library.SaScore import MoleculeProcessor
from rdkit import Chem

# Initialize the scorer with default weights
filepath = 'Your_path_to_freq_data.csv'
scorer = MoleculeProcessor(filepath)

# Load your molecule
smiles = "CC(=O)NC1=CC=C(C=C1)O" 
# Calculate SA Score

score = scorer.SaScorer(smiles)
print(f"SA Score: {score:.3f}")
```

### 2. Advanced Usage (Custom Dataset)
If your research involves a specific chemical space (e.g., natural products or specific macrocycles) and you wish to recalibrate the fragment penalties using your own dataset of 1 million+ molecules:

1.  **Prepare your dataset:** Ensure your molecules are in a format readable by RDKit (e.g., a `.smi` or `.csv` file).
2.  **Run the Penalty Calculator:** Use the training module to generate a new fragment penalty dataset (this will take ~ 5 minutes).

```python
from my_library.FragmentPenaltyFileGenerator import FragmentPenalty

# Initialize the calculator with your custom dataset
calculator = FragmentPenalty(path="your_custom_molecules.smi")

# Generate new fragment penalties
calculator.save_results(filename="your_filename")

# Main usage
from my_library.SaScore import MoleculeProcessor
custom_scorer = MoleculeProcessor(filepath="your_filepath")
```

---

### Interpretation
* **1.0:** Very easy to synthesize (simple, common fragments).
* **10.0:** Extremely difficult to synthesize (complex, rare fragments, high chirality).


## 📚 References
This implementation is based on the methodology described by **Ertl and Schuffenhauer**:
> Ertl, P., Schuffenhauer, A. "Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions." *Journal of Cheminformatics* 1, 8 (2009). [DOI: 10.1186/1758-2946-1-8](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8)

## 🤝 Contributing
Contributions are welcome! If you find a bug or have a feature request (e.g., support for additional fingerprints like MACCS), please open an issue or submit a pull request.

## 📜 License
This project is licensed under the MIT License - see the `LICENSE` file for details.
