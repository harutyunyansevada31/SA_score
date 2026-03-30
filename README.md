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
If you want to use the default fragment penalties and structural complexity logic, you can start scoring molecules immediately.

```python
from my_library.SaScore import MoleculeProcessor
from rdkit import Chem

# Initialize the scorer with default weights
scorer = MoleculeProcessor()

# Load your molecule
smiles = "CC(=O)NC1=CC=C(C=C1)O" 
# Calculate SA Score

score = scorer.SaScorer(smiles)
print(f"SA Score: {score:.3f}")
```

### 2. Advanced Usage (Custom Dataset)
If your research involves a specific chemical space (e.g., natural products or specific macrocycles) and you wish to recalibrate the fragment penalties using your own dataset of 1 million+ molecules:

1.  **Prepare your dataset:** Ensure your molecules are in a format readable by RDKit (e.g., a `.smi` or `.csv` file).
2.  **Run the Penalty Calculator:** Use the training module to generate a new fragment weight library.

```python
from my_library.FragmentPenaltyFileGeneration import FragmentPenalty

# Initialize the calculator with your custom dataset
calculator = FragmentPenalty(zip_path="your_custom_molecules.smi")

# Generate new fragment penalties
calculator.save_results(filename="your_filename")

# Main usage
from my_library.SaScore import MoleculeProcessor
custom_scorer = MoleculeProcessor(filepath="your_filepath")
```

---

### Key Improvements Made:
* **Hierarchical Structure:** Used clear subheadings to separate the "plug-and-play" user from the "power user."
* **Object-Oriented Logic:** I framed the code as if you have `SAScorer` and `PenaltyCalculator` classes, which is standard for Python packages. If your class names are different, simply swap them out.
* **The "Why":** I added a brief explanation of *why* someone would choose Path 2 (specialized chemical spaces), which adds professional context to your repository.

**Would you like me to help you write a "Technical Details" section that explains the math behind the fragment penalties or how the complexity score is weighted?**

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
