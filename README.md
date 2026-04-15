# SA-Score: Synthetic Accessibility Scorer

**SA-Score** is a Python tool for estimating the synthetic accessibility of drug-like molecules. By combining fragment-based statistics with structural complexity penalties, it produces a heuristic score (1–10) that helps prioritize compounds in virtual screening and generative design workflows.

Unlike standard implementations, this package allows custom recalibration of fragment penalties on user-defined datasets, making it adaptable to specialized chemical spaces.
## 🧪 Features
* **Fragment-Based Scoring:** Analyzes the frequency of molecular fragments (based on Morgan fingerprints) to identify common vs. "rare" substructures.
* **Complexity Penalties:** Accounts for structural features that complicate synthesis, including:
    * Chiral centers
    * Ring complexity (bridgeheads, spiro atoms, and macrocycles)
    * Molecular size and symmetry
* **RDKit Integration:** Seamlessly processes SMILES strings and RDKit molecule objects.
* **High Efficiency:** High-throughput processing capable of calculating scores for approximately 1,000 molecules per second.

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
If you want to use the default fragment penalties and structural complexity logic:

```python
from SA_score.SaScore import MoleculeProcessor

# Initialize the scorer with default penalties
scorer = MoleculeProcessor()

# Load your molecule
smiles = "CC(=O)NC1=CC=C(C=C1)O"

# Calculate SA Score
score = scorer.SaScorer(smiles)
print(f"SA Score: {score:.3f}")
```

### 2. Advanced Usage (Custom Dataset)
If your research involves a specific chemical space (e.g., natural products or specific macrocycles) and you wish to recalibrate the fragment penalties using your own dataset:

1.  **Prepare your dataset:** Ensure your molecules are in `.zip` or `.gz` format, or in an uncompressed file. The input file should consist of valid SMILES strings, with an optimal dataset size of one million molecules.
2.  **Run the Penalty Calculator:** Generate a new fragment penalty dataset, which will be exported as a CSV file. (this will take ~ 5 minutes) 

```python
from SA_score.FragmentPenaltyFileGenerator import FragmentPenalty

# Initialize the calculator with your custom dataset
calculator = FragmentPenalty(path="your_custom_molecules")

# Generate new fragment penalties
calculator.save_results(filename="your_filename")

# Main usage
from SA_score.SaScore import MoleculeProcessor

custom_scorer = MoleculeProcessor(filepath="your_filepath")

# Load your molecule
smiles = "CC(=O)NC1=CC=C(C=C1)O"

# Calculate SA Score
score = custom_scorer.SaScorer(smiles)
print(f"SA Score: {score:.3f}")
```
### 3. Accessing Detailed Metrics
Beyond the final SA-Score, you can extract individual components such as fragment contributions, complexity penalties, and specific structural counts (for further details regarding complexity and fragment scores, please refer to the [cited work](https://link.springer.com/article/10.1186/1758-2946-1-8):
```python
from SA_score.SaScore import MoleculeProcessor

scorer = MoleculeProcessor()
smiles = "CC(=O)NC1=CC=C(C=C1)O"

# Get individual components
complexity = scorer.ComplexityScore(smiles)
fragments = scorer.FragmentScore(smiles)

# Get structural counts (bridgeheads and spiro atoms)
bridgeheads, spiro_atoms = scorer.bridge(smiles)

print(f"Complexity: {complexity}, Fragments: {fragments}")
print(f"Bridgeheads: {bridgeheads}, Spiro Atoms: {spiro_atoms}")
```
---
## 📊 Validation & Performance
### Accuracy
This implementation has been validated against the original RDKit SA-Score contribution. The correlation between this class and the standard RDKit implementation is 0.975, ensuring high fidelity to the established methodology while providing additional flexibility.

### Benchmark Examples
Below are sample scores calculated using this tool:

| Molecule (SMILES)                                                            | SA-Score |
|:-----------------------------------------------------------------------------|:---------|
| `C1CCCCC1`                                                                   | 1.010    |
| `CC(=O)NC1=CC=C(C=C1)O`                                                      | 1.363    | 
| `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`                                               | 2.368    | 
| `CO[C@H]1CC[C@]2(CC1)Cc1ccccc1C21N=C(N)N(C)C1=O`                             | 4.012    |
| `C/C=C1/CN2CC[C@@]34C(=Nc5c(OC)cccc53)O/C=C(/C(=O)OC)[C@H]1CC24`             | 5.258    |

### Interpretation
* **1.0:** Very easy to synthesize (simple, common fragments).
* **10.0:** Extremely difficult to synthesize (complex, rare fragments, high chirality).


## 📚 References
This implementation is based on the methodology described by **Ertl and Schuffenhauer**:
> Ertl, P., Schuffenhauer, A. "Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions." *Journal of Cheminformatics* 1, 8 (2009). [DOI: 10.1186/1758-2946-1-8](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8)

## 🤝 Contributing
Contributions are welcome! If you find a bug or have a feature request, please open an issue or submit a pull request.

## 📜 License
This project is licensed under the MIT License - see the `LICENSE` file for details.
