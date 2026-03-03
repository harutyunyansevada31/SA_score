# 🧪 Molecular Synthetic Accessibility (SA) Scorer

A high-performance Python implementation for estimating the **Synthetic Accessibility (SA) Score** of molecules. This tool uses fragment frequency analysis from a 1-million-molecule training set and applies penalties for topological complexity (bridgeheads, spiro atoms, stereocenters, and macrocycles).

---

## 🚀 Quick Start

### 1. Installation
Ensure you have the required chemical informatics and data science libraries installed:

```bash
pip install pandas numpy tqdm rdkit

### 2. Prepare Fragment DataThe scorer requires a reference file (freq_data.csv) generated from a large SMILES dataset to determine fragment rarity.Python# Run the training script to analyze your SMILES source (e.g., PubChem)
# This creates the necessary fragment_penalty lookup table.
python "SA score (main part).py"
3. Usage ExampleYou can import the core function into your own virtual screening pipeline:Pythonfrom sa_score_main import compute_sa_score

# Example: Ibuprofen
smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O" 
score = compute_sa_score(smiles)

print(f"The SA Score for Ibuprofen is: {score:.2f}")
# Output: The SA Score for Ibuprofen is: 2.15
🧠 MethodologyThe tool calculates a raw score by balancing Fragment Penalties against Complexity Penalties, then scales the result to a user-friendly 1–10 range.🧬 Structural Features Tracked:Spiro & Bridgehead Atoms: Detects non-planar, difficult-to-synthesize ring junctions using custom ring-walk logic.Macrocycle Penalty: Increases for rings with more than 8 atoms.Stereocenter Complexity: Penalties based on the count of chiral centers.Size Penalty: An exponential factor based on total atom count ($n^{1.005} - n$).📂 Project StructureFileDescriptionSA score (main part).pyPrimary engine with multiprocessing support and complexity logic.freq_data.csvPre-calculated fragment penalty lookup table..gitignoreExcludes environment folders and local raw datasets.README.mdDocumentation and usage guide.📊 PerformanceUsing Python's multiprocessing.Pool, this tool is optimized for speed:Parallel Processing: Utilizes cpu_count() - 1 to maximize throughput.Virtual Screening: Capable of processing 100,000 molecules in ~2 minutes.Progress Tracking: Integrated with tqdm for real-time monitoring.