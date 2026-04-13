from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import time
import warnings
from rdkit import RDLogger

# Suppress warnings and logs for cleaner output
warnings.filterwarnings('ignore')
RDLogger.DisableLog('rdApp.*')

def load_fragment_data(filepath: str) -> tuple:
    """Load fragment frequency and penalty data from CSV file."""
    freq_df = pd.read_csv(filepath)

    freq_dict = dict(zip(freq_df['fragment'], freq_df['count']))
    penalty_dict = dict(zip(freq_df['fragment'], freq_df['fragment_penalty']))

    return freq_df, freq_dict, penalty_dict


# Load fragment data
freq_df, freq_dict, penalty_dict = load_fragment_data('freq_data.csv')
DEFAULT_PENALTY = freq_df.iloc[-1, 2]


def fragment_info_for_molecule(smiles: str, radius: int = 2) -> pd.DataFrame:
    """Calculate fragment information for a given SMILES string."""
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprint(mol, radius)
    local_frags = list(fp.GetNonzeroElements().keys())

    # Calculate counts and penalties
    counts = [freq_dict.get(frag, 0) for frag in local_frags]
    penalties = [penalty_dict.get(frag, 0) for frag in local_frags]

    # Apply default penalty for fragments with zero count
    for i, count in enumerate(counts):
        if count == 0:
            penalties[i] = DEFAULT_PENALTY

    return pd.DataFrame({
        'fragment': local_frags,
        'count': counts,
        'fragment_penalty': penalties
    })

def is_bridged(smiles: str) -> tuple:
    """Calculate number of bridgehead and spiro atoms in a molecule."""
    mol = Chem.MolFromSmiles(smiles)
    ring_info = mol.GetRingInfo()
    rings = list(ring_info.AtomRings())

    output = []

    # Analyze ring pairs
    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            shared_atoms = list(set(rings[i]) & set(rings[j]))
            atoms_in_ring = list(set(rings[i]) | set(rings[j]))

            if len(shared_atoms) == 2:
                # Check if rings are fused
                atom = mol.GetAtomWithIdx(shared_atoms[0])
                neighbor_indices = [nbr.GetIdx() for nbr in atom.GetNeighbors()]

                if shared_atoms[1] in neighbor_indices:
                    output.append('Fused')

            elif len(shared_atoms) == 1:
                output.append('Spiro')

            elif len(shared_atoms) == 0:
                output.append('No')

            else:
                # Multiple shared atoms - potential bridgeheads
                potential_bridgeheads = []
                for atom_idx in shared_atoms:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]

                    # Count neighbors in the combined ring system
                    neighbors_in_ring = sum(1 for nbr in neighbors if nbr in atoms_in_ring)
                    if neighbors_in_ring >= 3:
                        potential_bridgeheads.append(atom_idx)

                output.append(potential_bridgeheads)

    # Find atoms in 3 or more rings
    atom_ring_counts = [0] * mol.GetNumAtoms()
    atom_rings_map = [[] for _ in range(mol.GetNumAtoms())]

    for ring_idx, ring in enumerate(rings):
        for atom_idx in ring:
            atom_ring_counts[atom_idx] += 1
            atom_rings_map[atom_idx].append(ring_idx)

    high_coordination_atoms = [
        i for i, count in enumerate(atom_ring_counts)
        if count >= 3
    ]

    # Additional bridgehead analysis
    all_bridgehead_candidates = []
    additional_candidates = []

    for item in output:
        if isinstance(item, list):
            all_bridgehead_candidates.extend(item)
        else:
            all_bridgehead_candidates.append(item)

    for atom_idx in high_coordination_atoms:
        if atom_idx in all_bridgehead_candidates:
            continue

        additional_candidates.append(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbor_indices = [nbr.GetIdx() for nbr in atom.GetNeighbors()]

        # Check if all neighbors share the same rings
        same_ring_count = 0
        for neighbor_idx in neighbor_indices:
            if all(ring in atom_rings_map[atom_idx] for ring in atom_rings_map[neighbor_idx]):
                same_ring_count += 1

        if same_ring_count != len(neighbor_indices):
            continue
        else:
            additional_candidates.append(neighbor_indices)

    # Add additional candidates to output
    for candidate in additional_candidates:
        if isinstance(candidate, list):
            output.extend(candidate)
        else:
            output.append(candidate)

    # Final counts
    bridgehead_atoms = []
    spiro_count = 0

    for element in output:
        if isinstance(element, list):
            bridgehead_atoms.extend(element)
        elif element == 'Spiro':
            spiro_count += 1

    bridgehead_count = len(set(bridgehead_atoms))

    return bridgehead_count, spiro_count


def complexity_penalty_calculator(smiles: str) -> float:
    """Calculate molecular complexity penalty score."""
    mol = Chem.MolFromSmiles(smiles)

    # Stereocenter complexity
    stereo_centers = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    stereo_complexity = np.log10(stereo_centers + 1)

    # Macrocycle penalty
    ring_info = mol.GetRingInfo()
    rings = list(ring_info.AtomRings())

    macrocycle_count = sum(1 for ring in rings if len(ring) > 8)
    macrocycle_penalty = np.log10(macrocycle_count + 1)

    # Size penalty
    atom_count = mol.GetNumAtoms()
    size_penalty = atom_count ** 1.005 - atom_count

    # Ring complexity
    bridgehead_count, spiro_count = is_bridged(smiles)
    ring_complexity = np.log10(bridgehead_count + 1) + np.log10(spiro_count + 1)

    return macrocycle_penalty + size_penalty + ring_complexity + stereo_complexity


def compute_sa_score(smiles: str) -> float:
    """Calculate synthetic accessibility score for a molecule."""
    try:
        # Fragment penalty
        fragment_data = fragment_info_for_molecule(smiles)
        if fragment_data is None or 'fragment_penalty' not in fragment_data:
            return None

        fragment_penalty = fragment_data['fragment_penalty'].mean()

        # Complexity penalty
        complexity_penalty = complexity_penalty_calculator(smiles)

        # Raw score calculation
        raw_score = (fragment_penalty - complexity_penalty) * (-1) + 1.7

        # Scale to 1-10 range
        scaled_score = 1 + (raw_score / 8) * 9

        # Clamp to valid range
        if scaled_score < 1:
            return 1.01
        elif scaled_score > 10:
            return 9.99
        else:
            return scaled_score

    except Exception as e:
        return None


def process_molecules_in_parallel(smiles_list: list, num_workers: int = None) -> list:
    """Process molecules in parallel using multiprocessing."""
    if num_workers is None:
        num_workers = max(cpu_count() - 1, 1)

    with Pool(num_workers) as pool:
        results = list(tqdm(
            pool.imap(compute_sa_score, smiles_list),
            total=len(smiles_list),
            desc="Processing molecules",
            unit="mol"
        ))

    return results


def main():
    """Main execution function."""
    # Load SMILES data
    smiles_data = pd.read_csv('chembl_1000000_random.csv')
    smiles_list = smiles_data['smiles'].tolist()[:1000000]

    # Process molecules
    start_time = time.time()
    scores = process_molecules_in_parallel(smiles_list)
    end_time = time.time()

    # Output results
    processing_time = end_time - start_time
    valid_scores = [s for s in scores if s is not None]

    print(f"\n{'=' * 50}")
    print("Processing Complete!")
    print(f"{'=' * 50}")
    print(f"Total time: {processing_time:.2f} seconds")
    print(f"Processed molecules: {len(scores)}")
    print(f"Valid scores obtained: {len(valid_scores)}")
    print(f"Average score: {np.mean(valid_scores):.2f}" if valid_scores else "No valid scores")
    print(f"{'=' * 50}")


if __name__ == "__main__":
    main()