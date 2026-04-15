from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import pandas as pd
import numpy as np
import warnings
from rdkit import RDLogger
import os
warnings.filterwarnings('ignore')
RDLogger.DisableLog('rdApp.*')

class MoleculeProcessor:
    def __init__(self, filepath=None):
        if filepath is None:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            parent_dir = os.path.dirname(current_dir)
            self.filepath = os.path.join(parent_dir, 'freq_data.csv')
        else:
            self.filepath = filepath
        freq_df = pd.read_csv(self.filepath)
        self.freq_dict = dict(zip(freq_df['fragment'], freq_df['count']))
        self.penalty_dict = dict(zip(freq_df['fragment'], freq_df['fragment_penalty']))
        self.default_penalty = freq_df.iloc[-1, 2]
        self.radius = 2
        self.mol = None
        self.smiles = None
        self.bridgehead = 0
        self.spiro = 0

    def _prepare_mol(self, smiles):
        """Internal helper to set up the molecule and calculate bridge/spiro counts."""
        self.mol = Chem.MolFromSmiles(smiles)
        if self.mol:
            # We call bridge automatically so the counts are ready for ComplexityScore
            self.bridge(smiles)
        return self.mol

    def bridge(self, smiles):
        ring_info = self.mol.GetRingInfo()
        rings = list(ring_info.AtomRings())
        output = []
        # Analyze ring pairs
        for i in range(len(rings)):
            for j in range(i + 1, len(rings)):
                shared_atoms = list(set(rings[i]) & set(rings[j]))
                atoms_in_ring = list(set(rings[i]) | set(rings[j]))

                if len(shared_atoms) == 2:
                    # Check if rings are fused
                    atom = self.mol.GetAtomWithIdx(shared_atoms[0])
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
                        atom = self.mol.GetAtomWithIdx(atom_idx)
                        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]

                        # Count neighbors in the combined ring system
                        neighbors_in_ring = sum(1 for nbr in neighbors if nbr in atoms_in_ring)
                        if neighbors_in_ring >= 3:
                            potential_bridgeheads.append(atom_idx)

                    output.append(potential_bridgeheads)

        # Find atoms in 3 or more rings
        atom_ring_counts = [0] * self.mol.GetNumAtoms()
        atom_rings_map = [[] for _ in range(self.mol.GetNumAtoms())]

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
            atom = self.mol.GetAtomWithIdx(atom_idx)
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

        self.bridgehead = bridgehead_count
        self.spiro = spiro_count
        return bridgehead_count, spiro_count

    def ComplexityScore(self, smiles):
        if not self._prepare_mol(smiles):
            return 0.0
        # Stereocenter complexity
        stereo_centers = rdMolDescriptors.CalcNumAtomStereoCenters(self.mol)
        stereo_complexity = np.log10(stereo_centers + 1)
        # Macrocycle penalty
        ring_info = self.mol.GetRingInfo()
        rings = list(ring_info.AtomRings())
        macrocycle_count = sum(1 for ring in rings if len(ring) > 8)
        macrocycle_penalty = np.log10(macrocycle_count + 1)
        # Size penalty
        atom_count = self.mol.GetNumAtoms()
        size_penalty = atom_count ** 1.005 - atom_count
        # Ring complexity
        bridgehead_count, spiro_count = self.bridgehead, self.spiro
        ring_complexity = np.log10(bridgehead_count + 1) + np.log10(spiro_count + 1)

        return macrocycle_penalty + size_penalty + ring_complexity + stereo_complexity

    def FragmentScore(self, smiles):
        """Calculate fragment information for a given SMILES string."""
        if self.mol is None:
            return None
        fp = AllChem.GetMorganFingerprint(self.mol, self.radius)
        local_frags = list(fp.GetNonzeroElements().keys())

        counts = [self.freq_dict.get(frag, 0) for frag in local_frags]
        penalties = [self.penalty_dict.get(frag, 0) for frag in local_frags]

        for i, count in enumerate(counts):
            if count == 0:
                penalties[i] = self.default_penalty
        fragment_data = pd.DataFrame({
            'fragment': local_frags,
            'count': counts,
            'fragment_penalty': penalties
        })
        return fragment_data['fragment_penalty'].mean()

    def SaScorer(self, smiles):
        if not self._prepare_mol(smiles):
            return None
        try:
            f_penalty = self.FragmentScore(smiles)
            c_penalty = self.ComplexityScore(smiles)

            raw_score = (f_penalty - c_penalty) * (-1) + 1.7
            scaled_score = 1 + (raw_score / 8) * 9
            return max(1.01, min(9.99, scaled_score))
        except Exception:
            return None