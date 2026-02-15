"""Physicochemical feature computation."""

import math
import re
from collections import Counter
from typing import Dict

from ..models import PhysicochemicalFeatures
from ..utils.sequence import clean_sequence

# Kyte-Doolittle hydrophobicity
HYDROPHOBICITY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# Chou-Fasman β-sheet propensity
BETA_PROPENSITY = {
    'A': 0.83, 'R': 0.93, 'N': 0.89, 'D': 0.54, 'C': 1.19,
    'Q': 1.10, 'E': 0.37, 'G': 0.75, 'H': 0.87, 'I': 1.60,
    'L': 1.30, 'K': 0.74, 'M': 1.05, 'F': 1.38, 'P': 0.55,
    'S': 0.75, 'T': 1.19, 'W': 1.37, 'Y': 1.47, 'V': 1.70
}

# Residue charge at pH 7
CHARGE = {
    'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0,
    'Q': 0, 'E': -1, 'G': 0, 'H': 0.1, 'I': 0,
    'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
    'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0
}

# Molecular weight
MW = {
    'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
    'Q': 146.2, 'E': 147.1, 'G': 75.1, 'H': 155.2, 'I': 131.2,
    'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
    'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
}

# Aggregation propensity
AGGREGATION_PROPENSITY = {
    'A': 0.3, 'R': -1.0, 'N': -0.5, 'D': -1.0, 'C': 0.5,
    'Q': -0.3, 'E': -1.0, 'G': 0.0, 'H': -0.5, 'I': 1.5,
    'L': 1.2, 'K': -1.0, 'M': 0.8, 'F': 1.8, 'P': -2.0,
    'S': -0.2, 'T': 0.0, 'W': 1.5, 'Y': 1.2, 'V': 1.3
}


def compute_physicochemical_features(sequence: str) -> PhysicochemicalFeatures:
    """Compute physicochemical features for a sequence."""
    features = PhysicochemicalFeatures()
    
    seq = clean_sequence(sequence)
    if not seq:
        return features
    
    features.length = len(seq)
    features.molecular_weight = sum(MW.get(aa, 0) for aa in seq) - (len(seq) - 1) * 18.015
    
    # Hydrophobicity
    hydro_values = [HYDROPHOBICITY.get(aa, 0) for aa in seq]
    if hydro_values:
        features.hydrophobicity_mean = sum(hydro_values) / len(hydro_values)
        features.hydrophobicity_max = max(hydro_values)
        features.hydrophobicity_min = min(hydro_values)
        if len(hydro_values) > 1:
            mean = features.hydrophobicity_mean
            features.hydrophobicity_std = math.sqrt(
                sum((x - mean) ** 2 for x in hydro_values) / len(hydro_values)
            )
    
    # Charge
    charges = [CHARGE.get(aa, 0) for aa in seq]
    features.net_charge = sum(charges)
    features.positive_residues = sum(1 for aa in seq if aa in 'RKH')
    features.negative_residues = sum(1 for aa in seq if aa in 'DE')
    features.charge_density = features.net_charge / len(seq)
    
    # β-sheet propensity
    beta_values = [BETA_PROPENSITY.get(aa, 1.0) for aa in seq]
    if beta_values:
        features.beta_propensity_mean = sum(beta_values) / len(beta_values)
        features.beta_propensity_max = max(beta_values)
    
    # Aggregation propensity
    agg_values = [AGGREGATION_PROPENSITY.get(aa, 0) for aa in seq]
    if agg_values:
        features.aggregation_propensity_mean = sum(agg_values) / len(agg_values)
        features.aggregation_propensity_max = max(agg_values)
    
    # Composition fractions
    counter = Counter(seq)
    n = len(seq)
    features.aromatic_fraction = sum(counter.get(aa, 0) for aa in 'FWY') / n
    features.aliphatic_fraction = sum(counter.get(aa, 0) for aa in 'AVIL') / n
    features.polar_fraction = sum(counter.get(aa, 0) for aa in 'STNQ') / n
    features.charged_fraction = sum(counter.get(aa, 0) for aa in 'DEKRH') / n
    
    # Sequence patterns
    features.has_polyq = bool(re.search(r'Q{3,}', seq))
    features.has_glycine_rich = counter.get('G', 0) / n > 0.2
    features.has_proline = 'P' in seq
    
    return features
