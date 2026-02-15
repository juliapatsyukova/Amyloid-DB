"""Sequence composition feature computation."""

from collections import Counter
from typing import Dict

from ..models import SequenceComposition
from ..config import STANDARD_AA
from ..utils.sequence import clean_sequence


def compute_sequence_composition(sequence: str) -> SequenceComposition:
    """Compute sequence composition features."""
    comp = SequenceComposition()
    
    seq = clean_sequence(sequence)
    if not seq:
        return comp
    
    n = len(seq)
    counter = Counter(seq)
    
    # Single AA frequencies
    comp.aa_freq = {aa: counter.get(aa, 0) / n for aa in STANDARD_AA}
    
    # Dipeptide frequencies (top 25)
    if n >= 2:
        dipeptides = [seq[i:i+2] for i in range(n - 1)]
        dp_counter = Counter(dipeptides)
        total_dp = len(dipeptides)
        for dp, count in dp_counter.most_common(25):
            comp.dipeptide_freq[dp] = count / total_dp
    
    # Grouped compositions
    comp.tiny_fraction = sum(counter.get(aa, 0) for aa in 'AGS') / n
    comp.small_fraction = sum(counter.get(aa, 0) for aa in 'ACDGNPSTV') / n
    comp.large_fraction = sum(counter.get(aa, 0) for aa in 'FHKRWY') / n
    
    return comp
