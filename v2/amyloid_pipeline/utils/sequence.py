"""Sequence validation and utility functions."""

import re
from typing import List, Tuple, Optional

from ..config import STANDARD_AA, EVIDENCE_WEIGHTS, TIER_CONFIDENCE, METHODS_RULES

# Compile regex patterns once
_COMPILED_RULES = [(re.compile(r["pattern"], re.IGNORECASE), r) for r in METHODS_RULES]


def clean_value(val) -> str:
    """Clean and normalize a value."""
    if val is None:
        return ""
    val = str(val).strip()
    if val.lower() in ['na', 'n/a', 'n.a.', 'none', 'null', '-', '', 'nan']:
        return ""
    return val


def validate_sequence(seq: str) -> bool:
    """Validate amino acid sequence."""
    if not seq:
        return False
    return all(aa.upper() in STANDARD_AA for aa in seq)


def clean_sequence(seq: str) -> str:
    """Clean amino acid sequence, removing invalid characters."""
    if not seq:
        return ""
    return re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq.upper().strip())


def extract_uniprot_id(uniprot_str: str) -> str:
    """Extract clean UniProt ID from various formats."""
    if not uniprot_str:
        return ""
    match = re.match(r'^([A-Z][A-Z0-9]{5,9})', uniprot_str.strip().upper())
    return match.group(1) if match else clean_value(uniprot_str)


def parse_range(range_str: str) -> List[Tuple[Optional[int], Optional[int]]]:
    """Parse range string into list of (start, end) tuples."""
    if not range_str:
        return [(None, None)]
    
    ranges = []
    for part in re.split(r'[,\s]+', range_str.strip()):
        part = part.strip()
        if not part:
            continue
        match = re.match(r'(\d+)\s*[-–]\s*(\d+)', part)
        if match:
            start, end = int(match.group(1)), int(match.group(2))
            if start <= end:
                ranges.append((start, end))
        else:
            match = re.match(r'^(\d+)$', part)
            if match:
                pos = int(match.group(1))
                ranges.append((pos, pos))
    
    return ranges if ranges else [(None, None)]


def map_method_to_universal(raw_method: str) -> Tuple[str, str, float, int]:
    """
    Map raw method string to universal vocabulary.
    
    Returns: (universal_method, evidence_type, evidence_weight, confidence)
    """
    if not raw_method:
        return ("Unspecified", "literature_curated", 
                EVIDENCE_WEIGHTS["literature_curated"], TIER_CONFIDENCE[0])
    
    raw = raw_method.strip()
    
    for regex, rule in _COMPILED_RULES:
        if regex.search(raw):
            evidence_type = rule.get("evidence_type", "staining_binding")
            weight = EVIDENCE_WEIGHTS.get(evidence_type, 0.5)
            confidence = TIER_CONFIDENCE.get(rule["tier"], 10)
            return (rule["universal"], evidence_type, weight, confidence)
    
    return (raw, "literature_curated", EVIDENCE_WEIGHTS["literature_curated"], TIER_CONFIDENCE[0])
