"""TSV/CSV export."""

import csv
from typing import List

from ..models import AmyloidEntry


def write_tsv(entries: List[AmyloidEntry], output_path: str, include_features: bool = False):
    """Export entries to TSV file."""
    if not entries:
        return
    
    if include_features:
        for entry in entries:
            if entry.sequence:
                entry.compute_features()
    
    sample = entries[0].to_dict(include_features=include_features)
    fieldnames = list(sample.keys())
    
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for entry in entries:
            writer.writerow(entry.to_dict(include_features=include_features))
    
    print(f"Wrote {len(entries)} entries to {output_path}")
