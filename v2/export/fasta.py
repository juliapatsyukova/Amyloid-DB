"""FASTA export."""

from typing import List

from ..models import AmyloidEntry
from ..utils.sequence import clean_sequence


def write_fasta(entries: List[AmyloidEntry], output_path: str, 
                deduplicate: bool = True, min_length: int = 4):
    """Export entries to FASTA file."""
    sequences = {}
    
    for entry in entries:
        seq = clean_sequence(entry.sequence)
        if not seq or len(seq) < min_length:
            continue
        if deduplicate and seq in sequences:
            continue
        
        parts = [entry.source_db]
        if entry.record_id:
            parts.append(entry.record_id)
        if entry.protein_name:
            parts.append(entry.protein_name.replace(' ', '_').replace('|', '-'))
        if entry.uniprot_id:
            parts.append(f"UniProt:{entry.uniprot_id}")
        if entry.region_start and entry.region_end:
            parts.append(f"region:{entry.region_start}-{entry.region_end}")
        if entry.protein_family:
            parts.append(f"family:{entry.protein_family.replace(' ', '_')[:30]}")
        
        sequences[seq] = "|".join(parts)
    
    with open(output_path, 'w') as f:
        for seq, header in sorted(sequences.items(), key=lambda x: x[1]):
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(f"{seq[i:i+60]}\n")
    
    print(f"Wrote {len(sequences)} sequences to {output_path}")
