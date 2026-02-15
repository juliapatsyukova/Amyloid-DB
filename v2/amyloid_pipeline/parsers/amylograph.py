"""AmyloGraph parser."""

import csv
from typing import List

from .base import BaseParser
from ..models import AmyloidEntry
from ..config import EVIDENCE_WEIGHTS, TIER_CONFIDENCE
from ..utils.sequence import clean_value


class AmyloGraphParser(BaseParser):
    """Parser for AmyloGraph (CSV/TSV)."""
    
    def parse(self, filepath: str, **kwargs) -> List[AmyloidEntry]:
        self.logger.info(f"Parsing {self.db_name} from {filepath}...")
        
        seen_sequences = set()
        
        with open(filepath, 'r') as f:
            first_line = f.readline()
            delimiter = '\t' if '\t' in first_line else ','
        
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            
            for idx, row in enumerate(reader):
                for prefix in ['interactor', 'interactee']:
                    name = clean_value(row.get(f'{prefix}_name', ''))
                    seq = clean_value(row.get(f'{prefix}_sequence', ''))
                    
                    if seq and seq not in seen_sequences:
                        seen_sequences.add(seq)
                        
                        entry = AmyloidEntry(
                            record_id=f"AG_{name}_{len(seen_sequences)}",
                            source_db="AmyloGraph",
                            protein_name=name,
                            protein_family=self._infer_protein_family(name),
                            sequence=seq,
                            doi=clean_value(row.get('doi', '')),
                            category="cross_seeding",
                            is_amyloid=True,
                            experimental_label='amyloid',
                            structure_type='fibril',
                            aggregate_type='in_vitro',
                            pathogenicity='unknown',
                            raw_method="Cross-seeding assay",
                            method_universal="Aggregation kinetics",
                            evidence_type="kinetic",
                            evidence_weight=EVIDENCE_WEIGHTS["kinetic"],
                            confidence=TIER_CONFIDENCE[2],
                            notes=f"AGID: {row.get('AGID', '')}"
                        )
                        
                        if entry.is_valid():
                            self.entries.append(entry)
        
        self.log_result()
        return self.entries
