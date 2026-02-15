"""AmyLoad parser."""

import csv
from typing import List

from .base import BaseParser
from ..models import AmyloidEntry
from ..config import EVIDENCE_WEIGHTS, TIER_CONFIDENCE
from ..utils.sequence import clean_value, validate_sequence


class AmyLoadParser(BaseParser):
    """Parser for AmyLoad (CSV)."""
    
    def parse(self, filepath: str, **kwargs) -> List[AmyloidEntry]:
        self.logger.info(f"Parsing {self.db_name} from {filepath}...")
        
        with open(filepath, 'r', encoding='utf-8') as f:
            first_line = f.readline()
            has_header = 'sequence' in first_line.lower() or 'protein' in first_line.lower()
        
        with open(filepath, 'r', encoding='utf-8') as f:
            if has_header:
                reader = csv.DictReader(f)
                for idx, row in enumerate(reader):
                    try:
                        seq = clean_value(row.get('Sequence', '')).upper()
                        if not validate_sequence(seq):
                            continue
                        
                        amyloidogenic = clean_value(row.get('Amyloidogenic', '')).lower()
                        if amyloidogenic == 'yes':
                            label, is_amyloid = 'amyloid', True
                        elif amyloidogenic == 'no':
                            label, is_amyloid = 'non-amyloid', False
                        else:
                            continue
                        
                        protein_name = clean_value(row.get('Protein Name', ''))
                        
                        entry = AmyloidEntry(
                            record_id=clean_value(row.get('AMY_ID', f"AL_{idx}")),
                            source_db="AmyLoad",
                            protein_name=protein_name,
                            protein_family=self._infer_protein_family(protein_name),
                            sequence=seq,
                            is_amyloid=is_amyloid,
                            experimental_label=label,
                            category="experimental",
                            structure_type='unknown',
                            aggregate_type='unknown',
                            pathogenicity=self._infer_pathogenicity("", "", protein_name),
                            raw_method="Literature-curated",
                            method_universal="Literature-curated",
                            evidence_type="literature_curated",
                            evidence_weight=EVIDENCE_WEIGHTS["literature_curated"],
                            confidence=TIER_CONFIDENCE[0]
                        )
                        
                        if entry.is_valid():
                            self.entries.append(entry)
                    except Exception as e:
                        self.errors.append(f"Row {idx}: {str(e)}")
            else:
                reader = csv.reader(f)
                for idx, row in enumerate(reader):
                    if len(row) < 5:
                        continue
                    amy_id, protein, fragment, seq, is_amyloid_str = row[:5]
                    seq = clean_value(seq).upper()
                    if not validate_sequence(seq):
                        continue
                    
                    is_amyloid = clean_value(is_amyloid_str).lower() == 'yes'
                    protein_name = clean_value(protein)
                    
                    entry = AmyloidEntry(
                        record_id=clean_value(amy_id),
                        source_db="AmyLoad",
                        protein_name=protein_name,
                        protein_family=self._infer_protein_family(protein_name),
                        sequence=seq,
                        is_amyloid=is_amyloid,
                        experimental_label='amyloid' if is_amyloid else 'non-amyloid',
                        category="experimental",
                        structure_type='unknown',
                        aggregate_type='unknown',
                        pathogenicity=self._infer_pathogenicity("", "", protein_name),
                        raw_method="Literature-curated",
                        method_universal="Literature-curated",
                        evidence_type="literature_curated",
                        evidence_weight=EVIDENCE_WEIGHTS["literature_curated"],
                        confidence=TIER_CONFIDENCE[0],
                        notes=f"Fragment: {clean_value(fragment)}"
                    )
                    
                    if entry.is_valid():
                        self.entries.append(entry)
        
        self.log_result()
        return self.entries
