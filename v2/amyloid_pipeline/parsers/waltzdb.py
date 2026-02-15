"""WALTZ-DB parser."""

import csv
from typing import List

from .base import BaseParser
from ..models import AmyloidEntry
from ..config import EVIDENCE_WEIGHTS
from ..utils.sequence import clean_value, validate_sequence, extract_uniprot_id, map_method_to_universal


class WaltzDBParser(BaseParser):
    """Parser for WALTZ-DB 2.0 (CSV or TSV)."""
    
    def parse(self, filepath: str, **kwargs) -> List[AmyloidEntry]:
        self.logger.info(f"Parsing {self.db_name} from {filepath}...")
        
        with open(filepath, 'r') as f:
            first_line = f.readline()
            delimiter = '\t' if '\t' in first_line else ','
        
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            
            for idx, row in enumerate(reader):
                try:
                    seq = clean_value(row.get('Sequence', '') or row.get('sequence', '')).upper()
                    if not validate_sequence(seq):
                        continue
                    
                    classification = clean_value(row.get('Classification', '') or row.get('classification', '')).lower()
                    if 'non' in classification:
                        label, is_amyloid = 'non-amyloid', False
                    elif 'amyloid' in classification:
                        label, is_amyloid = 'amyloid', True
                    else:
                        continue
                    
                    methods = []
                    if clean_value(row.get('TEM Staining', '')):
                        methods.append('TEM')
                    if clean_value(row.get('Th-T Binding', '')):
                        methods.append('ThT')
                    if clean_value(row.get('FTIR peaks', '')):
                        methods.append('FTIR')
                    if clean_value(row.get('Proteostat binding', '')):
                        methods.append('Proteostat')
                    
                    raw_method = ', '.join(methods) if methods else 'Experimental'
                    universal, evidence_type, weight, confidence = map_method_to_universal(raw_method)
                    
                    position = clean_value(row.get('Position', ''))
                    region_start = int(position) if position.isdigit() else None
                    region_end = region_start + len(seq) - 1 if region_start else None
                    
                    uniprot_id = extract_uniprot_id(row.get('UniProt ID', ''))
                    protein_name = clean_value(row.get('UniProt AC', ''))
                    
                    entry = AmyloidEntry(
                        record_id=f"WDB_{seq[:10]}_{idx}",
                        source_db="WALTZ-DB",
                        protein_name=protein_name,
                        uniprot_id=uniprot_id,
                        protein_family=self._infer_protein_family(protein_name),
                        region_start=region_start,
                        region_end=region_end,
                        sequence=seq,
                        is_amyloid=is_amyloid,
                        experimental_label=label,
                        category="hexapeptide",
                        structure_type=self._infer_structure_type(raw_method),
                        aggregate_type='synthetic',
                        pathogenicity='unknown',
                        raw_method=raw_method,
                        method_universal=universal,
                        evidence_type=evidence_type,
                        evidence_weight=weight,
                        confidence=confidence,
                        pdb_id=clean_value(row.get('Structures', '')),
                    )
                    
                    if entry.is_valid():
                        self.entries.append(entry)
                except Exception as e:
                    self.errors.append(f"Row {idx}: {str(e)}")
        
        self.log_result()
        return self.entries
