"""Cross-Beta DB parser."""

import json
from typing import List

from .base import BaseParser
from ..models import AmyloidEntry
from ..utils.sequence import clean_value, validate_sequence, extract_uniprot_id, map_method_to_universal


class CrossBetaDBParser(BaseParser):
    """Parser for Cross-Beta DB (JSON)."""
    
    def parse(self, filepath: str, **kwargs) -> List[AmyloidEntry]:
        self.logger.info(f"Parsing {self.db_name} from {filepath}...")
        
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        for idx, record in enumerate(data):
            try:
                seq = clean_value(record.get('AR Sequence', '')).upper()
                if not validate_sequence(seq):
                    continue
                
                label_raw = clean_value(record.get('LABEL', '')).lower()
                if 'non' in label_raw:
                    label, is_amyloid = 'non-amyloid', False
                elif 'amyloid' in label_raw:
                    label, is_amyloid = 'amyloid', True
                else:
                    continue
                
                raw_method = clean_value(record.get('Method used', 'Experimental'))
                universal, evidence_type, weight, confidence = map_method_to_universal(raw_method)
                
                region_str = clean_value(record.get('Experimental Amyloid Region', ''))
                region_start, region_end = None, None
                if region_str and '-' in region_str:
                    parts = region_str.split('-')
                    try:
                        region_start, region_end = int(parts[0]), int(parts[1])
                    except:
                        pass
                
                organism = clean_value(record.get('Species / organism', ''))
                protein_name = clean_value(record.get('AR containing protein Accession Code(s)', ''))
                
                entry = AmyloidEntry(
                    record_id=f"CB_{idx}",
                    source_db="Cross-Beta DB",
                    protein_name=protein_name,
                    uniprot_id=extract_uniprot_id(protein_name),
                    organism=organism,
                    protein_family=self._infer_protein_family(protein_name),
                    region_start=region_start,
                    region_end=region_end,
                    sequence=seq,
                    is_amyloid=is_amyloid,
                    experimental_label=label,
                    structure_type=self._infer_structure_type(raw_method),
                    aggregate_type=self._infer_aggregate_type(""),
                    pathogenicity='unknown',
                    raw_method=raw_method,
                    method_universal=universal,
                    evidence_type=evidence_type,
                    evidence_weight=weight,
                    confidence=confidence,
                    doi=clean_value(record.get('Reference link / DOI', ''))
                )
                
                if entry.is_valid():
                    self.entries.append(entry)
            except Exception as e:
                self.errors.append(f"Record {idx}: {str(e)}")
        
        self.log_result()
        return self.entries
