"""AmyloidAtlas parser."""

import csv
from typing import List

from .base import BaseParser
from ..models import AmyloidEntry
from ..utils.sequence import clean_value, parse_range, map_method_to_universal


class AmyloidAtlasParser(BaseParser):
    """Parser for AmyloidAtlas (TSV)."""
    
    def parse(self, filepath: str, **kwargs) -> List[AmyloidEntry]:
        self.logger.info(f"Parsing {self.db_name} from {filepath}...")
        
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for idx, row in enumerate(reader):
                residues = clean_value(row.get('Residues_Ordered', ''))
                range_list = parse_range(residues)
                
                raw_method = clean_value(row.get('Method', ''))
                universal, evidence_type, weight, confidence = map_method_to_universal(raw_method)
                
                protein_name = clean_value(row.get('Protein', ''))
                notes = clean_value(row.get('Fibril_Origins', ''))
                
                for start, end in range_list:
                    entry = AmyloidEntry(
                        record_id=clean_value(row.get('PDB_ID', '')),
                        source_db="AmyloidAtlas",
                        protein_name=protein_name,
                        protein_family=self._infer_protein_family(protein_name),
                        region_start=start,
                        region_end=end,
                        pdb_id=clean_value(row.get('PDB_ID', '')),
                        raw_method=raw_method,
                        method_universal=universal,
                        evidence_type=evidence_type,
                        evidence_weight=weight,
                        confidence=confidence,
                        resolution=clean_value(row.get('Resolution_A', '')),
                        doi=clean_value(row.get('DOI_URL', '')).replace('https://doi.org/', ''),
                        reference=clean_value(row.get('Reference', '')),
                        notes=notes,
                        is_amyloid=True,
                        experimental_label='amyloid',
                        structure_type=self._infer_structure_type(raw_method, notes),
                        aggregate_type=self._infer_aggregate_type(notes),
                        pathogenicity=self._infer_pathogenicity("", "", notes),
                    )
                    
                    if entry.is_valid():
                        self.entries.append(entry)
        
        self.log_result()
        return self.entries
