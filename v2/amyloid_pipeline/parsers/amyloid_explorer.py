"""AmyloidExplorer parser."""

import csv
from typing import List

from .base import BaseParser
from ..models import AmyloidEntry
from ..utils.sequence import clean_value, extract_uniprot_id, map_method_to_universal


class AmyloidExplorerParser(BaseParser):
    """Parser for AmyloidExplorer (TSV)."""
    
    def parse(self, filepath: str, **kwargs) -> List[AmyloidEntry]:
        self.logger.info(f"Parsing {self.db_name} from {filepath}...")
        
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for idx, row in enumerate(reader):
                try:
                    start = int(row.get('resStart', 0)) if row.get('resStart') else None
                    end = int(row.get('resEnd', 0)) if row.get('resEnd') else None
                except:
                    start, end = None, None
                
                raw_method = clean_value(row.get('Method', ''))
                universal, evidence_type, weight, confidence = map_method_to_universal(raw_method)
                
                disease = clean_value(row.get('Disease', ''))
                protein_name = clean_value(row.get('Protein', ''))
                is_patient = row.get('Patient') == 'True'
                notes = clean_value(row.get('Fibril_Origins', '') or row.get('Notes', ''))
                
                entry = AmyloidEntry(
                    record_id=clean_value(row.get('Name', '')),
                    source_db="AmyloidExplorer",
                    protein_name=protein_name,
                    uniprot_id=extract_uniprot_id(row.get('UniprotId', '')),
                    protein_family=self._infer_protein_family(protein_name, disease),
                    region_start=start,
                    region_end=end,
                    pdb_id=clean_value(row.get('Name', '')),
                    emdb_id=clean_value(row.get('EMDB', '')),
                    raw_method=raw_method,
                    method_universal=universal,
                    evidence_type=evidence_type,
                    evidence_weight=weight,
                    confidence=confidence,
                    resolution=clean_value(row.get('Resolution', '')),
                    disease=disease,
                    tissue=clean_value(row.get('Tissue', '')),
                    mutation=clean_value(row.get('Mutant', '')),
                    doi=clean_value(row.get('DOI', '')),
                    pmid=clean_value(row.get('PMID', '')),
                    is_amyloid=True,
                    experimental_label='amyloid',
                    category="pathogenic" if is_patient else "in_vitro",
                    structure_type=self._infer_structure_type(raw_method, notes),
                    aggregate_type=self._infer_aggregate_type(notes, is_patient),
                    pathogenicity=self._infer_pathogenicity(disease, "", notes),
                )
                
                if entry.is_valid():
                    self.entries.append(entry)
        
        self.log_result()
        return self.entries
