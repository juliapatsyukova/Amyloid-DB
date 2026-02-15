"""AmyPro parser."""

import csv
from typing import List

from .base import BaseParser
from ..models import AmyloidEntry
from ..config import EVIDENCE_WEIGHTS, TIER_CONFIDENCE
from ..utils.sequence import clean_value, extract_uniprot_id, parse_range


class AmyProParser(BaseParser):
    """Parser for AmyPro (TSV)."""
    
    def parse(self, filepath: str, **kwargs) -> List[AmyloidEntry]:
        self.logger.info(f"Parsing {self.db_name} from {filepath}...")
        
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for idx, row in enumerate(reader):
                regions_str = clean_value(row.get('regions', ''))
                range_list = parse_range(regions_str)
                
                category = clean_value(row.get('category', ''))
                is_prion = row.get('prion', '').lower() == 'true'
                if is_prion:
                    category = f"prion/{category}" if category else "prion"
                
                seq = clean_value(row.get('sequence', ''))
                protein_name = clean_value(row.get('protein', ''))
                
                pathogenicity = 'functional' if 'functional' in category.lower() else 'unknown'
                if is_prion:
                    pathogenicity = 'pathogenic'
                
                for start, end in range_list:
                    entry = AmyloidEntry(
                        record_id=clean_value(row.get('ID', '')).replace('#', ''),
                        source_db="AmyPro",
                        protein_name=protein_name,
                        uniprot_id=extract_uniprot_id(row.get('uniprot', '')),
                        organism=clean_value(row.get('species', '')),
                        protein_family=self._infer_protein_family(protein_name),
                        region_start=start,
                        region_end=end,
                        sequence=seq,
                        full_protein_sequence=seq,
                        pdb_id=clean_value(row.get('pdb', '')),
                        pmid=clean_value(row.get('pubmed', '')),
                        category=category,
                        mutation=clean_value(row.get('mutations', '')),
                        is_amyloid=True,
                        experimental_label='amyloid',
                        structure_type='fibril' if is_prion else 'unknown',
                        aggregate_type='unknown',
                        pathogenicity=pathogenicity,
                        raw_method="Literature-curated",
                        method_universal="Literature-curated",
                        evidence_type="literature_curated",
                        evidence_weight=EVIDENCE_WEIGHTS["literature_curated"],
                        confidence=TIER_CONFIDENCE[0]
                    )
                    
                    if entry.is_valid():
                        self.entries.append(entry)
        
        self.log_result()
        return self.entries
