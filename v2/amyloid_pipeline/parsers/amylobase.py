"""Amylobase parser."""

import csv
from typing import List

from .base import BaseParser
from ..models import AmyloidEntry
from ..config import EVIDENCE_WEIGHTS, TIER_CONFIDENCE
from ..utils.sequence import clean_value, extract_uniprot_id


class AmylobaseParser(BaseParser):
    """Parser for Amylobase (TSV)."""
    
    def parse(self, filepath: str, **kwargs) -> List[AmyloidEntry]:
        self.logger.info(f"Parsing {self.db_name} from {filepath}...")
        
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for idx, row in enumerate(reader):
                try:
                    start = int(row.get('Start', 0)) if row.get('Start') else None
                    end = int(row.get('Stop', 0)) if row.get('Stop') else None
                except:
                    start, end = None, None
                
                mutation = clean_value(row.get('Mutation', ''))
                if mutation.lower() == 'na':
                    mutation = ""
                
                protein_name = clean_value(row.get('Protein name', ''))
                
                entry = AmyloidEntry(
                    record_id=f"AB_{extract_uniprot_id(row.get('Protein ID', ''))}_{idx}",
                    source_db="Amylobase",
                    protein_name=protein_name,
                    uniprot_id=extract_uniprot_id(row.get('Protein ID', '')),
                    organism=clean_value(row.get('Organism', '')),
                    protein_family=self._infer_protein_family(protein_name),
                    region_start=start,
                    region_end=end,
                    mutation=mutation,
                    category="kinetic_data",
                    is_amyloid=True,
                    experimental_label='amyloid',
                    structure_type='aggregate',
                    aggregate_type='in_vitro',
                    pathogenicity='unknown',
                    raw_method="Aggregation kinetics",
                    method_universal="Aggregation kinetics",
                    evidence_type="kinetic",
                    evidence_weight=EVIDENCE_WEIGHTS["kinetic"],
                    confidence=TIER_CONFIDENCE[2],
                    notes=f"pH={row.get('pH', '')}, T={row.get('Temperature (Kelvin Degrees)', '')}K"
                )
                
                if entry.is_valid():
                    self.entries.append(entry)
        
        self.log_result()
        return self.entries
