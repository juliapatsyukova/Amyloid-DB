"""Base parser class."""

import logging
from typing import List

from ..models import AmyloidEntry
from ..config import PROTEIN_FAMILIES


class BaseParser:
    """Base class for database parsers."""
    
    def __init__(self, db_name: str, logger=None):
        self.db_name = db_name
        self.logger = logger or logging.getLogger(__name__)
        self.entries: List[AmyloidEntry] = []
        self.errors: List[str] = []
    
    def log_result(self):
        valid = sum(1 for e in self.entries if e.is_valid())
        self.logger.info(f"{self.db_name}: {len(self.entries)} entries ({valid} valid), {len(self.errors)} errors")
    
    def parse(self, filepath: str, **kwargs) -> List[AmyloidEntry]:
        raise NotImplementedError
    
    def _infer_structure_type(self, method: str, notes: str = "") -> str:
        text = f"{method} {notes}".lower()
        if 'fibril' in text or 'fiber' in text:
            return 'fibril'
        if 'oligomer' in text:
            return 'oligomer'
        if 'crystal' in text:
            return 'crystal'
        if 'aggregate' in text:
            return 'aggregate'
        if 'monomer' in text:
            return 'monomer'
        return 'unknown'
    
    def _infer_aggregate_type(self, notes: str, patient: bool = False) -> str:
        if patient:
            return 'patient_derived'
        text = notes.lower()
        if 'patient' in text or 'ex vivo' in text or 'brain' in text:
            return 'patient_derived'
        if 'ex vivo' in text:
            return 'ex_vivo'
        if 'recombinant' in text:
            return 'recombinant'
        if 'synthetic' in text or 'peptide' in text:
            return 'synthetic'
        if 'in vitro' in text:
            return 'in_vitro'
        return 'unknown'
    
    def _infer_pathogenicity(self, disease: str, category: str, notes: str = "") -> str:
        text = f"{disease} {category} {notes}".lower()
        if any(d in text for d in ['alzheimer', 'parkinson', 'huntington', 'als', 'prion', 'diabetes', 'amyloidosis']):
            return 'pathogenic'
        if 'functional' in text or 'biofilm' in text:
            return 'functional'
        if 'non-pathogenic' in text or 'non pathogenic' in text:
            return 'non_pathogenic'
        if disease:
            return 'pathogenic'
        return 'unknown'
    
    def _infer_protein_family(self, protein_name: str, disease: str = "") -> str:
        text = f"{protein_name} {disease}".lower()
        for family, keywords in PROTEIN_FAMILIES.items():
            if any(kw in text for kw in keywords):
                return family
        return ""
