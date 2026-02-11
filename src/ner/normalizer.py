"""
Protein name normalization and UniProt enrichment
"""

import re
import logging
from typing import Dict, List, Optional, Tuple
import requests
from functools import lru_cache

logger = logging.getLogger(__name__)


class ProteinNormalizer:
    """Normalize protein names and resolve synonyms"""
    
    # Canonical names for common amyloid proteins
    CANONICAL_NAMES = {
        # Amyloid-beta variants
        'amyloid-beta': ['amyloid beta', 'amyloid-β', 'amyloid β', 'abeta', 'aβ', 
                         'a-beta', 'a-β', 'beta-amyloid', 'β-amyloid', 'beta amyloid',
                         'amyloid beta peptide', 'amyloid-beta peptide', 'abeta peptide',
                         'amyloid precursor protein', 'app'],
        'aβ40': ['abeta40', 'aβ-40', 'abeta-40', 'amyloid-beta-40', 'amyloid-beta(1-40)',
                 'aβ(1-40)', 'abeta(1-40)', 'aβ 40', 'abeta 40'],
        'aβ42': ['abeta42', 'aβ-42', 'abeta-42', 'amyloid-beta-42', 'amyloid-beta(1-42)',
                 'aβ(1-42)', 'abeta(1-42)', 'aβ 42', 'abeta 42'],
        
        # Alpha-synuclein
        'alpha-synuclein': ['α-synuclein', 'α synuclein', 'alpha synuclein', 'a-synuclein',
                           'αsyn', 'a-syn', 'asyn', 'snca', 'α-syn'],
        
        # Tau
        'tau': ['tau protein', 'mapt', 'tau-protein', 'microtubule-associated protein tau'],
        
        # Huntingtin
        'huntingtin': ['htt', 'mhtt', 'mutant huntingtin', 'huntingtin protein'],
        
        # Prion
        'prion': ['prp', 'prion protein', 'prpsc', 'prpc', 'prnp'],
        
        # TDP-43
        'tdp-43': ['tdp43', 'tardbp', 'tar dna-binding protein 43', 'tdp 43'],
        
        # SOD1
        'sod1': ['superoxide dismutase 1', 'cu/zn superoxide dismutase', 'sod-1'],
        
        # FUS
        'fus': ['fused in sarcoma', 'fus protein', 'tls'],
        
        # Transthyretin
        'transthyretin': ['ttr', 'prealbumin', 'attr', 'transthyretin amyloid'],
        
        # IAPP/Amylin
        'iapp': ['amylin', 'islet amyloid polypeptide', 'islet amyloid'],
        
        # Serum amyloid A
        'saa': ['serum amyloid a', 'saa1', 'saa2'],
        
        # Heat shock proteins
        'hsp70': ['heat shock protein 70', 'hspa1a', 'hsp 70', 'hsp-70'],
        'hsp90': ['heat shock protein 90', 'hsp 90', 'hsp-90'],
        
        # Other
        'lysozyme': ['lysozyme c', 'lyz'],
        'insulin': ['human insulin', 'ins'],
        'β2-microglobulin': ['b2m', 'beta-2-microglobulin', 'beta2-microglobulin', 'β2m'],
    }
    
    # Build reverse lookup
    _SYNONYM_TO_CANONICAL = {}
    for canonical, synonyms in CANONICAL_NAMES.items():
        _SYNONYM_TO_CANONICAL[canonical.lower()] = canonical
        for syn in synonyms:
            _SYNONYM_TO_CANONICAL[syn.lower()] = canonical
    
    @classmethod
    def normalize(cls, name: str) -> str:
        """Normalize protein name to canonical form"""
        if not name:
            return name
        
        # Clean whitespace and normalize
        cleaned = ' '.join(name.strip().split())
        lookup = cleaned.lower()
        
        # Check direct match
        if lookup in cls._SYNONYM_TO_CANONICAL:
            return cls._SYNONYM_TO_CANONICAL[lookup]
        
        # Check partial matches for long names
        for syn, canonical in cls._SYNONYM_TO_CANONICAL.items():
            if syn in lookup or lookup in syn:
                return canonical
        
        # Return cleaned lowercase if no match
        return cleaned.lower()
    
    @classmethod
    def are_same_protein(cls, name1: str, name2: str) -> bool:
        """Check if two names refer to the same protein"""
        return cls.normalize(name1) == cls.normalize(name2)
    
    @classmethod
    def deduplicate(cls, proteins: List[Dict]) -> List[Dict]:
        """Remove duplicate proteins, keeping best data"""
        seen = {}  # canonical_name -> best protein data
        
        for protein in proteins:
            name = protein.get('protein_name', '')
            canonical = cls.normalize(name)
            
            if canonical in seen:
                # Merge: keep higher confidence, combine data
                existing = seen[canonical]
                if protein.get('extraction_confidence', 0) > existing.get('extraction_confidence', 0):
                    # Update with better extraction but keep metadata
                    for key, val in protein.items():
                        if val and not existing.get(key):
                            existing[key] = val
                    existing['extraction_confidence'] = protein['extraction_confidence']
            else:
                protein['protein_name'] = canonical  # Use canonical name
                seen[canonical] = protein.copy()
        
        return list(seen.values())


class UniProtClient:
    """Fetch protein metadata from UniProt"""
    
    BASE_URL = "https://rest.uniprot.org"
    
    # Pre-defined mappings for common amyloid proteins (human)
    KNOWN_PROTEINS = {
        'amyloid-beta': {'uniprot_id': 'P05067', 'gene_name': 'APP', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'aβ40': {'uniprot_id': 'P05067', 'gene_name': 'APP', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'aβ42': {'uniprot_id': 'P05067', 'gene_name': 'APP', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'alpha-synuclein': {'uniprot_id': 'P37840', 'gene_name': 'SNCA', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'tau': {'uniprot_id': 'P10636', 'gene_name': 'MAPT', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'huntingtin': {'uniprot_id': 'P42858', 'gene_name': 'HTT', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'prion': {'uniprot_id': 'P04156', 'gene_name': 'PRNP', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'tdp-43': {'uniprot_id': 'Q13148', 'gene_name': 'TARDBP', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'sod1': {'uniprot_id': 'P00441', 'gene_name': 'SOD1', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'fus': {'uniprot_id': 'P35637', 'gene_name': 'FUS', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'transthyretin': {'uniprot_id': 'P02766', 'gene_name': 'TTR', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'iapp': {'uniprot_id': 'P10997', 'gene_name': 'IAPP', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'saa': {'uniprot_id': 'P0DJI8', 'gene_name': 'SAA1', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'hsp70': {'uniprot_id': 'P0DMV8', 'gene_name': 'HSPA1A', 'organism': 'Homo sapiens', 'forms_aggregates': False, 'functional_class': 'chaperone'},
        'hsp90': {'uniprot_id': 'P07900', 'gene_name': 'HSP90AA1', 'organism': 'Homo sapiens', 'forms_aggregates': False, 'functional_class': 'chaperone'},
        'lysozyme': {'uniprot_id': 'P61626', 'gene_name': 'LYZ', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
        'insulin': {'uniprot_id': 'P01308', 'gene_name': 'INS', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'functional'},
        'β2-microglobulin': {'uniprot_id': 'P61769', 'gene_name': 'B2M', 'organism': 'Homo sapiens', 'forms_aggregates': True, 'functional_class': 'pathological'},
    }
    
    @classmethod
    @lru_cache(maxsize=500)
    def lookup(cls, protein_name: str) -> Optional[Dict]:
        """Lookup protein in UniProt, returns metadata"""
        canonical = ProteinNormalizer.normalize(protein_name)
        
        # Check known proteins first
        if canonical in cls.KNOWN_PROTEINS:
            return cls.KNOWN_PROTEINS[canonical].copy()
        
        # Try UniProt API search
        try:
            query = f'(protein_name:"{protein_name}" OR gene:"{protein_name}") AND (organism_id:9606)'
            url = f"{cls.BASE_URL}/uniprotkb/search"
            params = {
                'query': query,
                'format': 'json',
                'size': 1,
                'fields': 'accession,gene_names,protein_name,organism_name'
            }
            
            response = requests.get(url, params=params, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if data.get('results'):
                    result = data['results'][0]
                    genes = result.get('genes', [])
                    gene_name = genes[0].get('geneName', {}).get('value') if genes else None
                    
                    return {
                        'uniprot_id': result.get('primaryAccession'),
                        'gene_name': gene_name,
                        'organism': result.get('organism', {}).get('scientificName'),
                    }
        except Exception as e:
            logger.debug(f"UniProt lookup failed for {protein_name}: {e}")
        
        return None
    
    @classmethod
    def enrich_proteins(cls, proteins: List[Dict]) -> List[Dict]:
        """Enrich protein list with UniProt metadata"""
        enriched = []
        for protein in proteins:
            name = protein.get('protein_name', '')
            metadata = cls.lookup(name)
            
            if metadata:
                # Only fill in missing fields
                for key, val in metadata.items():
                    if val and not protein.get(key):
                        protein[key] = val
            
            enriched.append(protein)
        
        return enriched


def normalize_and_enrich(proteins: List[Dict], use_uniprot: bool = True) -> List[Dict]:
    """Full normalization pipeline: normalize names, deduplicate, enrich"""
    # 1. Normalize and deduplicate
    deduped = ProteinNormalizer.deduplicate(proteins)
    
    # 2. Enrich with UniProt (optional)
    if use_uniprot:
        deduped = UniProtClient.enrich_proteins(deduped)
    
    return deduped


if __name__ == "__main__":
    # Test
    test_proteins = [
        {'protein_name': 'Amyloid-Beta', 'extraction_confidence': 0.8},
        {'protein_name': 'amyloid beta', 'extraction_confidence': 0.7},
        {'protein_name': 'Aβ42', 'extraction_confidence': 0.9},
        {'protein_name': 'Alpha-synuclein', 'extraction_confidence': 0.85},
        {'protein_name': 'α-synuclein', 'extraction_confidence': 0.75},
        {'protein_name': 'tau', 'extraction_confidence': 0.9},
        {'protein_name': 'HSP70', 'extraction_confidence': 0.6},
    ]
    
    print("Before normalization:")
    for p in test_proteins:
        print(f"  {p['protein_name']}")
    
    result = normalize_and_enrich(test_proteins, use_uniprot=True)
    
    print(f"\nAfter normalization ({len(result)} unique):")
    for p in result:
        print(f"  {p['protein_name']}: gene={p.get('gene_name')}, uniprot={p.get('uniprot_id')}")
