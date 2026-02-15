"""
Configuration constants for the amyloid pipeline.
"""

# Evidence weighting system
EVIDENCE_WEIGHTS = {
    "structural": 3.0,
    "kinetic": 2.0,
    "staining_binding": 1.0,
    "literature_curated": 0.5,
    "computational": 0.0,
}

# Tier-based confidence (0-100)
TIER_CONFIDENCE = {
    0: 10,   # computational / unspecified
    1: 35,   # staining / chemical binding
    2: 50,   # biological tags, kinetics
    3: 70,   # microscopy, NMR, structural
    4: 90,   # XRD / fiber diffraction
    5: 80,   # spectroscopy (FTIR/CD/Raman)
}

# Valid categories
STRUCTURE_TYPES = ['fibril', 'oligomer', 'crystal', 'aggregate', 'monomer', 'unknown']
AGGREGATE_TYPES = ['in_vitro', 'ex_vivo', 'patient_derived', 'recombinant', 'synthetic', 'unknown']
PATHOGENICITY_TYPES = ['pathogenic', 'functional', 'non_pathogenic', 'unknown']

# Standard amino acids
STANDARD_AA = set('ACDEFGHIKLMNPQRSTVWY')

# API endpoints
UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Request settings
REQUEST_DELAY = 0.4  # seconds between API requests
REQUEST_TIMEOUT = 30  # seconds

# Methods mapping rules
METHODS_RULES = [
    # XRD / diffraction (tier 4)
    {"pattern": r"\bxrd\b|\bx[-\s]?ray diffraction\b|\bfiber diffraction\b",
     "universal": "XRD (cross-β diffraction)", "tier": 4, "evidence_type": "structural"},
    
    # Microscopy (tier 3)
    {"pattern": r"\bcryo[-\s]?em\b", "universal": "Cryo-EM", "tier": 3, "evidence_type": "structural"},
    {"pattern": r"\btem\b|\btransmission electron microscopy\b|\belectron microscopy\b",
     "universal": "Electron microscopy (TEM/EM)", "tier": 3, "evidence_type": "structural"},
    {"pattern": r"\bafm\b|\batomic force microscopy\b", "universal": "AFM", "tier": 3, "evidence_type": "structural"},
    
    # NMR (tier 3)
    {"pattern": r"\bsolid[-\s]?state nmr\b|\bssnmr\b", "universal": "NMR (solid-state)", "tier": 3, "evidence_type": "structural"},
    {"pattern": r"\bsolution nmr\b", "universal": "NMR (solution)", "tier": 2, "evidence_type": "kinetic"},
    {"pattern": r"\bnmr\b", "universal": "NMR (unspecified)", "tier": 3, "evidence_type": "structural"},
    
    # X-ray crystallography (tier 3)
    {"pattern": r"\bcrystallograph\b|\bx-ray\b|\bxray\b", "universal": "X-ray (structural)", "tier": 3, "evidence_type": "structural"},
    {"pattern": r"\bmicroed\b|\belectron crystallography\b", "universal": "Electron crystallography / microED", "tier": 3, "evidence_type": "structural"},
    
    # Spectroscopy (tier 5)
    {"pattern": r"\bftir\b|\binfrared\b", "universal": "FTIR", "tier": 5, "evidence_type": "staining_binding"},
    {"pattern": r"\bcircular dichroism\b|\bcd\b", "universal": "CD", "tier": 5, "evidence_type": "staining_binding"},
    {"pattern": r"\braman\b", "universal": "Raman", "tier": 5, "evidence_type": "staining_binding"},
    
    # Staining / binding (tier 1)
    {"pattern": r"\bthioflavin[-\s]?t\b|\btht\b|\bth[-\s]?t\b", "universal": "ThT binding", "tier": 1, "evidence_type": "staining_binding"},
    {"pattern": r"\bcongo[-\s]?red\b", "universal": "Congo Red binding", "tier": 1, "evidence_type": "staining_binding"},
    {"pattern": r"\bproteostat\b", "universal": "Proteostat binding", "tier": 1, "evidence_type": "staining_binding"},
    
    # Kinetics (tier 2)
    {"pattern": r"\bkinetic\b|\blag\b|\baggregation rate\b|\bseeding\b|\baggregation assay\b",
     "universal": "Aggregation kinetics", "tier": 2, "evidence_type": "kinetic"},
    
    # Computational (tier 0)
    {"pattern": r"\bpasta\b", "universal": "PASTA 2.0 (prediction)", "tier": 0, "evidence_type": "computational"},
    {"pattern": r"\baggrescan\b", "universal": "AGGRESCAN (prediction)", "tier": 0, "evidence_type": "computational"},
    {"pattern": r"\btango\b", "universal": "TANGO (prediction)", "tier": 0, "evidence_type": "computational"},
    {"pattern": r"\bwaltz\b", "universal": "WALTZ (prediction)", "tier": 0, "evidence_type": "computational"},
    
    # Literature (tier 0)
    {"pattern": r"\bliterature[-\s]?curated\b|\bcurated\b", "universal": "Literature-curated", "tier": 0, "evidence_type": "literature_curated"},
]

# Protein family keywords
PROTEIN_FAMILIES = {
    'immunoglobulin': ['immunoglobulin', 'antibody', 'igg', 'light chain', 'kappa', 'lambda'],
    'amyloid-beta': ['amyloid beta', 'abeta', 'a-beta', 'aβ', 'app'],
    'tau': ['tau', 'mapt'],
    'alpha-synuclein': ['synuclein', 'alpha-syn', 'snca'],
    'prion': ['prion', 'prnp', 'prp'],
    'transthyretin': ['transthyretin', 'ttr'],
    'serum amyloid A': ['serum amyloid a', 'saa'],
    'lysozyme': ['lysozyme'],
    'insulin': ['insulin'],
    'islet amyloid': ['islet amyloid', 'iapp', 'amylin'],
    'huntingtin': ['huntingtin', 'htt', 'polyq'],
    'sod1': ['sod1', 'superoxide dismutase'],
    'tdp-43': ['tdp-43', 'tdp43', 'tardbp'],
    'fus': ['fus'],
    'apolipoprotein': ['apolipoprotein', 'apoa', 'apoe'],
}
