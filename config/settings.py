"""
Configuration settings for Amyloid Aggregation Database Pipeline
"""

# =============================================================================
# SEARCH TERMS
# =============================================================================

# Primary aggregation terms (broad definition)
AGGREGATION_TERMS = [
    "amyloid",
    "amyloid-like", 
    "amyloidogenic",
    "protein aggregation",
    "protein aggregates",
    "fibril formation",
    "fibrillar aggregates",
    "cross-beta",
    "cross-beta structure",
    "beta-sheet aggregation",
    "alpha-helical amyloid",  # newer structural diversity
    "functional amyloid",
    "pathological amyloid",
    "self-assembly protein",
    "protein self-assembly",
]

# Exclusion terms (prions and prion-like)
EXCLUSION_TERMS = [
    "prion",
    "prion-like",
    "prion disease",
    "prionoid",
    "transmissible spongiform",
    "PrP",
    "PrPSc",
    "PrPC",
]

# Experimental validation methods
EXPERIMENTAL_METHODS = {
    "biophysical": [
        "ThT", "Thioflavin T", "Thioflavin-T",
        "ThS", "Thioflavin S",
        "Congo red", "Congo Red",
        "circular dichroism", "CD spectroscopy",
        "FTIR", "infrared spectroscopy",
        "X-ray diffraction", "X-ray fiber diffraction",
        "solid-state NMR", "ssNMR",
    ],
    "microscopy": [
        "electron microscopy", "EM", "TEM", "transmission electron",
        "cryo-EM", "cryo-electron microscopy",
        "AFM", "atomic force microscopy",
        "super-resolution microscopy",
    ],
    "structural": [
        "crystal structure",
        "cryo-EM structure",
        "NMR structure",
        "fibril structure",
    ],
    "cellular": [
        "in vivo aggregation",
        "cell culture",
        "inclusion bodies",
        "puncta formation",
    ],
    "kinetic": [
        "aggregation kinetics",
        "seeding assay",
        "nucleation",
        "elongation rate",
        "lag phase",
    ]
}

# Assembly/disassembly indicators
ASSEMBLY_PATTERNS = {
    "reversible": [
        "reversible aggregation",
        "dynamic assembly",
        "regulated assembly",
        "disassembly",
        "disaggregation",
        "chaperone-mediated",
        "ATP-dependent",
    ],
    "irreversible": [
        "irreversible aggregation",
        "stable aggregates",
        "persistent aggregates",
        "insoluble aggregates",
        "resistant to",
    ],
    "autocatalytic": [
        "self-seeding",
        "autocatalytic",
        "template-assisted",
        "seeding-competent",
        "propagation",
    ],
    "cofactor_dependent": [
        "RNA-mediated",
        "chaperone-assisted",
        "cofactor",
        "metal-induced",
        "lipid-induced",
        "membrane-mediated",
    ]
}

# Functional vs Pathological context
FUNCTIONAL_CONTEXT = [
    "functional amyloid",
    "biofilm",
    "secretory granule",
    "hormone storage",
    "melanosome",
    "bacterial",
    "curli",
    "Pmel17",
    "CPEB",
]

PATHOLOGICAL_CONTEXT = [
    "neurodegenerative",
    "Alzheimer",
    "Parkinson",
    "Huntington",
    "ALS",
    "systemic amyloidosis",
    "AL amyloidosis",
    "AA amyloidosis",
    "ATTR",
    "pathological",
    "toxic",
    "disease-associated",
]

# =============================================================================
# API SETTINGS
# =============================================================================

PUBMED_API = {
    "base_url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
    "db": "pubmed",
    "retmax": 1000,  # max results per query
    "retmode": "xml",
    "email": "your_email@example.com",  # required for NCBI API
    "api_key": None,  # optional, increases rate limit
    "rate_limit_delay": 0.34,  # seconds between requests (3/sec without key)
}

STRING_DB_API = {
    "base_url": "https://string-db.org/api",
    "version": "12.0",
    "output_format": "json",
    "species": 9606,  # human by default
    "required_score": 700,  # high confidence
    "network_type": "physical",  # physical subnetwork
}

# LLM Settings
LLM_CONFIG = {
    "provider": "anthropic",  # or "openai"
    "model": "claude-sonnet-4-20250514",
    "max_tokens": 4096,
    "temperature": 0.1,  # low for structured extraction
}

# =============================================================================
# DATABASE SETTINGS
# =============================================================================

DATABASE = {
    "type": "sqlite",  # start simple, can migrate to PostgreSQL
    "path": "data/amyloid_aggregation.db",
}

# =============================================================================
# EXTRACTION SCHEMA
# =============================================================================

# What we want to extract from each paper
EXTRACTION_SCHEMA = {
    "protein": {
        "name": str,
        "uniprot_id": str,
        "gene_name": str,
        "organism": str,
    },
    "aggregation_properties": {
        "forms_aggregates": bool,
        "aggregate_structure": str,  # beta-sheet, alpha-helical, mixed, unknown
        "reversibility": str,  # reversible, irreversible, unknown
        "assembly_mechanism": str,  # autocatalytic, cofactor-dependent, unknown
    },
    "functional_classification": {
        "is_functional": bool,
        "is_pathological": bool,
        "biological_role": str,
    },
    "experimental_evidence": {
        "methods_used": list,
        "confidence_level": str,  # high, medium, low
    },
    "cofactors": {
        "assembly_helpers": list,  # proteins, RNAs, small molecules
        "disassembly_helpers": list,
    },
    "source": {
        "pmid": str,
        "doi": str,
        "title": str,
        "year": int,
    }
}
