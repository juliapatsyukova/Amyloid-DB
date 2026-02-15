"""
Unified Amyloid Database Pipeline
==================================

Modular pipeline for parsing, filtering, and exporting amyloid protein databases.

Usage:
    from amyloid_pipeline import run_pipeline, AmyloidFilter
    
    results = run_pipeline(
        input_files={'waltzdb': 'waltzdb.tsv'},
        output_dir='output'
    )
    
    filtered = AmyloidFilter().exclude_evidence_type('computational').apply(results['consensus'])
"""

from .models import AmyloidEntry, PhysicochemicalFeatures, SequenceComposition
from .filters import AmyloidFilter, filter_entries
from .unifier import DatasetUnifier
from .run_pipeline import run_pipeline

__version__ = "2.0.0"
__all__ = [
    'AmyloidEntry',
    'PhysicochemicalFeatures', 
    'SequenceComposition',
    'AmyloidFilter',
    'filter_entries',
    'DatasetUnifier',
    'run_pipeline',
]
