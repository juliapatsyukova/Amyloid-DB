#!/usr/bin/env python3
"""Command-line interface for amyloid pipeline."""

import argparse
import os

from .run_pipeline import run_pipeline
from .filters import AmyloidFilter
from .export import write_tsv


def main():
    parser = argparse.ArgumentParser(
        description='Unified Amyloid Database Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic usage with multiple databases
  python -m amyloid_pipeline --waltzdb waltzdb.tsv --crossbeta crossbeta.json -o output
  
  # Full pipeline with UniProt fetching
  python -m amyloid_pipeline --waltzdb waltzdb.tsv --amyloid-atlas atlas.tsv \\
      -o output --fetch-sequences -v
  
  # Filter results
  python -m amyloid_pipeline --waltzdb waltzdb.tsv -o output \\
      --exclude-computational --filter-organism "Homo sapiens"
'''
    )
    
    # Input databases
    inputs = parser.add_argument_group('Input databases')
    inputs.add_argument('--waltzdb', help='WALTZ-DB file (CSV/TSV)')
    inputs.add_argument('--crossbeta', help='Cross-Beta DB JSON')
    inputs.add_argument('--amyload', help='AmyLoad CSV')
    inputs.add_argument('--amyloid-explorer', help='AmyloidExplorer TSV')
    inputs.add_argument('--amyloid-atlas', help='AmyloidAtlas TSV')
    inputs.add_argument('--amylobase', help='Amylobase TSV')
    inputs.add_argument('--amypro', help='AmyPro TSV')
    inputs.add_argument('--amylograph', help='AmyloGraph CSV')
    inputs.add_argument('--cpad-peptides', help='CPAD peptides Excel')
    inputs.add_argument('--cpad-structures', help='CPAD structures Excel')
    
    # Output
    parser.add_argument('-o', '--output', default='amyloid_output', help='Output directory')
    
    # Options
    options = parser.add_argument_group('Options')
    options.add_argument('--fetch-sequences', action='store_true', help='Fetch from UniProt')
    options.add_argument('--no-features', action='store_true', help='Skip feature computation')
    options.add_argument('--no-sqlite', action='store_true', help='Skip SQLite export')
    options.add_argument('-v', '--verbose', action='store_true')
    
    # Filters
    filters = parser.add_argument_group('Filters')
    filters.add_argument('--filter-evidence', nargs='+', 
                        help='Include only these evidence types')
    filters.add_argument('--exclude-computational', action='store_true',
                        help='Exclude computational predictions')
    filters.add_argument('--filter-organism', help='Filter by organism (contains)')
    filters.add_argument('--filter-pathogenicity', nargs='+',
                        help='Filter by pathogenicity')
    filters.add_argument('--min-confidence', type=int, help='Minimum confidence')
    
    args = parser.parse_args()
    
    # Build input files dict
    input_files = {
        'waltzdb': args.waltzdb,
        'crossbeta': args.crossbeta,
        'amyload': args.amyload,
        'amyloid_explorer': args.amyloid_explorer,
        'amyloid_atlas': args.amyloid_atlas,
        'amylobase': args.amylobase,
        'amypro': args.amypro,
        'amylograph': args.amylograph,
        'cpad_peptides': args.cpad_peptides,
        'cpad_structures': args.cpad_structures,
    }
    input_files = {k: v for k, v in input_files.items() if v}
    
    if not input_files:
        parser.error("At least one input database file required")
    
    # Run pipeline
    results = run_pipeline(
        input_files=input_files,
        output_dir=args.output,
        fetch_sequences=args.fetch_sequences,
        compute_features=not args.no_features,
        export_sqlite=not args.no_sqlite,
        verbose=args.verbose
    )
    
    # Apply filters if specified
    if any([args.filter_evidence, args.filter_organism, args.filter_pathogenicity,
            args.min_confidence, args.exclude_computational]):
        
        print("\nApplying filters...")
        f = AmyloidFilter()
        
        if args.filter_evidence:
            f.evidence_type(args.filter_evidence)
        if args.exclude_computational:
            f.exclude_evidence_type('computational')
        if args.filter_organism:
            f.organism(args.filter_organism)
        if args.filter_pathogenicity:
            f.pathogenicity(args.filter_pathogenicity)
        if args.min_confidence:
            f.min_confidence(args.min_confidence)
        
        filtered = f.apply(results['consensus'])
        print(f"Filtered: {len(results['consensus'])} → {len(filtered)} entries")
        print(f"Filter: {f}")
        
        write_tsv(filtered, os.path.join(args.output, "filtered_entries.tsv"), include_features=True)


if __name__ == '__main__':
    main()
