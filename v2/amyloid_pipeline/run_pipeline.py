"""Main pipeline execution."""

import json
import logging
import os
import sys
from collections import defaultdict
from typing import Dict, List

from .models import AmyloidEntry
from .parsers import PARSERS
from .unifier import DatasetUnifier
from .utils.fetchers import SequenceCache, fetch_sequences_and_info
from .export import write_tsv, write_fasta, create_sqlite_database


def print_statistics(entries: List[AmyloidEntry], unifier: DatasetUnifier):
    """Print pipeline statistics."""
    print("\n" + "="*70)
    print("UNIFIED AMYLOID DATABASE STATISTICS")
    print("="*70)
    
    by_source = defaultdict(int)
    for e in entries:
        by_source[e.source_db] += 1
    print("\nEntries by source:")
    for db, count in sorted(by_source.items(), key=lambda x: -x[1]):
        print(f"  {db}: {count}")
    
    amyloid = sum(1 for e in entries if e.is_amyloid)
    print(f"\nClassification: {amyloid} amyloid / {len(entries) - amyloid} non-amyloid")
    
    by_evidence = defaultdict(int)
    for e in entries:
        by_evidence[e.evidence_type] += 1
    print(f"\nBy evidence type:")
    for ev, count in sorted(by_evidence.items(), key=lambda x: -x[1]):
        print(f"  {ev}: {count}")
    
    by_structure = defaultdict(int)
    for e in entries:
        by_structure[e.structure_type] += 1
    print(f"\nBy structure type:")
    for st, count in sorted(by_structure.items(), key=lambda x: -x[1]):
        print(f"  {st}: {count}")
    
    by_pathogenicity = defaultdict(int)
    for e in entries:
        by_pathogenicity[e.pathogenicity] += 1
    print(f"\nBy pathogenicity:")
    for p, count in sorted(by_pathogenicity.items(), key=lambda x: -x[1]):
        print(f"  {p}: {count}")
    
    with_seq = sum(1 for e in entries if e.sequence)
    with_uniprot = sum(1 for e in entries if e.uniprot_id)
    with_family = sum(1 for e in entries if e.protein_family)
    print(f"\nData completeness:")
    print(f"  With sequence: {with_seq}")
    print(f"  With UniProt ID: {with_uniprot}")
    print(f"  With protein family: {with_family}")
    
    stats = unifier.get_stats()
    print(f"\nDeduplication: {stats['total_entries']} → {stats['unique_keys']} unique")
    print("="*70)


def run_pipeline(input_files: Dict[str, str], output_dir: str,
                 fetch_sequences: bool = False, compute_features: bool = True,
                 export_sqlite: bool = True, verbose: bool = False) -> Dict:
    """
    Run the complete amyloid database pipeline.
    
    Args:
        input_files: Dict mapping database names to file paths
        output_dir: Output directory
        fetch_sequences: Whether to fetch sequences from UniProt
        compute_features: Whether to compute physicochemical features
        export_sqlite: Whether to export SQLite database
        verbose: Enable verbose logging
    
    Returns:
        Dict with results (all_entries, consensus, amyloid, non_amyloid, conflicts)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Setup logging
    log_file = os.path.join(output_dir, "pipeline.log")
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='[%(asctime)s] [%(levelname)s] %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)]
    )
    logger = logging.getLogger(__name__)
    
    logger.info("="*70)
    logger.info("UNIFIED AMYLOID DATABASE PIPELINE")
    logger.info("="*70)
    
    # Phase 1: Parse databases
    logger.info("\nPHASE 1: PARSING DATABASES")
    
    unifier = DatasetUnifier(logger)
    
    for db_key, filepath in input_files.items():
        if not filepath or not os.path.exists(filepath):
            logger.warning(f"Skipping {db_key}: file not found")
            continue
        
        db_key_norm = db_key.lower().replace('-', '_').replace(' ', '_')
        parser_class = PARSERS.get(db_key_norm)
        
        if parser_class:
            parser = parser_class(db_key, logger)
            entries = parser.parse(filepath)
            unifier.add_entries(entries, db_key)
        else:
            logger.warning(f"Unknown database: {db_key}")
    
    # Phase 2: Fetch sequences
    if fetch_sequences:
        logger.info("\nPHASE 2: FETCHING SEQUENCES")
        cache_path = os.path.join(output_dir, "sequence_cache.json")
        cache = SequenceCache(cache_path)
        fetch_sequences_and_info(unifier.all_entries, cache, logger)
    
    # Phase 3: Deduplication
    logger.info("\nPHASE 3: DEDUPLICATION & CONSENSUS")
    conflicts = unifier.detect_conflicts()
    consensus = unifier.build_consensus()
    
    # Phase 4: Compute features
    if compute_features:
        logger.info("\nPHASE 4: COMPUTING FEATURES")
        for entry in consensus:
            if entry.sequence:
                entry.compute_features()
        logger.info(f"Computed features for {sum(1 for e in consensus if e._physicochemical)} entries")
    
    # Phase 5: Export
    logger.info("\nPHASE 5: EXPORTING")
    
    amyloid_entries = [e for e in consensus if e.is_amyloid]
    non_amyloid_entries = [e for e in consensus if not e.is_amyloid]
    
    write_tsv(consensus, os.path.join(output_dir, "consensus_unified.tsv"), include_features=compute_features)
    write_tsv(amyloid_entries, os.path.join(output_dir, "amyloid_positive.tsv"), include_features=compute_features)
    write_tsv(non_amyloid_entries, os.path.join(output_dir, "non_amyloid.tsv"), include_features=compute_features)
    
    write_fasta(amyloid_entries, os.path.join(output_dir, "amyloid_sequences.fasta"))
    write_fasta(non_amyloid_entries, os.path.join(output_dir, "non_amyloid_sequences.fasta"))
    
    if export_sqlite:
        create_sqlite_database(
            consensus,
            os.path.join(output_dir, "amyloid_database.sqlite"),
            include_features=compute_features,
            logger=logger
        )
    
    if conflicts:
        with open(os.path.join(output_dir, "conflicts.json"), 'w') as f:
            conflict_data = [{
                'key': str(c['key']),
                'labels': list(c['labels']),
                'sources': c['sources'],
            } for c in conflicts]
            json.dump(conflict_data, f, indent=2)
    
    print_statistics(consensus, unifier)
    
    logger.info("\n✅ Pipeline completed!")
    logger.info(f"Output: {output_dir}")
    
    return {
        'all_entries': unifier.all_entries,
        'consensus': consensus,
        'amyloid': amyloid_entries,
        'non_amyloid': non_amyloid_entries,
        'conflicts': conflicts,
    }
