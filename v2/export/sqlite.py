"""SQLite database export."""

import sqlite3
from typing import List

from ..models import AmyloidEntry


def create_sqlite_database(entries: List[AmyloidEntry], db_path: str, 
                           include_features: bool = True, logger=None):
    """Export entries to SQLite database with separate tables for features."""
    import logging
    if logger is None:
        logger = logging.getLogger(__name__)
    
    logger.info(f"Creating SQLite database: {db_path}")
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Main entries table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS entries (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            record_id TEXT,
            source_db TEXT,
            protein_name TEXT,
            uniprot_id TEXT,
            organism TEXT,
            protein_family TEXT,
            region_start INTEGER,
            region_end INTEGER,
            sequence TEXT,
            sequence_length INTEGER,
            is_amyloid BOOLEAN,
            experimental_label TEXT,
            category TEXT,
            structure_type TEXT,
            aggregate_type TEXT,
            pathogenicity TEXT,
            pdb_id TEXT,
            emdb_id TEXT,
            raw_method TEXT,
            method_universal TEXT,
            evidence_type TEXT,
            evidence_weight REAL,
            confidence INTEGER,
            resolution TEXT,
            disease TEXT,
            tissue TEXT,
            mutation TEXT,
            doi TEXT,
            pmid TEXT,
            reference TEXT,
            notes TEXT,
            UNIQUE(record_id, source_db, sequence)
        )
    ''')
    
    if include_features:
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS physicochemical_features (
                entry_id INTEGER PRIMARY KEY,
                length INTEGER,
                molecular_weight REAL,
                hydrophobicity_mean REAL,
                hydrophobicity_std REAL,
                hydrophobicity_max REAL,
                hydrophobicity_min REAL,
                net_charge REAL,
                positive_residues INTEGER,
                negative_residues INTEGER,
                charge_density REAL,
                beta_propensity_mean REAL,
                beta_propensity_max REAL,
                aggregation_propensity_mean REAL,
                aggregation_propensity_max REAL,
                aromatic_fraction REAL,
                aliphatic_fraction REAL,
                polar_fraction REAL,
                charged_fraction REAL,
                has_polyq BOOLEAN,
                has_glycine_rich BOOLEAN,
                has_proline BOOLEAN,
                FOREIGN KEY (entry_id) REFERENCES entries(id)
            )
        ''')
        
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS sequence_composition (
                entry_id INTEGER PRIMARY KEY,
                aa_A REAL, aa_C REAL, aa_D REAL, aa_E REAL, aa_F REAL,
                aa_G REAL, aa_H REAL, aa_I REAL, aa_K REAL, aa_L REAL,
                aa_M REAL, aa_N REAL, aa_P REAL, aa_Q REAL, aa_R REAL,
                aa_S REAL, aa_T REAL, aa_V REAL, aa_W REAL, aa_Y REAL,
                tiny_fraction REAL,
                small_fraction REAL,
                large_fraction REAL,
                FOREIGN KEY (entry_id) REFERENCES entries(id)
            )
        ''')
    
    # Indices
    for col in ['uniprot_id', 'organism', 'is_amyloid', 'evidence_type', 
                'structure_type', 'pathogenicity', 'protein_family']:
        cursor.execute(f'CREATE INDEX IF NOT EXISTS idx_{col} ON entries({col})')
    
    # Insert data
    for entry in entries:
        if include_features and entry.sequence:
            entry.compute_features()
        
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO entries (
                    record_id, source_db, protein_name, uniprot_id, organism, protein_family,
                    region_start, region_end, sequence, sequence_length, is_amyloid, 
                    experimental_label, category, structure_type, aggregate_type, pathogenicity,
                    pdb_id, emdb_id, raw_method, method_universal, evidence_type, 
                    evidence_weight, confidence, resolution, disease, tissue, mutation,
                    doi, pmid, reference, notes
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                entry.record_id, entry.source_db, entry.protein_name, entry.uniprot_id,
                entry.organism, entry.protein_family, entry.region_start, entry.region_end,
                entry.sequence, len(entry.sequence) if entry.sequence else 0, entry.is_amyloid,
                entry.experimental_label, entry.category, entry.structure_type,
                entry.aggregate_type, entry.pathogenicity, entry.pdb_id, entry.emdb_id,
                entry.raw_method, entry.method_universal, entry.evidence_type,
                entry.evidence_weight, entry.confidence, entry.resolution, entry.disease,
                entry.tissue, entry.mutation, entry.doi, entry.pmid, entry.reference, entry.notes
            ))
            
            entry_id = cursor.lastrowid
            
            if include_features and entry._physicochemical:
                phys = entry._physicochemical
                cursor.execute('''
                    INSERT OR REPLACE INTO physicochemical_features VALUES (
                        ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                    )
                ''', (
                    entry_id, phys.length, phys.molecular_weight, phys.hydrophobicity_mean,
                    phys.hydrophobicity_std, phys.hydrophobicity_max, phys.hydrophobicity_min,
                    phys.net_charge, phys.positive_residues, phys.negative_residues,
                    phys.charge_density, phys.beta_propensity_mean, phys.beta_propensity_max,
                    phys.aggregation_propensity_mean, phys.aggregation_propensity_max,
                    phys.aromatic_fraction, phys.aliphatic_fraction, phys.polar_fraction,
                    phys.charged_fraction, phys.has_polyq, phys.has_glycine_rich, phys.has_proline
                ))
            
            if include_features and entry._composition:
                comp = entry._composition
                cursor.execute('''
                    INSERT OR REPLACE INTO sequence_composition VALUES (
                        ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                    )
                ''', (
                    entry_id,
                    comp.aa_freq.get('A', 0), comp.aa_freq.get('C', 0), comp.aa_freq.get('D', 0),
                    comp.aa_freq.get('E', 0), comp.aa_freq.get('F', 0), comp.aa_freq.get('G', 0),
                    comp.aa_freq.get('H', 0), comp.aa_freq.get('I', 0), comp.aa_freq.get('K', 0),
                    comp.aa_freq.get('L', 0), comp.aa_freq.get('M', 0), comp.aa_freq.get('N', 0),
                    comp.aa_freq.get('P', 0), comp.aa_freq.get('Q', 0), comp.aa_freq.get('R', 0),
                    comp.aa_freq.get('S', 0), comp.aa_freq.get('T', 0), comp.aa_freq.get('V', 0),
                    comp.aa_freq.get('W', 0), comp.aa_freq.get('Y', 0),
                    comp.tiny_fraction, comp.small_fraction, comp.large_fraction
                ))
        except sqlite3.IntegrityError:
            pass
    
    # Views
    cursor.execute('CREATE VIEW IF NOT EXISTS amyloid_entries AS SELECT * FROM entries WHERE is_amyloid = 1')
    cursor.execute('CREATE VIEW IF NOT EXISTS experimental_only AS SELECT * FROM entries WHERE evidence_type != "computational"')
    cursor.execute('CREATE VIEW IF NOT EXISTS high_confidence AS SELECT * FROM entries WHERE confidence >= 70')
    cursor.execute('''
        CREATE VIEW IF NOT EXISTS ml_features AS
        SELECT e.*, p.*, c.*
        FROM entries e
        LEFT JOIN physicochemical_features p ON e.id = p.entry_id
        LEFT JOIN sequence_composition c ON e.id = c.entry_id
        WHERE e.sequence IS NOT NULL AND length(e.sequence) >= 4
    ''')
    
    conn.commit()
    
    cursor.execute('SELECT COUNT(*) FROM entries')
    total = cursor.fetchone()[0]
    cursor.execute('SELECT COUNT(*) FROM entries WHERE is_amyloid = 1')
    amyloid = cursor.fetchone()[0]
    
    logger.info(f"SQLite database created: {total} entries ({amyloid} amyloid)")
    conn.close()
