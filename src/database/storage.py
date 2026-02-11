"""
Database Module
SQLite storage for extracted protein aggregation data
"""

import sqlite3
import json
from typing import List, Dict, Optional, Any
from datetime import datetime
from contextlib import contextmanager
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Database schema version for migrations
SCHEMA_VERSION = 1


class AmyloidDatabase:
    """SQLite database for amyloid/aggregation data"""
    
    def __init__(self, db_path: str = "data/amyloid_aggregation.db"):
        self.db_path = db_path
        self._init_database()
    
    @contextmanager
    def _connection(self):
        """Context manager for database connections"""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        try:
            yield conn
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            conn.close()
    
    def _init_database(self):
        """Initialize database schema"""
        with self._connection() as conn:
            cursor = conn.cursor()
            
            # Proteins table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS proteins (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    protein_name TEXT NOT NULL,
                    gene_name TEXT,
                    uniprot_id TEXT,
                    organism TEXT,
                    
                    -- Aggregation properties
                    forms_aggregates BOOLEAN DEFAULT 1,
                    aggregate_structure TEXT,
                    reversibility TEXT,
                    assembly_mechanism TEXT,
                    
                    -- Classification
                    functional_class TEXT,
                    biological_role TEXT,
                    associated_disease TEXT,
                    
                    -- Evidence
                    evidence_strength TEXT,
                    extraction_confidence REAL,
                    
                    -- Metadata
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    
                    UNIQUE(protein_name, organism)
                )
            ''')
            
            # Papers table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS papers (
                    pmid TEXT PRIMARY KEY,
                    title TEXT,
                    abstract TEXT,
                    authors TEXT,
                    journal TEXT,
                    year INTEGER,
                    doi TEXT,
                    relevance TEXT,
                    processed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            ''')
            
            # Protein-Paper association (many-to-many)
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS protein_papers (
                    protein_id INTEGER,
                    pmid TEXT,
                    PRIMARY KEY (protein_id, pmid),
                    FOREIGN KEY (protein_id) REFERENCES proteins(id),
                    FOREIGN KEY (pmid) REFERENCES papers(pmid)
                )
            ''')
            
            # Experimental methods
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS experimental_methods (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    protein_id INTEGER,
                    pmid TEXT,
                    method_name TEXT,
                    method_category TEXT,
                    FOREIGN KEY (protein_id) REFERENCES proteins(id),
                    FOREIGN KEY (pmid) REFERENCES papers(pmid)
                )
            ''')
            
            # Cofactors table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS cofactors (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    protein_id INTEGER,
                    cofactor_name TEXT,
                    cofactor_type TEXT,  -- 'assembly', 'disassembly'
                    source TEXT,  -- 'literature', 'string_db'
                    confidence REAL,
                    pmid TEXT,
                    FOREIGN KEY (protein_id) REFERENCES proteins(id)
                )
            ''')
            
            # STRING interactions table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS string_interactions (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    protein_id INTEGER,
                    partner_name TEXT,
                    partner_string_id TEXT,
                    combined_score REAL,
                    experimental_score REAL,
                    database_score REAL,
                    textmining_score REAL,
                    retrieved_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (protein_id) REFERENCES proteins(id)
                )
            ''')
            
            # Create indexes
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_proteins_name ON proteins(protein_name)')
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_proteins_uniprot ON proteins(uniprot_id)')
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_proteins_gene ON proteins(gene_name)')
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_papers_year ON papers(year)')
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_cofactors_type ON cofactors(cofactor_type)')
            
            logger.info(f"Database initialized: {self.db_path}")
    
    # ==========================================================================
    # Protein operations
    # ==========================================================================
    
    def add_protein(self, protein_data: Dict) -> int:
        """Add or update a protein record"""
        with self._connection() as conn:
            cursor = conn.cursor()
            
            # Check if exists
            cursor.execute('''
                SELECT id FROM proteins 
                WHERE protein_name = ? AND (organism = ? OR organism IS NULL)
            ''', (protein_data.get('protein_name'), protein_data.get('organism')))
            
            existing = cursor.fetchone()
            
            if existing:
                # Update existing
                protein_id = existing['id']
                cursor.execute('''
                    UPDATE proteins SET
                        gene_name = COALESCE(?, gene_name),
                        uniprot_id = COALESCE(?, uniprot_id),
                        aggregate_structure = COALESCE(?, aggregate_structure),
                        reversibility = COALESCE(?, reversibility),
                        assembly_mechanism = COALESCE(?, assembly_mechanism),
                        functional_class = COALESCE(?, functional_class),
                        biological_role = COALESCE(?, biological_role),
                        associated_disease = COALESCE(?, associated_disease),
                        evidence_strength = COALESCE(?, evidence_strength),
                        extraction_confidence = MAX(extraction_confidence, ?),
                        updated_at = CURRENT_TIMESTAMP
                    WHERE id = ?
                ''', (
                    protein_data.get('gene_name'),
                    protein_data.get('uniprot_id'),
                    protein_data.get('aggregate_structure'),
                    protein_data.get('reversibility'),
                    protein_data.get('assembly_mechanism'),
                    protein_data.get('functional_class'),
                    protein_data.get('biological_role'),
                    protein_data.get('associated_disease'),
                    protein_data.get('evidence_strength'),
                    protein_data.get('extraction_confidence', 0),
                    protein_id
                ))
            else:
                # Insert new
                cursor.execute('''
                    INSERT INTO proteins (
                        protein_name, gene_name, uniprot_id, organism,
                        forms_aggregates, aggregate_structure, reversibility,
                        assembly_mechanism, functional_class, biological_role,
                        associated_disease, evidence_strength, extraction_confidence
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    protein_data.get('protein_name'),
                    protein_data.get('gene_name'),
                    protein_data.get('uniprot_id'),
                    protein_data.get('organism'),
                    protein_data.get('forms_aggregates', True),
                    protein_data.get('aggregate_structure'),
                    protein_data.get('reversibility'),
                    protein_data.get('assembly_mechanism'),
                    protein_data.get('functional_class'),
                    protein_data.get('biological_role'),
                    protein_data.get('associated_disease'),
                    protein_data.get('evidence_strength'),
                    protein_data.get('extraction_confidence', 0),
                ))
                protein_id = cursor.lastrowid
            
            return protein_id
    
    def get_protein(self, protein_id: int = None, name: str = None) -> Optional[Dict]:
        """Get protein by ID or name"""
        with self._connection() as conn:
            cursor = conn.cursor()
            
            if protein_id:
                cursor.execute('SELECT * FROM proteins WHERE id = ?', (protein_id,))
            elif name:
                cursor.execute('SELECT * FROM proteins WHERE protein_name = ?', (name,))
            else:
                return None
            
            row = cursor.fetchone()
            return dict(row) if row else None
    
    def get_all_proteins(
        self,
        functional_class: str = None,
        min_confidence: float = None,
        limit: int = None
    ) -> List[Dict]:
        """Get proteins with optional filters"""
        with self._connection() as conn:
            cursor = conn.cursor()
            
            query = 'SELECT * FROM proteins WHERE 1=1'
            params = []
            
            if functional_class:
                query += ' AND functional_class = ?'
                params.append(functional_class)
            
            if min_confidence:
                query += ' AND extraction_confidence >= ?'
                params.append(min_confidence)
            
            query += ' ORDER BY extraction_confidence DESC'
            
            if limit:
                query += ' LIMIT ?'
                params.append(limit)
            
            cursor.execute(query, params)
            return [dict(row) for row in cursor.fetchall()]
    
    # ==========================================================================
    # Paper operations
    # ==========================================================================
    
    def add_paper(self, paper_data: Dict) -> str:
        """Add a paper record"""
        with self._connection() as conn:
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT OR REPLACE INTO papers (
                    pmid, title, abstract, authors, journal, year, doi, relevance
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                paper_data.get('pmid'),
                paper_data.get('title'),
                paper_data.get('abstract'),
                json.dumps(paper_data.get('authors', [])),
                paper_data.get('journal'),
                paper_data.get('year'),
                paper_data.get('doi'),
                paper_data.get('relevance'),
            ))
            
            return paper_data.get('pmid')
    
    def link_protein_paper(self, protein_id: int, pmid: str):
        """Create protein-paper association"""
        with self._connection() as conn:
            cursor = conn.cursor()
            cursor.execute('''
                INSERT OR IGNORE INTO protein_papers (protein_id, pmid)
                VALUES (?, ?)
            ''', (protein_id, pmid))
    
    def is_paper_processed(self, pmid: str) -> bool:
        """Check if paper was already processed"""
        with self._connection() as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT 1 FROM papers WHERE pmid = ?', (pmid,))
            return cursor.fetchone() is not None
    
    # ==========================================================================
    # Experimental methods
    # ==========================================================================
    
    def add_experimental_methods(
        self, 
        protein_id: int, 
        pmid: str, 
        methods: List[str]
    ):
        """Add experimental methods for a protein"""
        with self._connection() as conn:
            cursor = conn.cursor()
            
            for method in methods:
                # Categorize method
                category = self._categorize_method(method)
                
                cursor.execute('''
                    INSERT OR IGNORE INTO experimental_methods 
                    (protein_id, pmid, method_name, method_category)
                    VALUES (?, ?, ?, ?)
                ''', (protein_id, pmid, method, category))
    
    def _categorize_method(self, method: str) -> str:
        """Categorize experimental method"""
        method_lower = method.lower()
        
        if any(x in method_lower for x in ['tht', 'thioflavin', 'congo', 'cd ', 'ftir']):
            return 'biophysical'
        elif any(x in method_lower for x in ['em', 'microscopy', 'afm', 'cryo']):
            return 'microscopy'
        elif any(x in method_lower for x in ['crystal', 'nmr', 'structure']):
            return 'structural'
        elif any(x in method_lower for x in ['cell', 'vivo', 'culture']):
            return 'cellular'
        elif any(x in method_lower for x in ['kinetic', 'seeding', 'nucleation']):
            return 'kinetic'
        else:
            return 'other'
    
    # ==========================================================================
    # Cofactor operations
    # ==========================================================================
    
    def add_cofactor(
        self,
        protein_id: int,
        cofactor_name: str,
        cofactor_type: str,
        source: str = 'literature',
        confidence: float = 0.5,
        pmid: str = None
    ):
        """Add cofactor relationship"""
        with self._connection() as conn:
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT INTO cofactors 
                (protein_id, cofactor_name, cofactor_type, source, confidence, pmid)
                VALUES (?, ?, ?, ?, ?, ?)
            ''', (protein_id, cofactor_name, cofactor_type, source, confidence, pmid))
    
    def get_cofactors(
        self, 
        protein_id: int = None, 
        cofactor_type: str = None
    ) -> List[Dict]:
        """Get cofactors with optional filters"""
        with self._connection() as conn:
            cursor = conn.cursor()
            
            query = 'SELECT * FROM cofactors WHERE 1=1'
            params = []
            
            if protein_id:
                query += ' AND protein_id = ?'
                params.append(protein_id)
            
            if cofactor_type:
                query += ' AND cofactor_type = ?'
                params.append(cofactor_type)
            
            cursor.execute(query, params)
            return [dict(row) for row in cursor.fetchall()]
    
    # ==========================================================================
    # STRING interactions
    # ==========================================================================
    
    def add_string_interactions(
        self,
        protein_id: int,
        interactions: List[Dict]
    ):
        """Store STRING interactions for a protein"""
        with self._connection() as conn:
            cursor = conn.cursor()
            
            for interaction in interactions:
                cursor.execute('''
                    INSERT OR REPLACE INTO string_interactions
                    (protein_id, partner_name, partner_string_id, 
                     combined_score, experimental_score, database_score, textmining_score)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                ''', (
                    protein_id,
                    interaction.get('partner_name'),
                    interaction.get('partner_string_id'),
                    interaction.get('combined_score'),
                    interaction.get('experimental_score'),
                    interaction.get('database_score'),
                    interaction.get('textmining_score'),
                ))
    
    # ==========================================================================
    # Statistics and analysis
    # ==========================================================================
    
    def get_statistics(self) -> Dict:
        """Get database statistics"""
        with self._connection() as conn:
            cursor = conn.cursor()
            
            stats = {}
            
            cursor.execute('SELECT COUNT(*) FROM proteins')
            stats['total_proteins'] = cursor.fetchone()[0]
            
            cursor.execute('SELECT COUNT(*) FROM papers')
            stats['total_papers'] = cursor.fetchone()[0]
            
            cursor.execute('SELECT functional_class, COUNT(*) FROM proteins GROUP BY functional_class')
            stats['by_class'] = {row[0]: row[1] for row in cursor.fetchall()}
            
            cursor.execute('SELECT COUNT(DISTINCT cofactor_name) FROM cofactors')
            stats['unique_cofactors'] = cursor.fetchone()[0]
            
            cursor.execute('''
                SELECT cofactor_type, COUNT(DISTINCT cofactor_name) 
                FROM cofactors GROUP BY cofactor_type
            ''')
            stats['cofactors_by_type'] = {row[0]: row[1] for row in cursor.fetchall()}
            
            return stats
    
    def export_to_csv(self, output_path: str):
        """Export proteins to CSV"""
        import csv
        
        proteins = self.get_all_proteins()
        
        if not proteins:
            logger.warning("No proteins to export")
            return
        
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=proteins[0].keys())
            writer.writeheader()
            writer.writerows(proteins)
        
        logger.info(f"Exported {len(proteins)} proteins to {output_path}")


if __name__ == "__main__":
    # Test the database
    db = AmyloidDatabase("data/test_amyloid.db")
    
    # Add test protein
    protein_id = db.add_protein({
        'protein_name': 'Alpha-synuclein',
        'gene_name': 'SNCA',
        'uniprot_id': 'P37840',
        'organism': 'Homo sapiens',
        'aggregate_structure': 'beta-sheet',
        'reversibility': 'irreversible',
        'functional_class': 'pathological',
        'associated_disease': "Parkinson's disease",
        'extraction_confidence': 0.95,
    })
    
    # Add test paper
    db.add_paper({
        'pmid': '12345678',
        'title': 'Test paper about alpha-synuclein',
        'year': 2023,
    })
    
    # Link
    db.link_protein_paper(protein_id, '12345678')
    
    # Add methods
    db.add_experimental_methods(protein_id, '12345678', ['ThT fluorescence', 'TEM'])
    
    # Add cofactor
    db.add_cofactor(protein_id, 'HSP70', 'disassembly', 'literature', 0.8, '12345678')
    
    # Print stats
    print("\nDatabase statistics:")
    stats = db.get_statistics()
    for k, v in stats.items():
        print(f"  {k}: {v}")
