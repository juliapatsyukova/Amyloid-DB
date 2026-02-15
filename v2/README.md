# Unified Amyloid Database Pipeline

A modular Python pipeline for integrating, standardizing, and analyzing data from multiple amyloid protein databases.

## 1. Purpose

This pipeline addresses data heterogeneity in amyloid research by integrating **10 public databases** into a single, evidence-weighted dataset:

- Parsing different file types (CSV, TSV, JSON, Excel)
- Standardizing data into a common schema
- Resolving conflicting experimental labels using evidence weighting
- Computing physicochemical and sequence composition features
- Producing ML-ready outputs (TSV, FASTA, SQLite)

## 2. Modular Architecture

```
amyloid_pipeline/
├── __init__.py              # Package exports
├── __main__.py              # CLI entry point
├── config.py                # Constants, evidence weights, method rules
├── models.py                # AmyloidEntry, PhysicochemicalFeatures, SequenceComposition
├── filters.py               # AmyloidFilter chainable filtering system
├── unifier.py               # DatasetUnifier for deduplication/consensus
├── run_pipeline.py          # Main pipeline orchestration
│
├── parsers/                 # Database-specific parsers
│   ├── __init__.py          # PARSERS registry
│   ├── base.py              # BaseParser with inference methods
│   ├── waltzdb.py           # WALTZ-DB 2.0 (hexapeptides)
│   ├── crossbeta.py         # Cross-Beta DB (cross-β structures)
│   ├── amyload.py           # AmyLoad (experimental propensity)
│   ├── amyloid_explorer.py  # AmyloidExplorer (structural)
│   ├── amyloid_atlas.py     # AmyloidAtlas (cryo-EM/NMR)
│   ├── amylobase.py         # Amylobase (kinetic data)
│   ├── amypro.py            # AmyPro (functional amyloids)
│   ├── amylograph.py        # AmyloGraph (cross-seeding)
│   └── cpad.py              # CPAD 2.0 peptides & structures
│
├── features/                # Computed features for ML
│   ├── __init__.py
│   ├── physicochemical.py   # Hydrophobicity, charge, β-propensity
│   └── composition.py       # AA frequencies, dipeptides
│
├── export/                  # Output formats
│   ├── __init__.py
│   ├── tsv.py               # TSV/CSV export
│   ├── fasta.py             # FASTA sequences
│   └── sqlite.py            # SQLite database
│
└── utils/                   # Shared utilities
    ├── __init__.py
    ├── sequence.py          # Validation, cleaning, method mapping
    └── fetchers.py          # UniProt API, sequence cache
```

## 3. Supported Databases

| Database | Format | Content | Parser |
|----------|--------|---------|--------|
| WALTZ-DB 2.0 | CSV/TSV | Experimentally characterized hexapeptides | `WaltzDBParser` |
| Cross-Beta DB | JSON | High-quality cross-β amyloid structures | `CrossBetaDBParser` |
| AmyLoad | CSV | Peptides with verified amyloid propensity | `AmyLoadParser` |
| AmyloidExplorer | TSV | Structural data with disease context | `AmyloidExplorerParser` |
| AmyloidAtlas | TSV | Cryo-EM/NMR fibril structures | `AmyloidAtlasParser` |
| Amylobase | TSV | Aggregation kinetics data | `AmylobaseParser` |
| AmyPro | TSV | Functional and pathogenic amyloids | `AmyProParser` |
| AmyloGraph | CSV | Cross-seeding interactions | `AmyloGraphParser` |
| CPAD 2.0 (peptides) | Excel | Aggregating peptides | `CPADPeptideParser` |
| CPAD 2.0 (structures) | Excel | Amyloid structures | `CPADStructureParser` |

## 4. Pipeline Workflow

### Phase 1: Parsing and Standardization

Each database has a dedicated parser that transforms records into standardized `AmyloidEntry` objects:

```python
@dataclass
class AmyloidEntry:
    record_id: str
    source_db: str
    protein_name: str
    uniprot_id: str
    organism: str
    protein_family: str          # Inferred or fetched from UniProt
    region_start: Optional[int]
    region_end: Optional[int]
    sequence: str
    is_amyloid: bool
    experimental_label: str
    structure_type: str          # fibril, oligomer, crystal, etc.
    aggregate_type: str          # in_vitro, patient_derived, etc.
    pathogenicity: str           # pathogenic, functional, unknown
    evidence_type: str
    evidence_weight: float
    confidence: int
    # ... additional fields
```

### Phase 2: Methods Mapping and Evidence Weighting

Experimental methods are mapped to universal vocabulary using regex rules in `config.py`. Each method is assigned an evidence tier and weight:

| Tier | Evidence Type | Weight | Confidence | Examples |
|------|---------------|--------|------------|----------|
| 4 | Structural | 3.0 | 90 | XRD (cross-β diffraction) |
| 3 | Structural | 3.0 | 70 | Cryo-EM, ssNMR, TEM, AFM |
| 5 | Spectroscopy | 3.0 | 80 | FTIR, CD, Raman |
| 2 | Kinetic | 2.0 | 50 | Aggregation kinetics, seeding |
| 1 | Staining/Binding | 1.0 | 35 | ThT, Congo Red, Proteostat |
| 0 | Literature/Computational | 0.5/0.0 | 10 | Curated, PASTA, TANGO |

### Phase 3: Unification and Deduplication

The `DatasetUnifier` performs biologically meaningful deduplication using a compound key:

```python
key = (sequence, protein_id, region_start, region_end)
```

This ensures identical short peptides in different proteins or positions are treated as distinct biological entities.

### Phase 4: Conflict Detection and Resolution

When the same entry appears with conflicting labels:

1. Evidence weights for each label are summed
2. If one label's weight exceeds the other by **>50%**, it becomes consensus
3. Otherwise, the entry is flagged in `conflicts.json` for manual review

### Phase 5: Feature Computation

**Physicochemical** (22 features):
- Hydrophobicity (mean, std, max, min) — Kyte-Doolittle scale
- Charge (net, density, positive/negative residues)
- β-sheet propensity (mean, max) — Chou-Fasman scale
- Aggregation propensity (mean, max)
- Composition fractions (aromatic, aliphatic, polar, charged)
- Sequence patterns (polyQ, glycine-rich, proline)

**Composition** (23+ features):
- 20 amino acid frequencies
- Top 25 dipeptide frequencies
- Size-based fractions (tiny, small, large residues)

### Phase 6: Export

| File | Description |
|------|-------------|
| `consensus_unified.tsv` | All deduplicated entries with features |
| `amyloid_positive.tsv` | Amyloid-forming entries only |
| `non_amyloid.tsv` | Non-amyloid entries only |
| `amyloid_sequences.fasta` | Deduplicated FASTA |
| `amyloid_database.sqlite` | SQLite with indexed tables |
| `conflicts.json` | Unresolved conflicts |

## 5. Installation

```bash
pip install pandas openpyxl beautifulsoup4
```

## 6. Usage

### Command Line

```bash
# Basic
python -m amyloid_pipeline \
    --waltzdb waltzdb.csv \
    --crossbeta crossbetadb.json \
    -o output

# Full pipeline with all databases
python -m amyloid_pipeline \
    --waltzdb waltzdb.csv \
    --crossbeta crossbetadb.json \
    --amyload AmyLoad.csv \
    --amyloid-explorer AmyloidExplorer.tsv \
    --amyloid-atlas AmyloidAtlas.tsv \
    --amylobase amylobase.txt \
    --amypro amypro.txt \
    --amylograph amylograph.csv \
    --cpad-peptides aggregatingpeptides.xlsx \
    --cpad-structures amyloidstructure.xlsx \
    -o output \
    --fetch-sequences \
    --exclude-computational \
    -v
```

### Python API

```python
from amyloid_pipeline import run_pipeline, AmyloidFilter

results = run_pipeline(
    input_files={
        'waltzdb': 'waltzdb.csv',
        'crossbeta': 'crossbetadb.json',
    },
    output_dir='output',
    fetch_sequences=True,
    compute_features=True,
    export_sqlite=True
)

# Filter
f = AmyloidFilter()
f.exclude_evidence_type('computational')
f.organism('Homo sapiens')
f.min_confidence(50)
filtered = f.apply(results['consensus'])
```

## 7. Filtering System

```python
from amyloid_pipeline import AmyloidFilter

f = AmyloidFilter()
f.evidence_type(['structural', 'kinetic'])   # Include only these
f.exclude_evidence_type('computational')      # Exclude predictions
f.organism('Homo sapiens')                    # Human proteins
f.pathogenicity(['pathogenic'])               # Disease-associated
f.structure_type(['fibril'])                  # Fibrillar only
f.has_sequence(min_length=6)                  # Min 6 residues
f.has_pdb()                                   # Has structure

filtered = f.apply(entries)
print(f"Applied: {f}")
```

### Available Filters

| Method | Description |
|--------|-------------|
| `.evidence_type(types)` | Include only specified evidence types |
| `.exclude_evidence_type(types)` | Exclude specified types |
| `.organism(name)` | Filter by organism (contains match) |
| `.min_confidence(n)` | Minimum confidence score |
| `.is_amyloid(bool)` | Amyloid status |
| `.structure_type(types)` | fibril, oligomer, crystal, aggregate |
| `.aggregate_type(types)` | in_vitro, patient_derived, synthetic |
| `.pathogenicity(types)` | pathogenic, functional, non_pathogenic |
| `.source_db(dbs)` | Specific databases only |
| `.protein_family(names)` | Family keyword match |
| `.disease(names)` | Disease keyword match |
| `.has_sequence(min_length)` | Minimum sequence length |
| `.has_pdb()` | Has PDB structure |
| `.has_uniprot()` | Has UniProt ID |
| `.custom(func, desc)` | Custom filter function |

## 8. SQLite Database

### Schema

```sql
-- Main entries table (indexed on key columns)
entries (
    id, record_id, source_db, protein_name, uniprot_id, organism,
    protein_family, region_start, region_end, sequence, sequence_length,
    is_amyloid, experimental_label, category, structure_type,
    aggregate_type, pathogenicity, pdb_id, evidence_type, evidence_weight,
    confidence, disease, mutation, doi, pmid, ...
)

-- Physicochemical features (1:1 with entries)
physicochemical_features (
    entry_id, length, molecular_weight, hydrophobicity_mean, net_charge,
    beta_propensity_mean, aggregation_propensity_mean, ...
)

-- Sequence composition (1:1 with entries)
sequence_composition (
    entry_id, aa_A, aa_C, ..., aa_Y, tiny_fraction, small_fraction, ...
)
```

### Pre-built Views

```sql
SELECT * FROM amyloid_entries;      -- is_amyloid = 1
SELECT * FROM experimental_only;    -- evidence_type != 'computational'
SELECT * FROM high_confidence;      -- confidence >= 70
SELECT * FROM ml_features;          -- Joined with all features
```

### Example Queries

```python
import sqlite3
import pandas as pd

conn = sqlite3.connect('output/amyloid_database.sqlite')

# Human pathogenic fibrils with features
df = pd.read_sql_query('''
    SELECT e.protein_name, e.sequence, p.hydrophobicity_mean, p.beta_propensity_mean
    FROM entries e
    JOIN physicochemical_features p ON e.id = p.entry_id
    WHERE e.organism LIKE '%Homo sapiens%'
    AND e.pathogenicity = 'pathogenic'
    AND e.structure_type = 'fibril'
''', conn)
```

## 9. Adding a New Parser

1. Create `parsers/newdb.py`:

```python
from .base import BaseParser
from ..models import AmyloidEntry
from ..utils.sequence import clean_value, validate_sequence, map_method_to_universal

class NewDBParser(BaseParser):
    def parse(self, filepath: str, **kwargs):
        self.logger.info(f"Parsing {self.db_name} from {filepath}...")
        
        with open(filepath) as f:
            for idx, record in enumerate(parse_records(f)):
                seq = clean_value(record.get('sequence', '')).upper()
                if not validate_sequence(seq):
                    continue
                
                raw_method = record.get('method', '')
                universal, evidence_type, weight, confidence = map_method_to_universal(raw_method)
                
                entry = AmyloidEntry(
                    record_id=f"NEW_{idx}",
                    source_db="NewDB",
                    sequence=seq,
                    is_amyloid=True,
                    evidence_type=evidence_type,
                    evidence_weight=weight,
                    confidence=confidence,
                )
                
                if entry.is_valid():
                    self.entries.append(entry)
        
        self.log_result()
        return self.entries
```

2. Register in `parsers/__init__.py`:

```python
from .newdb import NewDBParser
PARSERS['newdb'] = NewDBParser
```

3. Add CLI argument in `__main__.py`:

```python
inputs.add_argument('--newdb', help='NewDB file')
```

## 10. ML Integration Example

```python
from amyloid_pipeline import run_pipeline, AmyloidFilter
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier

# 1. Run pipeline
results = run_pipeline(
    input_files={'waltzdb': 'waltzdb.csv', 'crossbeta': 'data.json'},
    output_dir='output',
    compute_features=True
)

# 2. Filter (experimental evidence only)
f = AmyloidFilter()
f.exclude_evidence_type('computational')
f.has_sequence(min_length=6)
ml_entries = f.apply(results['consensus'])

# 3. Convert to DataFrame
data = [e.to_dict(include_features=True) for e in ml_entries]
df = pd.DataFrame(data)

# 4. Prepare features
feature_cols = [c for c in df.columns if c.startswith('phys_') or c.startswith('comp_')]
X = df[feature_cols].fillna(0)
y = df['is_amyloid'].astype(int)

# 5. Train model
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y)
clf = RandomForestClassifier(n_estimators=100)
clf.fit(X_train, y_train)
print(f"Accuracy: {clf.score(X_test, y_test):.3f}")
```

## 11. Output Schema Reference

### AmyloidEntry Fields

| Field | Type | Description |
|-------|------|-------------|
| `record_id` | str | Unique identifier within source |
| `source_db` | str | Source database name |
| `protein_name` | str | Protein name |
| `uniprot_id` | str | UniProt accession |
| `organism` | str | Source organism |
| `protein_family` | str | Protein family (inferred/fetched) |
| `region_start` | int | Start position in protein |
| `region_end` | int | End position in protein |
| `sequence` | str | Amino acid sequence |
| `is_amyloid` | bool | Amyloid-forming status |
| `experimental_label` | str | 'amyloid' or 'non-amyloid' |
| `structure_type` | str | fibril, oligomer, crystal, aggregate, monomer |
| `aggregate_type` | str | in_vitro, ex_vivo, patient_derived, synthetic |
| `pathogenicity` | str | pathogenic, functional, non_pathogenic |
| `evidence_type` | str | structural, kinetic, staining_binding, etc. |
| `evidence_weight` | float | Numeric weight (0.0-3.0) |
| `confidence` | int | Confidence score (0-100) |
| `pdb_id` | str | PDB structure ID |
| `disease` | str | Associated disease |
| `mutation` | str | Mutation information |
| `doi` | str | Publication DOI |
| `pmid` | str | PubMed ID |

## 12. Authors

- Xenia Sukhanova[^1]
- Julia Patsiukova[^2]

[^1]: Laboratory of Amyloid Studies, Saint Petersburg State University
[^2]: ITMO University
