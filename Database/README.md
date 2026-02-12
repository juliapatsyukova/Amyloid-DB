# Data Unification Pipeline and Source Files

This directory contains the core data processing pipeline (`amyloid_pipeline.py`) and the raw source data files used to build the unified amyloid dataset.

## 1. Purpose

The primary goal of this pipeline is to address the significant challenge of data heterogeneity in amyloid research. It systematically integrates data from multiple public databases, each with its own format and standards, into a single, clean, and evidence-weighted dataset. This process involves parsing different file types, standardizing data into a common schema, resolving conflicting experimental labels, and producing ready-to-use files for machine learning and analysis.

## 2. File Descriptions

-   `amyloid_pipeline.py`: The main Python script that orchestrates the entire data unification process. It contains the logic for parsing, standardizing, deduplicating, and resolving conflicts.

-   **Source Data Files**:
    -   `waltzdb.csv`: Data from the WALTZ-DB 2.0 database, focusing on experimentally characterized amyloid-forming hexapeptides.
    -   `crossbetadb.json`: Data from the Cross-Beta DB, containing high-quality, naturally occurring cross-β amyloid structures.
    -   `AmyLoad_session_unlogged_unlogged_.csv`: Data from the AmyLoad database, a repository of peptides and proteins with experimentally verified amyloid propensity.
    -   `aggregatingpeptides.xlsx`: Data from the CPAD 2.0 database, covering aggregating peptides.
    -   `amyloidstructure.xlsx`: Data from the CPAD 2.0 database, covering amyloid structures.

## 3. Pipeline Workflow

The `amyloid_pipeline.py` script executes a multi-phase process to ensure data quality and biological validity.

### Phase 1: Parsing and Standardization

For each source database, a dedicated parser class (e.g., `WaltzDBParser`, `CrossBetaDBParser`) reads the raw file. Each record is then transformed into a standardized `AmyloidSequenceEntry` or `AmyloidStructureEntry` object. This ensures that all data, regardless of its origin, conforms to a unified schema.

### Phase 2: Evidence Weighting

During parsing, each entry is assigned an `evidence_weight` based on the experimental method used. This weight is determined by a predefined dictionary (`EVIDENCE_WEIGHTS`) that prioritizes higher-quality evidence:

-   **Structural data** (Cryo-EM, ssNMR, X-ray): Highest weight (3.0)
-   **Kinetic data**: Medium weight (2.0)
-   **Staining/Binding assays** (ThT, TEM): Lower weight (1.0)
-   **Literature-curated**: Lowest weight (0.5)

### Phase 3: Unification and Deduplication

All standardized sequence entries are collected by the `SequenceDatasetUnifier`. A key feature of this pipeline is its biologically meaningful deduplication strategy. Instead of simply removing duplicate sequences, it uses a compound key to define uniqueness:

`key = (sequence, protein_id, region_start, region_end)`

This ensures that identical short peptide sequences found in different proteins or at different locations are treated as distinct biological entities, preventing incorrect data merging.

### Phase 4: Conflict Detection and Resolution

The pipeline identifies cases where the same unique entry (based on the key above) is reported with conflicting labels (e.g., as both 'amyloid' and 'non-amyloid' in different databases). These conflicts are resolved using a weighted consensus algorithm:

1.  The evidence weights for each label ('amyloid' vs. 'non-amyloid') are summed up.
2.  If one label's total weight is **more than 50% greater** than the other's, it is chosen as the consensus label.
3.  If the evidence is not decisive, the entry is removed from the consensus set and stored in a separate `conflicts_sequences.csv` file for manual review.

### Phase 5: Output Generation

The pipeline generates four primary output files:

-   `consensus_sequences.csv`: The main dataset containing clean, deduplicated entries with a single, high-confidence label.
-   `non_amyloid_sequences.csv`: A subset of the consensus dataset containing only the 'non-amyloid' entries.
-   `conflicts_sequences.csv`: A dataset containing all entries that had conflicting labels and could not be resolved automatically.
-   `unified_structures.csv`: A separate dataset containing all entries related to protein structures (from PDB).

## 4. How to Run the Pipeline

To execute the pipeline, you need to configure and run the `amyloid_pipeline.py` script.

1.  **Set File Paths**: Open `amyloid_pipeline.py` and modify the file paths in the `run_pipeline` function call at the end of the script to point to the correct locations of your input data files.

    ```python
    # Example configuration at the end of the script
    if __name__ == "__main__":
        run_pipeline(
            waltzdb_path="/path/to/your/waltzdb.csv",
            crossbeta_path="/path/to/your/crossbetadb.json",
            amyload_path="/path/to/your/AmyLoad_session_unlogged_unlogged_.csv",
            cpad_peptide_path="/path/to/your/aggregatingpeptides.xlsx",
            cpad_structure_path="/path/to/your/amyloidstructure.xlsx"
        )
    ```

2.  **Execute the Script**:

    ```bash
    python amyloid_pipeline.py
    ```

    The script will log its progress to the console and to a `pipeline.log` file, and the final CSV files will be saved in the specified output directory.

## 5. Output Schema

The main output file, `consensus_sequences.csv`, will have the following columns, corresponding to the `AmyloidSequenceEntry` data structure:

| Column                 | Description                                                                 |
| ---------------------- | --------------------------------------------------------------------------- |
| `sequence`             | The amino acid sequence of the peptide or protein fragment.                 |
| `protein_id`           | The UniProt ID of the parent protein, if available.                         |
| `species`              | The source organism of the protein.                                         |
| `region_start`         | The starting position of the sequence within the parent protein.            |
| `region_end`           | The ending position of the sequence within the parent protein.              |
| `experimental_label`   | The final consensus label: 'amyloid' or 'non-amyloid'.                      |
| `evidence_type`        | The category of the primary experimental evidence (e.g., 'structural').     |
| `experimental_methods` | A comma-separated list of the specific experimental methods used.           |
| `evidence_weight`      | The numeric weight of the evidence supporting the final label.              |
| `source_database`      | The original database from which the highest-evidence record was sourced.   |
| `pmid_or_doi`          | The publication identifier (PMID or DOI) for the supporting evidence.       |
