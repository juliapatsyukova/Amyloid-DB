# Methods Mapping: From Raw to Universal

This directory contains the script `build_universal_methods_mapping.py`, which is a critical first step in the data unification pipeline. Its purpose is to scan all the raw source database files, extract the various names used for experimental methods, and map them to a standardized, universal vocabulary. This process is essential for the evidence-weighting system used in the main `amyloid_pipeline.py`.

## 1. The Problem: Inconsistent Method Names

Different databases report experimental methods in different ways:

-   **WALTZ-DB** uses separate columns for each assay (e.g., `TEM Staining`, `Th-T Binding`). The method is implied by the presence of a value in that column.
-   **Cross-Beta DB** uses a text field (`Method used`) with free-form descriptions like "Cryo-EM" or "X-ray diffraction".
-   **CPAD** lists computational predictors in separate columns and also has a "Method" column for its structural data.

To treat this data systematically, we must first standardize these varied and inconsistent names into a single, clean vocabulary.

## 2. The Solution: `build_universal_methods_mapping.py`

This script automates the creation of a master mapping table. It works in three main stages:

### Stage 1: Extraction

The script contains dedicated extractor functions (e.g., `extract_waltz_methods`, `extract_crossbeta_methods`) for each source database. Each function knows how to find and extract the raw method names from its specific file format. It scans all files and compiles a comprehensive list of every unique method string found.

### Stage 2: Rule-Based Mapping

A central component of the script is the `RULES` list. This is a list of dictionaries, where each dictionary defines a regular expression `pattern` and the corresponding `new` universal name. The script iterates through every raw method string and applies these regex rules to find a match.

**Example Rule**:
```python
{
    "pattern": r"\bcryo[-\s]?em\b", 
    "new": "Cryo-EM", 
    "tier": 3, 
    "category": "microscopy_structural",
    "notes": "microscopy; fibrils/structure; binary-ish evidence"
}
```
This rule will match strings like "cryo-em", "cryo em", or "Cryo-EM" and map them all to the universal name "Cryo-EM".

### Stage 3: Tiering and Confidence Scoring

Once a method is mapped to a universal name, it is assigned a **tier** and a **confidence level**. This is the basis for the evidence-weighting system in the main pipeline.

The tier system is defined in the `TIER_CONF` dictionary and categorizes methods based on the strength of the evidence they provide:

| Tier | Category                  | Confidence | Description                                       |
| :--- | :------------------------ | :--------- | :------------------------------------------------ |
| 5    | Spectroscopy              | 80         | Secondary structure signatures (FTIR, CD)         |
| 4    | XRD (Cross-β)             | 90         | Definitive cross-β diffraction pattern            |
| 3    | Structural Microscopy/NMR | 70         | Direct visualization of fibrils (Cryo-EM, TEM, AFM) |
| 2    | Kinetics / Biological     | 50         | Mid-strength assays (aggregation kinetics)        |
| 1    | Staining / Binding        | 35         | Dye binding, suggestive evidence (ThT, Congo Red) |
| 0    | Computational/Unspecified | 10         | Predictions or literature-curated entries         |

## 3. Output Files

The script generates two CSV files:

1.  `methods_mapping_table_raw_to_universal_FULL.csv`: A complete log of every raw method string found and its corresponding universal mapping.
2.  `methods_mapping_table_raw_to_universal_CORE.csv`: A deduplicated version of the full table, showing each unique mapping once.

**Output Schema**:

| Column                 | Description                                                              |
| ---------------------- | ------------------------------------------------------------------------ |
| `database_name`        | The source database where the raw name was found.                        |
| `raw_name_of_method`   | The original, unmodified method string from the source file.             |
| `new_name_universal`   | The standardized, clean name for the method.                             |
| `confidence_level`     | The numeric confidence score (0-100) associated with the method's tier.  |
| `tier`                 | The evidence tier (0-5) assigned to the method.                          |
| `category`             | A broader category for the method (e.g., `microscopy_structural`).       |
| `source_column`        | The name of the column in the source file where the method was found.    |
| `notes`                | Explanatory notes about the mapping rule.                                |

## 4. How to Run the Script

This script should be run **before** the main `amyloid_pipeline.py`.

1.  **Set File Paths**: Open `build_universal_methods_mapping.py` and update the `INPUT_FILES` dictionary to point to the correct locations of your source data files.

    ```python
    INPUT_FILES = {
        "WALTZ-DB 2.0": "/path/to/your/waltzdb.csv",
        "Cross-Beta DB": "/path/to/your/crossbetadb.json",
        # ... and so on for all files
    }
    ```

2.  **Execute the Script**:

    ```bash
    python build_universal_methods_mapping.py
    ```

    The script will print its progress and save the output CSV files to the specified `OUTPUT_DIR`.
