#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build a unified methods mapping table directly from SOURCE database files.

Goal output (main columns):
database_name | raw_name_of_method | new_name_universal | confidence_level

We also keep helpful debug columns:
tier | category | source_column | notes

Sources supported:
- WALTZ-DB 2.0 (CSV): method presence stored across assay columns (native)
- Cross-Beta DB (JSON): "Method used" values (native)
- AmyLoad (CSV): usually literature-curated (may have extra columns; we handle both)
- CPAD peptides (Excel): predictor columns (native) + PDB links via "Amyloid Structure" (native)
- CPAD structures (Excel): "Method" column (native)

Author intent:
- Extract "raw method names" from the original sources (their own columns/values),
  then map to a universal vocabulary + confidence (0..100) using your tier logic.
"""

import os
import re
import sys
import json
import pandas as pd
import openpyxl
from typing import List, Dict, Tuple, Optional, Set

# =========================
# CONFIG: set your source paths
# =========================
INPUT_FILES = {
    "WALTZ-DB 2.0": "/Users/julia_patsiukova/Downloads/Amyloid/waltzdb.csv",
    "Cross-Beta DB": "/Users/julia_patsiukova/Downloads/Amyloid/crossbetadb.json",
    "AmyLoad": "/Users/julia_patsiukova/Downloads/Amyloid/AmyLoad_session_unlogged_unlogged_.csv",
    "CPAD 2.0 peptides": "/Users/julia_patsiukova/Downloads/Amyloid/CPAD/aggregating peptides.xlsx",
    "CPAD 2.0 structures": "/Users/julia_patsiukova/Downloads/Amyloid/CPAD/amyloid structure.xlsx",
}

OUTPUT_DIR = "/Users/julia_patsiukova/Downloads/Amyloid/amyloid_datasets_corrected"
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUT_CSV_FULL = os.path.join(OUTPUT_DIR, "methods_mapping_table_raw_to_universal_FULL.csv")
OUT_CSV_CORE = os.path.join(OUTPUT_DIR, "methods_mapping_table_raw_to_universal_CORE.csv")

# CPAD sheet names (as in your files)
CPAD_PEPTIDE_SHEET = "peptide"
CPAD_STRUCT_SHEET = "final data"

# =========================
# Your tier logic (0..5) + confidence (0..100)
# You can tune these numbers; these are priors per method class.
# =========================
# 0: suspicion marker / computational / unspecified
# 1: staining / chemical binding
# 2: biological tags + fluorescence (and similar mid-strength assays)
# 3: physical structural (microscopy + NMR + "structural X-ray" if not XRD cross-beta)
# 4: XRD / fiber diffraction cross-beta pattern (binary)
# 5: spectroscopy (FTIR/CD/Raman) – you listed it as highest tier; keep as such if you mean it.
TIER_CONF = {
    0: 10,
    1: 35,
    2: 50,
    3: 70,
    4: 90,
    5: 80,  # adjust if you truly want spectroscopy > XRD
}

# =========================
# Universal mapping rules (editable)
# =========================
# This is intentionally lightweight: it creates a *starting* universal vocabulary.
# You will expand it as you discover more raw strings from Cross-Beta, etc.

RULES: List[Dict] = [
    # --- XRD / diffraction cross-beta (tier 4) ---
    {"pattern": r"\bxrd\b|\bx[-\s]?ray diffraction\b|\bx-ray diffraction\b|\bfiber diffraction\b",
     "new": "XRD (cross-β diffraction)", "tier": 4, "category": "xrd_cross_beta",
     "notes": "binary: cross-β diffraction pattern yes/no"},

    # --- microscopy (tier 3) ---
    {"pattern": r"\bcryo[-\s]?em\b", "new": "Cryo-EM", "tier": 3, "category": "microscopy_structural",
     "notes": "microscopy; fibrils/structure; binary-ish evidence"},
    {"pattern": r"\btem\b|\btransmission electron microscopy\b|\belectron microscopy\b",
     "new": "Electron microscopy (TEM/EM)", "tier": 3, "category": "microscopy_structural",
     "notes": "microscopy; fibrils visible; binary fibrils yes/no"},
    {"pattern": r"\bafm\b|\batomic force microscopy\b",
     "new": "AFM", "tier": 3, "category": "microscopy_structural",
     "notes": "microscopy; fibrils morphology"},

    # --- NMR (split inside) (tier 3 default; you can adjust solution vs solid-state) ---
    {"pattern": r"\bsolid[-\s]?state nmr\b|\bssnmr\b",
     "new": "NMR (solid-state)", "tier": 3, "category": "structural",
     "notes": "structural evidence for fibrils"},
    {"pattern": r"\bsolution nmr\b",
     "new": "NMR (solution)", "tier": 2, "category": "structural",
     "notes": "often not direct fibril evidence; tune if needed"},
    {"pattern": r"\bnmr\b",
     "new": "NMR (unspecified)", "tier": 3, "category": "structural",
     "notes": "structural; consider splitting later"},

    # --- structural X-ray / crystallography (not diffraction cross-beta unless explicitly XRD) ---
    {"pattern": r"\bcrystallograph\b|\bx-ray\b|\bxray\b",
     "new": "X-ray (structural)", "tier": 3, "category": "structural",
     "notes": "structural X-ray; NOT assumed cross-β diffraction unless XRD matched"},
    {"pattern": r"\bmicroed\b|\belectron crystallography\b",
     "new": "Electron crystallography / microED", "tier": 3, "category": "structural",
     "notes": "structural method"},

    # --- spectroscopy (tier 5) ---
    {"pattern": r"\bftir\b|\binfrared\b",
     "new": "FTIR", "tier": 5, "category": "spectroscopy",
     "notes": "spectroscopy; secondary structure signatures"},
    {"pattern": r"\bcircular dichroism\b|\bcd\b",
     "new": "CD", "tier": 5, "category": "spectroscopy",
     "notes": "spectroscopy; secondary structure"},
    {"pattern": r"\braman\b",
     "new": "Raman", "tier": 5, "category": "spectroscopy",
     "notes": "spectroscopy"},

    # --- staining / chemical binding (tier 1) ---
    {"pattern": r"\bthioflavin[-\s]?t\b|\btht\b|\bth[-\s]?t\b",
     "new": "ThT binding", "tier": 1, "category": "staining_chemical",
     "notes": "dye binding; propensity evidence"},
    {"pattern": r"\bcongo[-\s]?red\b",
     "new": "Congo Red binding", "tier": 1, "category": "staining_chemical",
     "notes": "dye binding; propensity evidence"},
    {"pattern": r"\bproteostat\b",
     "new": "Proteostat binding", "tier": 1, "category": "staining_chemical",
     "notes": "dye binding; propensity evidence"},

    # --- kinetics-ish (you placed kinetics earlier in your pipeline; here map to tier 2) ---
    {"pattern": r"\bkinetic\b|\blag\b|\baggregation rate\b|\bseeding\b|\baggregation assay\b",
     "new": "Aggregation kinetics", "tier": 2, "category": "kinetics",
     "notes": "propensity / intermediate evidence"},

    # --- computational predictors (tier 0) ---
    {"pattern": r"\bpasta\b",
     "new": "PASTA 2.0 (prediction)", "tier": 0, "category": "computational_prediction",
     "notes": "computational propensity predictor"},
    {"pattern": r"\baggrescan\b",
     "new": "AGGRESCAN (prediction)", "tier": 0, "category": "computational_prediction",
     "notes": "computational aggregation predictor"},
    {"pattern": r"\btango\b",
     "new": "TANGO (prediction)", "tier": 0, "category": "computational_prediction",
     "notes": "computational aggregation predictor"},
    {"pattern": r"\bnuaprpred\b",
     "new": "NuAPRpred (prediction)", "tier": 0, "category": "computational_prediction",
     "notes": "computational amyloid propensity predictor"},

    # --- literature curated / unspecified (tier 0) ---
    {"pattern": r"\bliterature[-\s]?curated\b|\bcurated\b",
     "new": "Literature-curated", "tier": 0, "category": "literature",
     "notes": "no raw experiment in source file"},
    {"pattern": r"^\s*experimental\s*$|^\s*assay\s*$|^\s*unknown\s*$",
     "new": "Unspecified experimental", "tier": 0, "category": "unspecified",
     "notes": "suspicion marker; needs follow-up"},
]

_COMPILED = [(re.compile(r["pattern"], flags=re.IGNORECASE), r) for r in RULES]


def norm(s: str) -> str:
    return re.sub(r"\s+", " ", str(s).strip())


def apply_universal(raw_method: str) -> Tuple[str, int, str, int, str]:
    """
    Returns (new_name, tier, category, confidence, notes)
    """
    raw = norm(raw_method)
    if not raw or raw.lower() in {"nan", "na", "n.a.", "none"}:
        return ("(missing)", 0, "missing", TIER_CONF[0], "missing raw method")

    for rx, rule in _COMPILED:
        if rx.search(raw):
            tier = int(rule["tier"])
            conf = int(TIER_CONF[tier])
            return (rule["new"], tier, rule["category"], conf, rule.get("notes", ""))

    # fallback: keep raw but mark unmapped
    return (raw, 0, "unmapped", TIER_CONF[0], "no rule matched; keep raw")


def add_record(records: List[Dict], db: str, raw: str, source_column: str):
    raw = norm(raw)
    if not raw:
        return
    new_name, tier, category, conf, notes = apply_universal(raw)
    records.append({
        "database_name": db,
        "raw_name_of_method": raw,
        "new_name_universal": new_name,
        "confidence_level": conf,
        "tier": tier,
        "category": category,
        "source_column": source_column,
        "notes": notes,
    })


# =========================
# EXTRACTORS
# =========================

def extract_waltz_methods(path: str) -> List[Dict]:
    """
    WALTZ-DB: methods are represented by assay-specific columns.
    We treat the *column names* as native method indicators if they occur.
    """
    df = pd.read_csv(path)
    records: List[Dict] = []

    # native columns -> raw method token
    col_map = {
        "TEM Staining": "TEM",
        "Th-T Binding": "ThT",
        "FTIR peaks": "FTIR",
        "Proteostat binding": "Proteostat",
    }

    for col, raw_method in col_map.items():
        if col not in df.columns:
            continue
        # presence: any row not n.a./nan/empty
        s = df[col].astype(str).str.strip().str.lower()
        used = ~s.isin({"nan", "n.a.", "na", "", "none"})
        if used.any():
            add_record(records, "WALTZ-DB 2.0", raw_method, source_column=col)

    return records


def extract_crossbeta_methods(path: str) -> List[Dict]:
    """
    Cross-Beta DB: raw methods are stored as unique text values in 'Method used'.
    """
    records: List[Dict] = []
    with open(path, "r") as f:
        data = json.load(f)

    field = "Method used"
    values: Set[str] = set()
    for rec in data:
        v = norm(rec.get(field, "") or "")
        if not v or v.lower() in {"nan", "n/a", "na", "none"}:
            continue
        values.add(v)

    for v in sorted(values):
        add_record(records, "Cross-Beta DB", v, source_column=field)

    return records


def extract_amyload_methods(path: str) -> List[Dict]:
    """
    AmyLoad: often literature curated.
    If there is an explicit method-like column, we capture it, otherwise we record Literature-curated.
    """
    records: List[Dict] = []
    df = pd.read_csv(path)

    # Heuristic: find a column that looks like method/assay/experimental
    method_cols = [c for c in df.columns if re.search(r"(method|assay|experiment)", c, flags=re.IGNORECASE)]
    if method_cols:
        col = method_cols[0]
        vals = set(norm(v) for v in df[col].dropna().astype(str).tolist())
        vals = {v for v in vals if v and v.lower() not in {"nan", "n/a", "na", "none"}}
        if vals:
            for v in sorted(vals):
                add_record(records, "AmyLoad", v, source_column=col)
            return records

    # fallback: your stated interpretation
    add_record(records, "AmyLoad", "Literature-curated", source_column="(implicit)")
    return records


def parse_pdb_list(cell: str) -> List[str]:
    """
    CPAD peptides: 'Amyloid Structure' contains 'No structures' OR comma-separated PDB IDs.
    """
    s = norm(cell)
    if not s or s.lower() in {"nan", "none"}:
        return []
    if s.lower().startswith("no structures"):
        return []
    parts = [p.strip() for p in s.split(",")]
    pdbs = []
    for p in parts:
        p2 = re.sub(r"[^A-Za-z0-9]", "", p).upper()
        if len(p2) == 4:
            pdbs.append(p2)
    return pdbs


def extract_cpad_methods(peptides_xlsx: str, structures_xlsx: str) -> List[Dict]:
    """
    CPAD:
    - peptides: native predictor columns (PASTA/AGGRESCAN/Tango/NuAPRpred) + PDB links via 'Amyloid Structure'
    - structures: native 'Method' column values
    """
    records: List[Dict] = []

    # --- structures ---
    struct = pd.read_excel(structures_xlsx, sheet_name=CPAD_STRUCT_SHEET)
    if "PDB-ID" not in struct.columns or "Method" not in struct.columns:
        raise ValueError("CPAD structures must have columns 'PDB-ID' and 'Method' in sheet 'final data'.")

    struct["PDB-ID"] = struct["PDB-ID"].astype(str).str.strip().str.upper()
    struct["Method"] = struct["Method"].astype(str).map(norm)

    # add unique methods from structure file directly
    for m in sorted(set(struct["Method"].dropna().tolist())):
        if m and m.lower() not in {"nan", "none"}:
            add_record(records, "CPAD 2.0 structures", m, source_column="Method")

    # map PDB -> methods (some PDB may have multiple method strings)
    pdb_to_methods = (
        struct.loc[struct["PDB-ID"].str.len() == 4, ["PDB-ID", "Method"]]
        .dropna()
        .drop_duplicates()
        .groupby("PDB-ID")["Method"]
        .apply(lambda s: sorted(set(s.tolist())))
        .to_dict()
    )

    # --- peptides ---
    pep = pd.read_excel(peptides_xlsx, sheet_name=CPAD_PEPTIDE_SHEET)
    if "Amyloid Structure" not in pep.columns:
        raise ValueError("CPAD peptides must have column 'Amyloid Structure' in sheet 'peptide'.")

    # Detect native predictor columns (based on what we observed in your file)
    predictor_cols = [
        "% Helix (PASTA 2.0)",
        "% Beta Strand (PASTA 2.0)",
        "% Coil (PASTA 2.0)",
        "% Disorder (PASTA 2.0)",
        "Best Energy Score (PASTA 2.0)",
        "Aggregate Orientation (PASTA 2.0)",
        "Normalized Aggregation Propensity (AGGRESCAN)",
        "Area of the profile Above Threshold (AGGRESCAN):",
        "NuAPRpred",
        "Tango",
    ]
    present = [c for c in predictor_cols if c in pep.columns]

    # Add computational “methods” by native column families
    if any("PASTA 2.0" in c for c in present):
        add_record(records, "CPAD 2.0 peptides", "PASTA 2.0", source_column="; ".join([c for c in present if "PASTA 2.0" in c]))
    if any("AGGRESCAN" in c for c in present):
        add_record(records, "CPAD 2.0 peptides", "AGGRESCAN", source_column="; ".join([c for c in present if "AGGRESCAN" in c]))
    if "NuAPRpred" in present:
        add_record(records, "CPAD 2.0 peptides", "NuAPRpred", source_column="NuAPRpred")
    if "Tango" in present:
        add_record(records, "CPAD 2.0 peptides", "Tango", source_column="Tango")

    # Link peptide PDBs to structure methods
    all_pdbs: Set[str] = set()
    for v in pep["Amyloid Structure"].dropna().astype(str).tolist():
        for pdb in parse_pdb_list(v):
            all_pdbs.add(pdb)

    linked_methods: Set[str] = set()
    for pdb in sorted(all_pdbs):
        for m in pdb_to_methods.get(pdb, []):
            if m and m.lower() not in {"nan", "none"}:
                linked_methods.add(m)

    for m in sorted(linked_methods):
        add_record(records, "CPAD 2.0 (linked structures via peptides)", m,
                   source_column="peptides.Amyloid Structure -> structures.Method")

    return records


# =========================
# MAIN
# =========================

def main():
    # check inputs
    for name, p in INPUT_FILES.items():
        if not os.path.exists(p):
            print(f"ERROR: missing input file for {name}: {p}", file=sys.stderr)
            sys.exit(1)

    records: List[Dict] = []
    records.extend(extract_waltz_methods(INPUT_FILES["WALTZ-DB 2.0"]))
    records.extend(extract_crossbeta_methods(INPUT_FILES["Cross-Beta DB"]))
    records.extend(extract_amyload_methods(INPUT_FILES["AmyLoad"]))
    records.extend(extract_cpad_methods(INPUT_FILES["CPAD 2.0 peptides"], INPUT_FILES["CPAD 2.0 structures"]))

    df = pd.DataFrame(records)
    if df.empty:
        print("No methods extracted. Check source files/columns.", file=sys.stderr)
        sys.exit(2)

    # Deduplicate at the level you asked: per database + raw method name
    df["raw_name_of_method"] = df["raw_name_of_method"].astype(str).str.strip()
    df = df.drop_duplicates(subset=["database_name", "raw_name_of_method"], keep="first").copy()

    # Sort for readability
    df = df.sort_values(["database_name", "tier", "new_name_universal", "raw_name_of_method"]).reset_index(drop=True)

    # Save full with debug columns
    df.to_csv(OUT_CSV_FULL, index=False)

    # Save core 4 columns (exactly as you requested)
    core = df[["database_name", "raw_name_of_method", "new_name_universal", "confidence_level"]].copy()
    core.to_csv(OUT_CSV_CORE, index=False)

    print(f"Saved FULL mapping table: {OUT_CSV_FULL}  (rows={len(df)})")
    print(f"Saved CORE mapping table: {OUT_CSV_CORE}  (rows={len(core)})")
    print("\nPreview (CORE):")
    print(core.head(40).to_string(index=False))


if __name__ == "__main__":
    main()
