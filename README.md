# Amyloid Aggregation Database Pipeline

Automated pipeline for DB creation of proteins capable to aggregate, using LLM for information mining from published data.

## Defyning therminology

We use **broad therminology** of amyloidogenic protein:
- Any protein able to form aggregates
- Includes β-sheets, α-helix и mixed cross-structures
- **Exclusions:** prion and prion-like proteins

### Types of aggregation

1. **Reversible (functional amyloids)**
   - Self-induced assembly/disassembly
   - Using cofactors (chaperons, RNA, etc.)

2. **Irreversible (pathological amyloids)**
   - Autocatalysed assembly → stabilisation
   - Assembly with cofactors → stabilised by other agents

## Architecture

```
                    ┌─────────────────────────────────────┐
                    │         HYBRID MODE (default)       │
                    │  ┌─────────┐    ┌─────────────────┐ │
PubMed Search ──────┼─▶│scispaCy │───▶│ Pattern Matching│ │
                    │  │  NER    │    │   (context)     │ │
                    │  └─────────┘    └────────┬────────┘ │
                    │                          │          │
                    │         confidence < 0.5?│          │
                    │                    ┌─────▼─────┐    │
                    │                    │ LLM (opt) │    │
                    │                    └───────────┘    │
                    └─────────────────────────────────────┘
                                       │
                                       ▼
                              ┌─────────────────┐
                              │    Database     │
                              └────────┬────────┘
                                       │
                              ┌────────▼────────┐
                              │   STRING DB     │
                              │  (enrichment)   │
                              └─────────────────┘
```

### Mode of extraction

| Mode | Description | GPU | Speed | Accuracy |
|-------|----------|-----|----------|----------|
| `hybrid` | NER + Patterns + LLM if necessary | No | ⚡⚡ | High |
| `ner_only` | Only NER + patterns | No | ⚡⚡⚡ | Moderate |
| `llm_only` | Only LLM (old mode) | No | ⚡ | High |

### NER models

| Model | Speed (CPU) | RAM | Quality | When to use |
|--------|----------------|-----|----------|-------------------|
| `scispacy` | ~0.1 sec/abstract | ~200MB | Moderate | Big datasets (>500 статей) |
| `pubmedbert` | ~3-5 sec/abstract | ~2-4GB | High | High-quality processing (<200 статей) |

## Installation

```bash
cd amyloid_pipeline
pip install -r requirements.txt
```

## Usage

### Get API keys

**Anthropic (Claude):**
1. https://console.anthropic.com → Sign up
2. Settings → API Keys → Create Key

**OpenAI (GPT-4):**
1. https://platform.openai.com → Sign up
2. API Keys → Create new secret key

**Ollama (free, locally):**
- API key not needed!
- Installation: https://ollama.ai/download

### Launch full pipeline

```bash
# RECOMMENDED: Hybrid mode with scispaCy (fast)
python src/pipeline.py \
    --email your@email.com \
    --mode hybrid \
    --ner-model scispacy \
    --max-papers 500

# High-quality processing with PubMedBERT (slower, more accurate)
python src/pipeline.py \
    --email your@email.com \
    --mode hybrid \
    --ner-model pubmedbert \
    --max-papers 100

# Only NER without LLM (the fastest)
python src/pipeline.py \
    --email your@email.com \
    --mode ner_only \
    --ner-model scispacy \
    --max-papers 2000

# With Anthropic API for LLM
python src/pipeline.py \
    --email your@email.com \
    --mode hybrid \
    --ner-model pubmedbert \
    --llm-provider anthropic \
    --llm-key sk-ant-api03-xxxxx \
    --profile quality
```

### Options

| Flag | Description | Default |
|------|----------|--------------|
| `--email` | Email for PubMed API (required) | - |
| `--mode` | `hybrid`, `ner_only`, `llm_only` | `hybrid` |
| `--ner-model` | `scispacy` (fast) или `pubmedbert` (high-quality) | `scispacy` |
| `--profile` | `fast`, `balanced`, `quality` | `balanced` |
| `--llm-provider` | `anthropic`, `openai`, `ollama` | `ollama` |
| `--llm-key` | API key (optional for ollama) | - |
| `--llm-model` | Model (e.g. `qwen2.5:7b`) | auto |
| `--max-papers` | Max articles | 100 |
| `--date-from` | Starting data of search | 2020/01/01 |
| `--dry-run` | Only search, without extraction | False |

### Recommended models Ollama (CPU)

| Model | RAM | Purpose |
|--------|-----|----------|
| `qwen2.5:3b` | ~4GB | Fast testing |
| `qwen2.5:7b` | ~8GB | **Recommended** |
| `qwen2.5:14b` | ~16GB | Maximum of quality |
| `llama3.2:3b` | ~4GB | Qwen alternative |

### Dry run (test request)

```bash
python src/pipeline.py --email your@email.com --dry-run --max-papers 10
```

## Project structure

```
amyloid_pipeline/
├── config/
│   └── settings.py          # Configuration, search patterns
├── src/
│   ├── pubmed/
│   │   └── searcher.py      # PubMed API client
│   ├── llm/
│   │   └── extractor.py     # LLM term mining
│   ├── string_db/
│   │   └── client.py        # STRING DB integration
│   ├── database/
│   │   └── storage.py       # SQLite storage
│   └── pipeline.py          # Main orchestrator
├── data/
│   ├── raw/                 # Raw data
│   └── processed/           # Processed data
└── tests/
```

## Database

SQLite scheme:

- `proteins` — proteins' main information
- `papers` — articles metadata
- `protein_papers` — protein-article linkage
- `experimental_methods` — experimental methods
- `cofactors` — assembly/disassembly cofactors
- `string_interactions` — interactions from STRING

### CSV export

```python
from src.database.storage import AmyloidDatabase
db = AmyloidDatabase("data/amyloid_aggregation.db")
db.export_to_csv("proteins_export.csv")
```

## Search patterns configuration

Edit `config/settings.py`:

```python
AGGREGATION_TERMS = [
    "amyloid",
    "protein aggregation",
    "fibril formation",
    # ...
]

EXCLUSION_TERMS = [
    "prion",
    "prion-like",
    # ...
]
```

## Example of usage

### Search for functional amyloids

```python
from src.pubmed.searcher import search_functional_amyloids

for article in search_functional_amyloids(max_results=100):
    print(f"{article.pmid}: {article.title}")
```

### Cofactors' analysis

```python
from src.pipeline import AmyloidPipeline

pipeline = AmyloidPipeline()
cofactor_analysis = pipeline.run_cofactor_analysis()

print("Shared cofactors:", cofactor_analysis['shared'])
print("Functional only:", cofactor_analysis['functional_only'])
print("Pathological only:", cofactor_analysis['pathological_only'])
```

## TODO:

- [ ] Add full-text extraction support (PubMed Central)
- [ ] Integration with UniProt for ID validation
- [ ] Web-interface for DB search
- [ ] Export into the format for AmyloidBench
- [ ] Analysis of cofactors' intersection between functional/pathological amyloids

## License

MIT
