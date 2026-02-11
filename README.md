# Amyloid Aggregation Database Pipeline

Автоматизированный пайплайн для создания базы данных белков, способных образовывать агрегаты, с использованием LLM для извлечения информации из научных статей.

## Расширенное определение

Мы используем **широкое определение** амилоидогенного белка:
- Любой белок, способный образовывать агрегаты
- Включает β-складчатые, α-спиральные и смешанные структуры
- **Исключены:** прионы и прион-подобные белки

### Типы агрегации

1. **Обратимая (функциональные амилоиды)**
   - Автокаталитическая сборка/разборка
   - С участием кофакторов (шапероны, РНК и др.)

2. **Необратимая (патологические амилоиды)**
   - Автокаталитическая сборка → стабилизация
   - Сборка с кофакторами → стабилизация другими агентами

## Архитектура

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

### Режимы извлечения

| Режим | Описание | GPU | Скорость | Точность |
|-------|----------|-----|----------|----------|
| `hybrid` | NER + Patterns + LLM при необходимости | No | ⚡⚡ | High |
| `ner_only` | Только NER + patterns | No | ⚡⚡⚡ | Moderate |
| `llm_only` | Только LLM (старый режим) | No | ⚡ | High |

### NER модели

| Модель | Скорость (CPU) | RAM | Качество | Когда использовать |
|--------|----------------|-----|----------|-------------------|
| `scispacy` | ~0.1 сек/абстракт | ~200MB | Moderate | Большие датасеты (>500 статей) |
| `pubmedbert` | ~3-5 сек/абстракт | ~2-4GB | High | Качественная обработка (<200 статей) |

## Установка

```bash
cd amyloid_pipeline
pip install -r requirements.txt
```

## Использование

### Получение API ключей

**Anthropic (Claude):**
1. https://console.anthropic.com → Sign up
2. Settings → API Keys → Create Key

**OpenAI (GPT-4):**
1. https://platform.openai.com → Sign up
2. API Keys → Create new secret key

**Ollama (бесплатно, локально):**
- API ключ не нужен!
- Установка: https://ollama.ai/download

### Запуск полного пайплайна

```bash
# РЕКОМЕНДУЕМЫЙ: Гибридный режим с scispaCy (быстро)
python src/pipeline.py \
    --email your@email.com \
    --mode hybrid \
    --ner-model scispacy \
    --max-papers 500

# Качественная обработка с PubMedBERT (медленнее, точнее)
python src/pipeline.py \
    --email your@email.com \
    --mode hybrid \
    --ner-model pubmedbert \
    --max-papers 100

# Только NER без LLM (самый быстрый)
python src/pipeline.py \
    --email your@email.com \
    --mode ner_only \
    --ner-model scispacy \
    --max-papers 2000

# С Anthropic API для LLM
python src/pipeline.py \
    --email your@email.com \
    --mode hybrid \
    --ner-model pubmedbert \
    --llm-provider anthropic \
    --llm-key sk-ant-api03-xxxxx \
    --profile quality
```

### Опции

| Флаг | Описание | По умолчанию |
|------|----------|--------------|
| `--email` | Email для PubMed API (обязателен) | - |
| `--mode` | `hybrid`, `ner_only`, `llm_only` | `hybrid` |
| `--ner-model` | `scispacy` (быстро) или `pubmedbert` (качественно) | `scispacy` |
| `--profile` | `fast`, `balanced`, `quality` | `balanced` |
| `--llm-provider` | `anthropic`, `openai`, `ollama` | `ollama` |
| `--llm-key` | API ключ (не нужен для ollama) | - |
| `--llm-model` | Модель (например, `qwen2.5:7b`) | auto |
| `--max-papers` | Максимум статей | 100 |
| `--date-from` | Начальная дата поиска | 2020/01/01 |
| `--dry-run` | Только поиск, без извлечения | False |

### Рекомендуемые модели Ollama (CPU)

| Модель | RAM | Для чего |
|--------|-----|----------|
| `qwen2.5:3b` | ~4GB | Быстрое тестирование |
| `qwen2.5:7b` | ~8GB | **Рекомендуется** |
| `qwen2.5:14b` | ~16GB | Максимальное качество |
| `llama3.2:3b` | ~4GB | Альтернатива Qwen |

### Dry run (тестирование запроса)

```bash
python src/pipeline.py --email your@email.com --dry-run --max-papers 10
```

## Структура проекта

```
amyloid_pipeline/
├── config/
│   └── settings.py          # Конфигурация, поисковые термины
├── src/
│   ├── pubmed/
│   │   └── searcher.py      # PubMed API клиент
│   ├── llm/
│   │   └── extractor.py     # LLM экстракция сущностей
│   ├── string_db/
│   │   └── client.py        # STRING DB интеграция
│   ├── database/
│   │   └── storage.py       # SQLite хранилище
│   └── pipeline.py          # Главный оркестратор
├── data/
│   ├── raw/                 # Сырые данные
│   └── processed/           # Обработанные данные
└── tests/
```

## База данных

SQLite схема:

- `proteins` — основная информация о белках
- `papers` — метаданные статей
- `protein_papers` — связь белок-статья
- `experimental_methods` — экспериментальные методы
- `cofactors` — кофакторы сборки/разборки
- `string_interactions` — взаимодействия из STRING

### Экспорт в CSV

```python
from src.database.storage import AmyloidDatabase
db = AmyloidDatabase("data/amyloid_aggregation.db")
db.export_to_csv("proteins_export.csv")
```

## Конфигурация поисковых терминов

Редактируйте `config/settings.py`:

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

## Примеры использования

### Поиск функциональных амилоидов

```python
from src.pubmed.searcher import search_functional_amyloids

for article in search_functional_amyloids(max_results=100):
    print(f"{article.pmid}: {article.title}")
```

### Анализ кофакторов

```python
from src.pipeline import AmyloidPipeline

pipeline = AmyloidPipeline()
cofactor_analysis = pipeline.run_cofactor_analysis()

print("Shared cofactors:", cofactor_analysis['shared'])
print("Functional only:", cofactor_analysis['functional_only'])
print("Pathological only:", cofactor_analysis['pathological_only'])
```

## Следующие шаги

- [ ] Добавить поддержку full-text извлечения (PubMed Central)
- [ ] Интеграция с UniProt для валидации ID
- [ ] Веб-интерфейс для просмотра базы
- [ ] Экспорт в формат для AmyloidBench
- [ ] Анализ пересечения кофакторов между функциональными/патологическими

## Лицензия

MIT
