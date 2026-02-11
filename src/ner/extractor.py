"""
Biomedical NER Module (CPU-friendly)
Supports:
- scispaCy: Fast, lightweight (default)
- PubMedBERT: Higher quality, slower on CPU
"""

import re
from typing import List, Dict, Set, Tuple, Optional, Union
from dataclasses import dataclass, field
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class BioEntity:
    """Extracted biomedical entity"""
    text: str
    label: str  # PROTEIN, GENE, CHEMICAL, DISEASE, ORGANISM
    start: int
    end: int
    normalized_id: Optional[str] = None  # UniProt, MESH, etc.
    confidence: float = 1.0


@dataclass
class NERResult:
    """NER extraction result for a document"""
    proteins: List[BioEntity] = field(default_factory=list)
    genes: List[BioEntity] = field(default_factory=list)
    chemicals: List[BioEntity] = field(default_factory=list)
    diseases: List[BioEntity] = field(default_factory=list)
    organisms: List[BioEntity] = field(default_factory=list)
    
    def all_entities(self) -> List[BioEntity]:
        return self.proteins + self.genes + self.chemicals + self.diseases + self.organisms
    
    def protein_names(self) -> List[str]:
        """Get unique protein names"""
        return list(set(e.text for e in self.proteins))


# =============================================================================
# PubMedBERT NER Extractor
# =============================================================================

class PubMedBERTExtractor:
    """
    NER using PubMedBERT (transformers)
    Higher quality, slower on CPU
    
    d4data/biomedical-ner-all labels:
    - Anatomy, Chemicals, Disease, Gene/Protein, Species
    """
    
    def __init__(
        self, 
        model: str = "d4data/biomedical-ner-all",
        device: str = "cpu",
        batch_size: int = 1
    ):
        self.model_name = model
        self.device = device
        self.batch_size = batch_size
        self.pipeline = None
        self._load_model()
    
    def _load_model(self):
        """Load transformers NER pipeline"""
        try:
            from transformers import pipeline
            
            logger.info(f"Loading {self.model_name} (this may take a moment)...")
            
            self.pipeline = pipeline(
                "ner",
                model=self.model_name,
                tokenizer=self.model_name,
                aggregation_strategy="simple",
                device=-1 if self.device == "cpu" else 0
            )
            
            logger.info(f"Loaded NER model: {self.model_name}")
            
        except ImportError:
            logger.error("transformers not installed. Run: pip install transformers torch")
            self.pipeline = None
        except Exception as e:
            logger.error(f"Failed to load model: {e}")
            self.pipeline = None
    
    def extract(self, text: str) -> NERResult:
        """Extract entities using PubMedBERT"""
        if self.pipeline is None:
            logger.warning("PubMedBERT not available, using fallback")
            return self._fallback_extract(text)
        
        result = NERResult()
        
        # Split long texts (BERT has 512 token limit)
        chunks = self._split_text(text, max_length=400)
        
        for chunk_start, chunk in chunks:
            try:
                entities = self.pipeline(chunk)
                
                for ent in entities:
                    entity = BioEntity(
                        text=ent['word'],
                        label=ent['entity_group'],
                        start=chunk_start + ent['start'],
                        end=chunk_start + ent['end'],
                        confidence=ent['score']
                    )
                    
                    # Categorize - d4data model uses these labels:
                    # Anatomy, Chemicals, Disease, Gene/Protein, Species
                    label = ent['entity_group'].upper()
                    
                    # FIX: Handle 'GENE/PROTEIN' label from d4data model
                    if label in ('PROTEIN', 'GENE/PROTEIN', 'GENE_OR_GENE_PRODUCT', 'GGP'):
                        result.proteins.append(entity)
                    elif label in ('GENE', 'DNA', 'RNA'):
                        result.genes.append(entity)
                    elif label in ('CHEMICAL', 'CHEMICALS', 'DRUG', 'SIMPLE_CHEMICAL'):
                        result.chemicals.append(entity)
                    elif label in ('DISEASE', 'DISORDER'):
                        result.diseases.append(entity)
                    elif label in ('ORGANISM', 'SPECIES', 'TAXON'):
                        result.organisms.append(entity)
                        
            except Exception as e:
                logger.warning(f"NER failed on chunk: {e}")
                continue
        
        logger.debug(f"PubMedBERT found {len(result.proteins)} proteins")
        return result
    
    def _split_text(self, text: str, max_length: int = 400) -> List[Tuple[int, str]]:
        """Split text into chunks preserving sentence boundaries"""
        if len(text) <= max_length:
            return [(0, text)]
        
        chunks = []
        sentences = re.split(r'(?<=[.!?])\s+', text)
        
        current_chunk = ""
        current_start = 0
        
        for sentence in sentences:
            if len(current_chunk) + len(sentence) > max_length:
                if current_chunk:
                    chunks.append((current_start, current_chunk.strip()))
                current_start += len(current_chunk)
                current_chunk = sentence + " "
            else:
                current_chunk += sentence + " "
        
        if current_chunk.strip():
            chunks.append((current_start, current_chunk.strip()))
        
        return chunks
    
    def _looks_like_protein(self, text: str) -> bool:
        """Heuristic check if text looks like a protein name"""
        patterns = [
            r'^[A-Z][a-z]+[A-Z]',  # CamelCase
            r'^[A-Z]{2,}[0-9]*$',  # HSP70
            r'ase$|in$',          # Enzyme/protein suffix
        ]
        return any(re.search(p, text) for p in patterns)
    
    def _fallback_extract(self, text: str) -> NERResult:
        """Regex fallback when model unavailable"""
        logger.debug("Using fallback regex extraction")
        result = NERResult()
        
        # Specific amyloid-related protein patterns
        protein_patterns = [
            r'\b((?:alpha|α|beta|β|gamma|γ)[- ]?(?:synuclein|amyloid|casein|crystallin))\b',
            r'\b(A[βb](?:40|42|38)?(?:\(1-4[02]\))?)\b',
            r'\b(amyloid[- ]?(?:beta|β)(?:[- ]?\(1-4[02]\))?)\b',
            r'\b([Tt]au(?:[- ]?(?:protein|P301[LS]|441|383))?)\b',
            r'\b(SNCA|APP|PSEN[12]|MAPT|HTT|SOD1|FUS|TARDBP|C9orf72)\b',
            r'\b(TDP-?43|PrP(?:Sc|C)?|TTR|IAPP|SAA[12]?|ATTR)\b',
            r'\b(HSP(?:70|90|40|27)|HSC70)\b',
            r'\b(BACE[12]?|GSK3[βB]?|CDK5)\b',
            r'\b(huntingtin|transthyretin|prion|synuclein|lysozyme|insulin|amylin)\b',
            r'\b(ubiquitin|parkin|LRRK2|PINK1|DJ-1)\b',
            r'\b(apolipoprotein\s*E|ApoE[234]?)\b',
        ]
        
        skip_words = {
            'the', 'and', 'for', 'with', 'from', 'that', 'this', 'these', 'those',
            'method', 'study', 'result', 'data', 'figure', 'protein', 'peptide',
            'fibril', 'aggregate', 'oligomer', 'amyloid', 'plaque', 'tangle',
            'disease', 'patient', 'cell', 'brain', 'analysis', 'structure',
        }
        
        seen = set()
        for pattern in protein_patterns:
            for match in re.finditer(pattern, text, re.IGNORECASE):
                name = match.group(1) if match.lastindex else match.group(0)
                name = name.strip()
                name_lower = name.lower()
                if name_lower not in skip_words and name_lower not in seen and len(name) > 2:
                    seen.add(name_lower)
                    result.proteins.append(BioEntity(
                        text=name, label='PROTEIN',
                        start=match.start(), end=match.end()
                    ))
        return result


# =============================================================================
# scispaCy NER Extractor
# =============================================================================

class BioNERExtractor:
    """
    Biomedical Named Entity Recognition using scispaCy
    CPU-optimized, fast
    """
    
    def __init__(self, model: str = "en_core_sci_sm", use_linker: bool = False):
        self.model_name = model
        self.use_linker = use_linker
        self.nlp = None
        self._load_model()
    
    def _load_model(self):
        """Load spaCy model"""
        try:
            import spacy
            self.nlp = spacy.load(self.model_name)
            logger.info(f"Loaded scispaCy model: {self.model_name}")
            
            if self.use_linker:
                self._add_entity_linker()
                
        except ImportError:
            logger.warning("spaCy not installed, will use regex fallback")
            self.nlp = None
        except OSError:
            logger.warning(f"Model not found: {self.model_name}")
            logger.info("Install with:")
            logger.info(f"  pip install scispacy")
            logger.info(f"  pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/{self.model_name}-0.5.4.tar.gz")
            self.nlp = None
    
    def _add_entity_linker(self):
        """Add UMLS entity linker"""
        try:
            from scispacy.linking import EntityLinker
            self.nlp.add_pipe("scispacy_linker", config={
                "resolve_abbreviations": True,
                "linker_name": "umls"
            })
            logger.info("Entity linker added")
        except Exception as e:
            logger.warning(f"Could not add entity linker: {e}")
    
    def extract(self, text: str) -> NERResult:
        """Extract biomedical entities from text"""
        if self.nlp is None:
            return self._fallback_extract(text)
        
        doc = self.nlp(text)
        result = NERResult()
        
        for ent in doc.ents:
            entity = BioEntity(
                text=ent.text,
                label=ent.label_,
                start=ent.start_char,
                end=ent.end_char
            )
            
            if self.use_linker and hasattr(ent, '_') and hasattr(ent._, 'kb_ents'):
                if ent._.kb_ents:
                    entity.normalized_id = ent._.kb_ents[0][0]
            
            label_upper = ent.label_.upper()
            
            if label_upper in ('PROTEIN', 'GENE_OR_GENE_PRODUCT', 'GGP'):
                result.proteins.append(entity)
            elif label_upper in ('GENE', 'DNA', 'RNA'):
                result.genes.append(entity)
            elif label_upper in ('CHEMICAL', 'SIMPLE_CHEMICAL', 'CHEBI'):
                result.chemicals.append(entity)
            elif label_upper in ('DISEASE', 'DISORDER'):
                result.diseases.append(entity)
            elif label_upper in ('ORGANISM', 'SPECIES', 'TAXON'):
                result.organisms.append(entity)
            elif label_upper == 'ENTITY':
                if self._looks_like_protein(ent.text):
                    result.proteins.append(entity)
        
        return result
    
    def _looks_like_protein(self, text: str) -> bool:
        """Heuristic check if text looks like a protein name"""
        patterns = [
            r'^[A-Z][a-z]+[A-Z]',
            r'^[A-Z]{2,}[0-9]*$',
            r'ase$',
            r'in$',
            r'^(alpha|beta|gamma)-',
        ]
        return any(re.search(p, text) for p in patterns)
    
    def _fallback_extract(self, text: str) -> NERResult:
        """Regex-based fallback when scispaCy not available"""
        logger.debug("Using fallback regex extraction")
        result = NERResult()
        
        # Specific protein patterns only - no broad patterns
        protein_patterns = [
            # Greek prefix proteins
            r'\b((?:alpha|α|beta|β|gamma|γ)[- ]?synuclein)\b',
            r'\b((?:alpha|β|β)[- ]?amyloid)\b',
            # Aβ variants
            r'\b(A[βb](?:40|42|38)?(?:\s*\(1-4[02]\))?)\b',
            r'\b(amyloid[- ]?beta(?:[- ]?\(1-4[02]\))?)\b',
            r'\b(amyloid[- ]?β(?:\s*\(1-4[02]\))?)\b',
            # Tau variants
            r'\b([Tt]au(?:[- ]?(?:protein|P301[LS]|441|383))?)\b',
            # Specific gene/protein names (uppercase)
            r'\b(SNCA|APP|PSEN[12]|MAPT|HTT|mHTT|SOD1|FUS|TARDBP|C9orf72)\b',
            r'\b(TDP-?43|PrP(?:Sc|C)?|TTR|IAPP|SAA[12]?|ATTR)\b',
            r'\b(HSP(?:70|90|40|27)|HSC70|GRP78)\b',
            r'\b(BACE[12]?|GSK3[βBb]?|CDK5|TREM2|CLU)\b',
            # Common amyloid proteins (lowercase)
            r'\b(huntingtin|transthyretin|prion|synuclein|lysozyme|insulin|amylin)\b',
            r'\b(ubiquitin|parkin|LRRK2|PINK1|DJ-?1)\b',
            r'\b(apolipoprotein\s*E|ApoE[234]?)\b',
            r'\b(α-?syn|β-?amy)\b',
            # Specific proteins with -in/-ase suffix that are amyloid-related
            r'\b(cystatin|gelsolin|prolactin|calcitonin|fibrinogen)\b',
            r'\b(lactoferrin|immunoglobulin|keratin|semenogelin)\b',
        ]
        
        seen = set()
        for pattern in protein_patterns:
            for match in re.finditer(pattern, text, re.IGNORECASE):
                name = match.group(1) if match.lastindex else match.group(0)
                name = name.strip()
                name_lower = name.lower()
                
                if len(name) < 2 or name_lower in seen:
                    continue
                
                seen.add(name_lower)
                result.proteins.append(BioEntity(
                    text=name,
                    label='PROTEIN',
                    start=match.start(),
                    end=match.end()
                ))
        
        logger.debug(f"Fallback found {len(result.proteins)} proteins")
        return result


# =============================================================================
# Aggregation Context Extractor
# =============================================================================

class AggregationContextExtractor:
    """Extract aggregation-related context from text"""
    
    AGGREGATION_PATTERNS = {
        'functional': [
            r'functional\s+amyloid',
            r'protective\s+aggregat',
            r'beneficial\s+fibril',
            r'physiological\s+aggregat',
            r'native\s+fibril',
        ],
        'pathological': [
            r'toxic\s+(?:oligomer|aggregate|fibril)',
            r'pathological\s+aggregat',
            r'disease[- ]?associated\s+fibril',
            r'neurotoxic',
            r'cytotoxic\s+aggregat',
        ],
        'experimental_methods': [
            r'thioflavin\s*[- ]?[TtSs]',
            r'ThT\s+(?:fluorescence|assay|binding)',
            r'congo\s*red',
            r'cryo[- ]?EM',
            r'solid[- ]?state\s+NMR',
            r'X[- ]?ray\s+(?:diffraction|fiber)',
            r'AFM|atomic\s+force\s+microscop',
            r'transmission\s+electron\s+microscop|TEM',
            r'circular\s+dichroism|CD\s+spectro',
            r'mass\s+spectrometry|MS[/ ]MS',
            r'SAXS|small[- ]?angle',
            r'fluorescence\s+microscop',
            r'Western\s+blot',
            r'ELISA',
            r'immunohistochem',
            r'pull[- ]?down\s+assay',
            r'cross[- ]?link',
        ],
        'aggregation_type': [
            r'amyloid\s+fibril',
            r'oligomer',
            r'protofibril',
            r'aggregate',
            r'inclusion\s+bod',
            r'plaques?',
            r'tangle',
            r'deposit',
        ],
    }
    
    @classmethod
    def extract_context(cls, text: str, protein_name: str = None) -> Dict[str, List[str]]:
        """Extract aggregation context from text, optionally focused on a specific protein"""
        results = {}
        
        # If protein_name provided, look for context near that protein
        search_text = text
        if protein_name:
            # Find sentences containing the protein
            import re
            sentences = re.split(r'(?<=[.!?])\s+', text)
            relevant = [s for s in sentences if protein_name.lower() in s.lower()]
            if relevant:
                search_text = ' '.join(relevant)
        
        for category, patterns in cls.AGGREGATION_PATTERNS.items():
            matches = []
            for pattern in patterns:
                for match in re.finditer(pattern, search_text, re.IGNORECASE):
                    matches.append(match.group(0))
            results[category] = list(set(matches))
        
        # Add boolean flags for easier processing
        results['forms_aggregates'] = bool(results.get('aggregation_type'))
        results['has_experimental_evidence'] = bool(results.get('experimental_methods'))
        
        return results
        
        return results
    
    @classmethod
    def classify_aggregation_type(cls, text: str) -> str:
        """Classify whether aggregation is functional or pathological"""
        functional_count = 0
        pathological_count = 0
        
        for pattern in cls.AGGREGATION_PATTERNS['functional']:
            functional_count += len(re.findall(pattern, text, re.IGNORECASE))
        
        for pattern in cls.AGGREGATION_PATTERNS['pathological']:
            pathological_count += len(re.findall(pattern, text, re.IGNORECASE))
        
        if functional_count > pathological_count:
            return 'functional'
        elif pathological_count > functional_count:
            return 'pathological'
        else:
            return 'unknown'


# =============================================================================
# Factory function
# =============================================================================

def create_ner_extractor(
    model_type: str = "scispacy",
    use_gpu: bool = False
) -> Union[BioNERExtractor, PubMedBERTExtractor]:
    """
    Factory to create appropriate NER extractor
    
    Args:
        model_type: "scispacy" (fast) or "pubmedbert" (quality)
        use_gpu: Whether GPU is available
        
    Returns:
        NER extractor instance
    """
    if model_type == "pubmedbert":
        device = "cuda" if use_gpu else "cpu"
        return PubMedBERTExtractor(
            model="d4data/biomedical-ner-all",
            device=device
        )
    else:  # scispacy (default)
        if use_gpu:
            return BioNERExtractor(
                model="en_ner_bionlp13cg_md",
                use_linker=True
            )
        else:
            return BioNERExtractor(
                model="en_core_sci_sm",
                use_linker=False
            )


if __name__ == "__main__":
    # Test
    test_text = """
    Amyloid-β (Aβ) protein aggregates are the main component of senile plaques 
    in Alzheimer's disease. Alpha-synuclein forms toxic oligomers in Parkinson's.
    HSP70 and tau protein are also involved. ThT fluorescence shows fibril formation.
    """
    
    print("Testing PubMedBERT extractor...")
    extractor = create_ner_extractor(model_type="pubmedbert")
    result = extractor.extract(test_text)
    
    print(f"\nFound {len(result.proteins)} proteins:")
    for p in result.proteins:
        print(f"  - {p.text} ({p.label}, conf={p.confidence:.2f})")
    
    print(f"\nFound {len(result.chemicals)} chemicals:")
    for c in result.chemicals:
        print(f"  - {c.text}")
    
    print(f"\nFound {len(result.diseases)} diseases:")
    for d in result.diseases:
        print(f"  - {d.text}")
