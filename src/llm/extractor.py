"""
LLM-based Entity Extraction Module
Extracts structured information about aggregating proteins from paper abstracts
"""

import json
import re
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict
from enum import Enum
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class AggregateStructure(Enum):
    BETA_SHEET = "beta-sheet"
    ALPHA_HELICAL = "alpha-helical"
    MIXED = "mixed"
    AMORPHOUS = "amorphous"
    UNKNOWN = "unknown"


class Reversibility(Enum):
    REVERSIBLE = "reversible"
    IRREVERSIBLE = "irreversible"
    CONDITIONALLY_REVERSIBLE = "conditionally_reversible"
    UNKNOWN = "unknown"


class AssemblyMechanism(Enum):
    AUTOCATALYTIC = "autocatalytic"
    COFACTOR_DEPENDENT = "cofactor_dependent"
    MIXED = "mixed"
    UNKNOWN = "unknown"


class FunctionalClass(Enum):
    FUNCTIONAL = "functional"
    PATHOLOGICAL = "pathological"
    BOTH = "both"  # some proteins have dual roles
    UNKNOWN = "unknown"


@dataclass
class ExtractedProtein:
    """Structured extraction result for a protein"""
    protein_name: str
    gene_name: Optional[str] = None
    uniprot_id: Optional[str] = None
    organism: Optional[str] = None
    
    # Aggregation properties
    forms_aggregates: bool = True
    aggregate_structure: str = "unknown"
    reversibility: str = "unknown"
    assembly_mechanism: str = "unknown"
    
    # Classification
    functional_class: str = "unknown"
    biological_role: Optional[str] = None
    associated_disease: Optional[str] = None
    
    # Experimental evidence
    experimental_methods: List[str] = None
    evidence_strength: str = "unknown"  # high, medium, low
    
    # Cofactors/helpers
    assembly_cofactors: List[str] = None
    disassembly_cofactors: List[str] = None
    
    # Source
    pmid: str = ""
    extraction_confidence: float = 0.0
    
    def __post_init__(self):
        if self.experimental_methods is None:
            self.experimental_methods = []
        if self.assembly_cofactors is None:
            self.assembly_cofactors = []
        if self.disassembly_cofactors is None:
            self.disassembly_cofactors = []
    
    def to_dict(self) -> Dict:
        return asdict(self)


# =============================================================================
# PROMPT TEMPLATES
# =============================================================================

EXTRACTION_SYSTEM_PROMPT = """You are an expert bioinformatics assistant specialized in protein aggregation and amyloid biology.

Your task is to extract structured information about proteins that form aggregates from scientific paper abstracts.

IMPORTANT DEFINITIONS:
- We use a BROAD definition of amyloidogenic proteins: any protein capable of forming aggregates
- This includes beta-sheet amyloids, alpha-helical amyloids, and other aggregate types
- We are NOT interested in prions or prion-like proteins - SKIP these entirely
- We want both functional (beneficial) and pathological (disease-related) aggregating proteins

EXTRACTION RULES:
1. Only extract proteins with EXPERIMENTAL EVIDENCE of aggregation
2. Do not infer or guess - only extract what is explicitly stated
3. If information is not available, use "unknown" or null
4. Extract ALL proteins mentioned that meet criteria, not just the main subject
5. Be precise with protein names - use standard nomenclature when possible
6. For UniProt IDs: only include if explicitly mentioned or clearly identifiable

OUTPUT FORMAT:
Return a JSON object with a "proteins" array containing extracted protein information.
"""

EXTRACTION_USER_PROMPT = """Analyze this scientific abstract and extract information about proteins that form aggregates.

ABSTRACT:
Title: {title}
Text: {abstract}
PMID: {pmid}

Extract the following for EACH aggregating protein mentioned:

1. **Protein identification**:
   - protein_name: Standard protein name
   - gene_name: Gene symbol (if mentioned)
   - uniprot_id: UniProt accession (if mentioned)
   - organism: Species (if mentioned)

2. **Aggregation properties**:
   - forms_aggregates: true/false
   - aggregate_structure: "beta-sheet", "alpha-helical", "mixed", "amorphous", or "unknown"
   - reversibility: "reversible", "irreversible", "conditionally_reversible", or "unknown"
   - assembly_mechanism: "autocatalytic", "cofactor_dependent", "mixed", or "unknown"

3. **Classification**:
   - functional_class: "functional", "pathological", "both", or "unknown"
   - biological_role: Brief description of biological function (if functional)
   - associated_disease: Disease name (if pathological)

4. **Evidence**:
   - experimental_methods: List of methods used (e.g., ["ThT fluorescence", "TEM", "CD spectroscopy"])
   - evidence_strength: "high" (multiple methods), "medium" (single method), "low" (indirect)

5. **Cofactors/helpers** (proteins, RNAs, or small molecules that assist):
   - assembly_cofactors: List of factors that help aggregation
   - disassembly_cofactors: List of factors that help disaggregation

CRITICAL:
- SKIP any proteins related to prions (PrP, prion, prion-like)
- Only include proteins with experimental evidence of aggregation
- Set extraction_confidence (0.0-1.0) based on how clear the information is

Return JSON:
{{
  "proteins": [
    {{
      "protein_name": "...",
      "gene_name": "...",
      ...
    }}
  ],
  "paper_relevance": "high/medium/low",
  "notes": "any important observations"
}}
"""


class LLMExtractor:
    """
    Extracts structured protein information using LLM
    Supports Anthropic, OpenAI, and Ollama (local) backends
    """
    
    def __init__(
        self, 
        provider: str = "anthropic",
        model: str = None,
        api_key: str = None,
        ollama_base_url: str = "http://localhost:11434"
    ):
        self.provider = provider
        self.api_key = api_key
        self.ollama_base_url = ollama_base_url
        
        if provider == "anthropic":
            self.model = model or "claude-sonnet-4-20250514"
            self._init_anthropic()
        elif provider == "openai":
            self.model = model or "gpt-4-turbo-preview"
            self._init_openai()
        elif provider == "ollama":
            # Good models for extraction: llama3.1, mistral, mixtral, qwen2.5
            self.model = model or "llama3.1:8b"
            self._init_ollama()
        else:
            raise ValueError(f"Unknown provider: {provider}. Use: anthropic, openai, ollama")
    
    def _init_anthropic(self):
        """Initialize Anthropic client"""
        try:
            import anthropic
            self.client = anthropic.Anthropic(api_key=self.api_key)
        except ImportError:
            logger.warning("anthropic package not installed. Run: pip install anthropic")
            self.client = None
    
    def _init_openai(self):
        """Initialize OpenAI client"""
        try:
            import openai
            self.client = openai.OpenAI(api_key=self.api_key)
        except ImportError:
            logger.warning("openai package not installed. Run: pip install openai")
            self.client = None
    
    def _init_ollama(self):
        """Initialize Ollama client (local LLM)"""
        # Ollama uses REST API, no special client needed
        self.client = True  # Mark as initialized
        
        # Check if Ollama is running
        try:
            import requests
            response = requests.get(f"{self.ollama_base_url}/api/tags", timeout=5)
            if response.status_code == 200:
                models = [m['name'] for m in response.json().get('models', [])]
                if models:
                    logger.info(f"Ollama available. Models: {models}")
                    if self.model not in models and not any(self.model in m for m in models):
                        logger.warning(f"Model {self.model} not found. Available: {models}")
                        logger.info(f"Pull it with: ollama pull {self.model}")
            else:
                logger.warning("Ollama server not responding")
                self.client = None
        except Exception as e:
            logger.warning(f"Ollama not available: {e}")
            logger.info("Install Ollama: https://ollama.ai/download")
            self.client = None
    
    def extract_from_abstract(
        self, 
        title: str, 
        abstract: str, 
        pmid: str
    ) -> Dict[str, Any]:
        """
        Extract protein information from a single abstract
        
        Args:
            title: Paper title
            abstract: Paper abstract
            pmid: PubMed ID
            
        Returns:
            Dictionary with extracted proteins and metadata
        """
        if not abstract or len(abstract) < 50:
            return {"proteins": [], "paper_relevance": "low", "error": "Abstract too short"}
        
        prompt = EXTRACTION_USER_PROMPT.format(
            title=title,
            abstract=abstract,
            pmid=pmid
        )
        
        try:
            if self.provider == "anthropic":
                response = self._call_anthropic(prompt)
            elif self.provider == "openai":
                response = self._call_openai(prompt)
            elif self.provider == "ollama":
                response = self._call_ollama(prompt)
            else:
                return {"proteins": [], "error": f"Unknown provider: {self.provider}"}
            
            # Parse JSON response
            result = self._parse_response(response, pmid)
            return result
            
        except Exception as e:
            logger.error(f"Extraction failed for PMID {pmid}: {e}")
            return {"proteins": [], "paper_relevance": "unknown", "error": str(e)}
    
    def _call_anthropic(self, prompt: str) -> str:
        """Call Anthropic API"""
        if self.client is None:
            raise RuntimeError("Anthropic client not initialized")
        
        message = self.client.messages.create(
            model=self.model,
            max_tokens=4096,
            temperature=0.1,
            system=EXTRACTION_SYSTEM_PROMPT,
            messages=[
                {"role": "user", "content": prompt}
            ]
        )
        return message.content[0].text
    
    def _call_openai(self, prompt: str) -> str:
        """Call OpenAI API"""
        if self.client is None:
            raise RuntimeError("OpenAI client not initialized")
        
        response = self.client.chat.completions.create(
            model=self.model,
            temperature=0.1,
            messages=[
                {"role": "system", "content": EXTRACTION_SYSTEM_PROMPT},
                {"role": "user", "content": prompt}
            ]
        )
        return response.choices[0].message.content
    
    def _call_ollama(self, prompt: str) -> str:
        """Call Ollama local API"""
        import requests
        
        if self.client is None:
            raise RuntimeError("Ollama not available. Install: https://ollama.ai/download")
        
        response = requests.post(
            f"{self.ollama_base_url}/api/generate",
            json={
                "model": self.model,
                "prompt": f"{EXTRACTION_SYSTEM_PROMPT}\n\n{prompt}",
                "stream": False,
                "options": {
                    "temperature": 0.1,
                    "num_predict": 4096,
                }
            },
            timeout=120  # Local models can be slow
        )
        response.raise_for_status()
        return response.json().get("response", "")
    
    def _parse_response(self, response: str, pmid: str) -> Dict[str, Any]:
        """Parse LLM response and validate structure"""
        
        # Try to extract JSON from response
        json_match = re.search(r'\{[\s\S]*\}', response)
        if not json_match:
            return {"proteins": [], "error": "No JSON found in response"}
        
        try:
            data = json.loads(json_match.group())
        except json.JSONDecodeError as e:
            return {"proteins": [], "error": f"JSON parse error: {e}"}
        
        # Validate and clean proteins
        proteins = []
        for p in data.get("proteins", []):
            # Skip prion-related (safety check)
            name_lower = p.get("protein_name", "").lower()
            if any(term in name_lower for term in ["prion", "prp"]):
                continue
            
            # Add PMID
            p["pmid"] = pmid
            
            # Validate required fields
            if p.get("protein_name"):
                # Create ExtractedProtein to validate structure
                try:
                    protein = ExtractedProtein(
                        protein_name=p.get("protein_name"),
                        gene_name=p.get("gene_name"),
                        uniprot_id=p.get("uniprot_id"),
                        organism=p.get("organism"),
                        forms_aggregates=p.get("forms_aggregates", True),
                        aggregate_structure=p.get("aggregate_structure", "unknown"),
                        reversibility=p.get("reversibility", "unknown"),
                        assembly_mechanism=p.get("assembly_mechanism", "unknown"),
                        functional_class=p.get("functional_class", "unknown"),
                        biological_role=p.get("biological_role"),
                        associated_disease=p.get("associated_disease"),
                        experimental_methods=p.get("experimental_methods", []),
                        evidence_strength=p.get("evidence_strength", "unknown"),
                        assembly_cofactors=p.get("assembly_cofactors", []),
                        disassembly_cofactors=p.get("disassembly_cofactors", []),
                        pmid=pmid,
                        extraction_confidence=p.get("extraction_confidence", 0.5)
                    )
                    proteins.append(protein.to_dict())
                except Exception as e:
                    logger.warning(f"Invalid protein data: {e}")
                    continue
        
        return {
            "proteins": proteins,
            "paper_relevance": data.get("paper_relevance", "unknown"),
            "notes": data.get("notes", "")
        }
    
    def extract_batch(
        self, 
        articles: List[Dict[str, str]],
        progress_callback=None
    ) -> List[Dict[str, Any]]:
        """
        Extract from multiple articles
        
        Args:
            articles: List of {"title": ..., "abstract": ..., "pmid": ...}
            progress_callback: Optional callback(current, total)
            
        Returns:
            List of extraction results
        """
        results = []
        total = len(articles)
        
        for i, article in enumerate(articles):
            result = self.extract_from_abstract(
                title=article.get("title", ""),
                abstract=article.get("abstract", ""),
                pmid=article.get("pmid", "")
            )
            results.append(result)
            
            if progress_callback:
                progress_callback(i + 1, total)
        
        return results


# =============================================================================
# Rule-based pre-filtering (to reduce LLM calls)
# =============================================================================

class RuleBasedFilter:
    """
    Pre-filter abstracts using keyword rules before LLM extraction
    Reduces API costs by filtering obviously irrelevant papers
    """
    
    # Keywords indicating experimental aggregation evidence
    POSITIVE_INDICATORS = [
        r'\b(amyloid|fibril|aggregate|oligomer)\b',
        r'\b(ThT|thioflavin|congo red)\b',
        r'\b(electron microscopy|TEM|AFM|cryo-EM)\b',
        r'\b(aggregation kinetics|fibril formation)\b',
        r'\b(cross-beta|beta-sheet)\b',
        r'\b(self-assembl|self assembl)\b',
    ]
    
    # Keywords to exclude (prions)
    NEGATIVE_INDICATORS = [
        r'\bprion\b',
        r'\bPrP\b',
        r'\bprion-like\b',
        r'\btransmissible spongiform\b',
    ]
    
    @classmethod
    def should_process(cls, title: str, abstract: str) -> tuple[bool, float]:
        """
        Determine if abstract should be sent to LLM
        
        Returns:
            (should_process, confidence_score)
        """
        text = f"{title} {abstract}".lower()
        
        # Check negative indicators first
        for pattern in cls.NEGATIVE_INDICATORS:
            if re.search(pattern, text, re.IGNORECASE):
                return False, 0.0
        
        # Count positive indicators
        positive_matches = 0
        for pattern in cls.POSITIVE_INDICATORS:
            if re.search(pattern, text, re.IGNORECASE):
                positive_matches += 1
        
        # Need at least 2 positive indicators
        if positive_matches >= 2:
            confidence = min(positive_matches / 4.0, 1.0)  # Max at 4 matches
            return True, confidence
        elif positive_matches == 1:
            return True, 0.3  # Low confidence, but worth checking
        else:
            return False, 0.0


if __name__ == "__main__":
    # Test the module with a sample abstract
    test_abstract = """
    Alpha-synuclein aggregation is a hallmark of Parkinson's disease. 
    We investigated the aggregation kinetics of alpha-synuclein using 
    ThT fluorescence assays and transmission electron microscopy (TEM). 
    Our results show that alpha-synuclein forms beta-sheet rich fibrils 
    through a nucleation-dependent mechanism. The aggregation was 
    irreversible under physiological conditions. We also found that 
    Hsp70 chaperone can prevent aggregation when present during the 
    lag phase.
    """
    
    # Test rule-based filter
    should_process, confidence = RuleBasedFilter.should_process(
        "Alpha-synuclein aggregation in Parkinson's disease",
        test_abstract
    )
    print(f"Should process: {should_process}, confidence: {confidence}")
    
    # Test LLM extractor (without actual API call)
    print("\nExtraction prompt would be:")
    print(EXTRACTION_USER_PROMPT.format(
        title="Alpha-synuclein aggregation in Parkinson's disease",
        abstract=test_abstract[:200] + "...",
        pmid="12345678"
    )[:500])
