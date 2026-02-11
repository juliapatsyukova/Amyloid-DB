"""
Hybrid Extraction Pipeline (CPU-optimized)
Combines:
1. scispaCy NER for protein/entity extraction (fast, accurate)
2. Pattern matching for aggregation context
3. LLM only for complex classification tasks

This reduces LLM calls by 70-80% while improving accuracy on biomedical terms.
"""

import json
from typing import List, Dict, Optional, Any
from dataclasses import dataclass, asdict
import logging
import os

import sys
project_dir = os.path.dirname(os.path.abspath("amyloid_pipeline"))
sys.path.append(project_dir)

from src.ner.extractor import (
    BioNERExtractor, 
    PubMedBERTExtractor,
    AggregationContextExtractor, 
    NERResult,
    create_ner_extractor
)
from src.llm.extractor import LLMExtractor

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class HybridExtractionResult:
    """Combined result from NER + pattern matching + LLM"""
    proteins: List[Dict]
    paper_relevance: str
    ner_proteins_count: int
    llm_called: bool
    extraction_method: str  # 'ner_only', 'ner+patterns', 'hybrid'


# Simplified LLM prompt for classification only (not extraction)
CLASSIFICATION_PROMPT = """Analyze these protein aggregation findings and classify them.

PROTEIN: {protein_name}
CONTEXT FROM PAPER:
- Forms aggregates: {forms_aggregates}
- Structure hints: {aggregate_structure}
- Reversibility hints: {reversibility}
- Functional/pathological hints: {functional_class}
- Experimental methods: {methods}
- Potential cofactors: {cofactors}

ABSTRACT EXCERPT: {abstract_excerpt}

Classify and fill gaps. Return JSON only:
{{
  "aggregate_structure": "beta-sheet|alpha-helical|mixed|amorphous|unknown",
  "reversibility": "reversible|irreversible|conditionally_reversible|unknown",
  "assembly_mechanism": "autocatalytic|cofactor_dependent|mixed|unknown",
  "functional_class": "functional|pathological|both|unknown",
  "biological_role": "string or null",
  "associated_disease": "string or null",
  "evidence_strength": "high|medium|low",
  "confidence": 0.0-1.0,
  "assembly_cofactors": ["list"],
  "disassembly_cofactors": ["list"]
}}
"""


class HybridExtractor:
    """
    Hybrid extraction: NER + Patterns + LLM
    
    Strategy:
    1. NER extracts protein names (fast, accurate)
    2. Pattern matching extracts aggregation context
    3. LLM called ONLY when:
       - Classification is ambiguous
       - Multiple proteins need relationship analysis
       - High-value papers with partial info
    """
    
    def __init__(
        self,
        ner_model: str = "scispacy",  # "scispacy" or "pubmedbert"
        llm_provider: str = "ollama",
        llm_model: str = None,
        llm_api_key: str = None,
        ollama_url: str = "http://localhost:11434",
        use_llm: bool = True,
        llm_threshold: float = 0.5  # Call LLM if pattern confidence < threshold
    ):
        # NER extractor
        self.ner_model_type = ner_model
        self.ner = create_ner_extractor(model_type=ner_model, use_gpu=False)
        
        # LLM extractor (optional)
        self.use_llm = use_llm
        self.llm_threshold = llm_threshold
        
        if use_llm:
            self.llm = LLMExtractor(
                provider=llm_provider,
                model=llm_model,
                api_key=llm_api_key,
                ollama_base_url=ollama_url
            )
        else:
            self.llm = None
        
        # Stats
        self.stats = {
            'papers_processed': 0,
            'ner_extractions': 0,
            'pattern_extractions': 0,
            'llm_calls': 0,
            'llm_skipped': 0,
            'ner_model': ner_model,
        }
    
    def extract(
        self,
        title: str,
        abstract: str,
        pmid: str,
        force_llm: bool = False
    ) -> HybridExtractionResult:
        """
        Extract protein aggregation info using hybrid approach
        
        Args:
            title: Paper title
            abstract: Paper abstract
            pmid: PubMed ID
            force_llm: Force LLM call even if patterns are confident
            
        Returns:
            HybridExtractionResult
        """
        self.stats['papers_processed'] += 1
        
        full_text = f"{title}. {abstract}"
        
        # Step 1: NER extraction
        ner_result = self.ner.extract(full_text)
        self.stats['ner_extractions'] += 1
        
        protein_names = ner_result.protein_names()
        
        if not protein_names:
            return HybridExtractionResult(
                proteins=[],
                paper_relevance='low',
                ner_proteins_count=0,
                llm_called=False,
                extraction_method='ner_only'
            )
        
        # Step 2: Pattern-based context extraction for each protein
        proteins = []
        needs_llm = []
        
        for protein_name in protein_names:
            context = AggregationContextExtractor.extract_context(
                full_text, 
                protein_name
            )
            self.stats['pattern_extractions'] += 1
            
            # Calculate pattern confidence
            confidence = self._calculate_pattern_confidence(context)
            
            protein_data = {
                'protein_name': protein_name,
                'pmid': pmid,
                'extraction_confidence': confidence,
                **context
            }
            
            # Decide if LLM needed
            if confidence < self.llm_threshold or force_llm:
                needs_llm.append(protein_data)
            else:
                proteins.append(protein_data)
        
        # Step 3: LLM for ambiguous cases
        llm_called = False
        
        if needs_llm and self.use_llm and self.llm:
            for protein_data in needs_llm:
                enriched = self._enrich_with_llm(protein_data, abstract)
                proteins.append(enriched)
                llm_called = True
                self.stats['llm_calls'] += 1
        elif needs_llm:
            # No LLM available, use pattern results
            proteins.extend(needs_llm)
            self.stats['llm_skipped'] += len(needs_llm)
        
        # Determine paper relevance
        if any(p.get('forms_aggregates') for p in proteins):
            relevance = 'high' if len(proteins) > 1 else 'medium'
        else:
            relevance = 'low'
        
        return HybridExtractionResult(
            proteins=proteins,
            paper_relevance=relevance,
            ner_proteins_count=len(protein_names),
            llm_called=llm_called,
            extraction_method='hybrid' if llm_called else 'ner+patterns'
        )
    
    def _calculate_pattern_confidence(self, context: Dict) -> float:
        """Calculate confidence based on pattern matches"""
        score = 0.0
        
        if context.get('forms_aggregates'):
            score += 0.3
        
        if context.get('aggregate_structure') != 'unknown':
            score += 0.15
        
        if context.get('reversibility') != 'unknown':
            score += 0.15
        
        if context.get('functional_class') != 'unknown':
            score += 0.15
        
        if context.get('experimental_methods'):
            score += min(0.15, len(context['experimental_methods']) * 0.05)
        
        if context.get('cofactors'):
            score += 0.1
        
        return min(score, 1.0)
    
    def _enrich_with_llm(self, protein_data: Dict, abstract: str) -> Dict:
        """Use LLM to fill gaps in pattern extraction"""
        prompt = CLASSIFICATION_PROMPT.format(
            protein_name=protein_data['protein_name'],
            forms_aggregates=protein_data.get('forms_aggregates', 'unknown'),
            aggregate_structure=protein_data.get('aggregate_structure', 'unknown'),
            reversibility=protein_data.get('reversibility', 'unknown'),
            functional_class=protein_data.get('functional_class', 'unknown'),
            methods=', '.join(protein_data.get('experimental_methods', [])) or 'unknown',
            cofactors=protein_data.get('cofactors', []),
            abstract_excerpt=abstract[:500]
        )
        
        try:
            if self.llm.provider == "anthropic":
                response = self.llm._call_anthropic(prompt)
            elif self.llm.provider == "openai":
                response = self.llm._call_openai(prompt)
            elif self.llm.provider == "ollama":
                response = self.llm._call_ollama(prompt)
            else:
                return protein_data
            
            # Parse JSON response
            import re
            json_match = re.search(r'\{[\s\S]*\}', response)
            if json_match:
                llm_data = json.loads(json_match.group())
                
                # Merge LLM results with pattern results (LLM fills gaps)
                for key in ['aggregate_structure', 'reversibility', 'assembly_mechanism',
                           'functional_class', 'biological_role', 'associated_disease',
                           'evidence_strength', 'assembly_cofactors', 'disassembly_cofactors']:
                    if key in llm_data and llm_data[key] and llm_data[key] != 'unknown':
                        if protein_data.get(key) in (None, 'unknown', []):
                            protein_data[key] = llm_data[key]
                
                # Update confidence
                if 'confidence' in llm_data:
                    protein_data['extraction_confidence'] = max(
                        protein_data.get('extraction_confidence', 0),
                        llm_data['confidence']
                    )
                
        except Exception as e:
            logger.warning(f"LLM enrichment failed: {e}")
        
        return protein_data
    
    def get_stats(self) -> Dict:
        """Get extraction statistics"""
        stats = self.stats.copy()
        if stats['papers_processed'] > 0:
            stats['llm_call_rate'] = stats['llm_calls'] / stats['papers_processed']
        return stats


class CPUOptimizedPipeline:
    """
    Complete CPU-optimized extraction pipeline
    No GPU required!
    """
    
    # Recommended Ollama models for CPU
    RECOMMENDED_LLM_MODELS = {
        'fast': 'qwen2.5:3b',
        'balanced': 'qwen2.5:7b',
        'quality': 'qwen2.5:14b',
    }
    
    # NER model comparison
    NER_MODELS = {
        'scispacy': {
            'speed': '⚡⚡⚡ (~0.1 sec/abstract)',
            'ram': '~200MB',
            'quality': '⭐⭐',
        },
        'pubmedbert': {
            'speed': '⚡ (~3-5 sec/abstract on CPU)',
            'ram': '~2-4GB',
            'quality': '⭐⭐⭐',
        }
    }
    
    @classmethod
    def create(
        cls,
        profile: str = 'balanced',
        ner_model: str = 'scispacy',
        custom_llm_model: str = None,
        use_llm: bool = True
    ) -> HybridExtractor:
        """
        Create CPU-optimized extractor
        
        Args:
            profile: 'fast', 'balanced', or 'quality'
            ner_model: 'scispacy' (fast) or 'pubmedbert' (quality)
            custom_llm_model: Override LLM model selection
            use_llm: Enable LLM (False = NER+patterns only)
            
        Returns:
            Configured HybridExtractor
        """
        llm_model = custom_llm_model or cls.RECOMMENDED_LLM_MODELS.get(profile, 'qwen2.5:7b')
        
        # Adjust threshold based on profile
        thresholds = {
            'fast': 0.7,
            'balanced': 0.5,
            'quality': 0.3,
        }
        
        return HybridExtractor(
            ner_model=ner_model,
            llm_provider='ollama',
            llm_model=llm_model,
            use_llm=use_llm,
            llm_threshold=thresholds.get(profile, 0.5)
        )


if __name__ == "__main__":
    # Test hybrid extraction
    test_abstract = """
    Alpha-synuclein aggregation is a hallmark of Parkinson's disease. 
    We investigated the aggregation kinetics of alpha-synuclein using 
    ThT fluorescence assays and transmission electron microscopy (TEM). 
    Our results show that alpha-synuclein forms β-sheet rich fibrils 
    through a nucleation-dependent mechanism. The aggregation was 
    irreversible under physiological conditions. HSP70 chaperone can 
    prevent aggregation. We also observed that Tau showed similar behavior.
    """
    
    print("Testing Hybrid Extractor (CPU-only mode)...")
    print("=" * 60)
    
    # Test without LLM (NER + patterns only)
    extractor = HybridExtractor(use_llm=False)
    
    result = extractor.extract(
        title="Alpha-synuclein aggregation in Parkinson's disease",
        abstract=test_abstract,
        pmid="12345678"
    )
    
    print(f"\nExtraction method: {result.extraction_method}")
    print(f"Proteins found: {result.ner_proteins_count}")
    print(f"LLM called: {result.llm_called}")
    print(f"Paper relevance: {result.paper_relevance}")
    
    print("\nExtracted proteins:")
    for p in result.proteins:
        print(f"\n  {p['protein_name']}:")
        print(f"    Forms aggregates: {p.get('forms_aggregates')}")
        print(f"    Structure: {p.get('aggregate_structure')}")
        print(f"    Reversibility: {p.get('reversibility')}")
        print(f"    Class: {p.get('functional_class')}")
        print(f"    Methods: {p.get('experimental_methods')}")
        print(f"    Confidence: {p.get('extraction_confidence'):.2f}")
    
    print("\nStats:", extractor.get_stats())
