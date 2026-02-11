"""
Main Pipeline Orchestrator
Coordinates PubMed search, LLM extraction, and STRING DB enrichment
"""

import argparse
import logging
from typing import List, Dict, Optional, Generator
from dataclasses import dataclass
from datetime import datetime
import json
import os

import sys
project_dir = os.path.dirname(os.path.abspath("amyloid_pipeline"))
sys.path.append(project_dir)

from src.pubmed.searcher import PubMedSearcher, Article
from src.llm.extractor import LLMExtractor, RuleBasedFilter
from src.hybrid.extractor import HybridExtractor, CPUOptimizedPipeline
from src.string_db.client import StringDBClient, identify_potential_cofactors
from src.database.storage import AmyloidDatabase
from src.ner.normalizer import normalize_and_enrich, ProteinNormalizer
from config.settings import PUBMED_API, LLM_CONFIG, STRING_DB_API

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class PipelineConfig:
    """Pipeline configuration"""
    # PubMed settings
    pubmed_email: str = "your_email@example.com"
    pubmed_api_key: Optional[str] = None
    max_papers: int = 1000
    date_from: str = "2000/01/01"
    
    # Extraction mode
    extraction_mode: str = "hybrid"  # "hybrid", "llm_only", "ner_only"
    
    # NER settings
    ner_model: str = "scispacy"  # "scispacy" (fast) or "pubmedbert" (quality)
    
    # LLM settings (for hybrid and llm_only modes)
    llm_provider: str = "ollama"  # anthropic, openai, ollama
    llm_model: Optional[str] = None  # None = use provider default
    llm_api_key: Optional[str] = None
    ollama_url: str = "http://localhost:11434"
    
    # Hybrid mode settings
    hybrid_profile: str = "balanced"  # fast, balanced, quality
    llm_threshold: float = 0.5  # Call LLM if pattern confidence < threshold
    
    # STRING DB settings
    species: int = 9606  # human
    string_score_threshold: int = 700
    
    # Processing settings
    use_prefilter: bool = True
    skip_processed: bool = True
    batch_size: int = 10
    
    # Database
    db_path: str = "data/amyloid_aggregation.db"


class AmyloidPipeline:
    """
    Main pipeline for building amyloid/aggregation protein database
    
    Steps:
    1. Search PubMed for relevant papers
    2. Pre-filter using rules (optional)
    3. Extract protein information using LLM
    4. Store in database
    5. Enrich with STRING DB interactions
    6. Secondary LLM analysis for partner proteins
    """
    
    def __init__(self, config: PipelineConfig = None):
        self.config = config or PipelineConfig()
        
        # Initialize components
        self.pubmed = PubMedSearcher(
            email=self.config.pubmed_email,
            api_key=self.config.pubmed_api_key
        )
        
        # Initialize extractor based on mode
        if self.config.extraction_mode == "hybrid":
            self.extractor = HybridExtractor(
                ner_model=self.config.ner_model,
                llm_provider=self.config.llm_provider,
                llm_model=self.config.llm_model,
                llm_api_key=self.config.llm_api_key,
                ollama_url=self.config.ollama_url,
                use_llm=True,
                llm_threshold=self.config.llm_threshold
            )
            self.use_hybrid = True
        elif self.config.extraction_mode == "ner_only":
            self.extractor = HybridExtractor(
                ner_model=self.config.ner_model,
                use_llm=False
            )
            self.use_hybrid = True
        else:  # llm_only
            self.extractor = LLMExtractor(
                provider=self.config.llm_provider,
                model=self.config.llm_model,
                api_key=self.config.llm_api_key,
                ollama_base_url=self.config.ollama_url
            )
            self.use_hybrid = False
        
        self.string_client = StringDBClient(species=self.config.species)
        
        self.db = AmyloidDatabase(self.config.db_path)
        
        # Statistics
        self.stats = {
            'papers_searched': 0,
            'papers_filtered': 0,
            'papers_processed': 0,
            'proteins_extracted': 0,
            'cofactors_found': 0,
            'string_enriched': 0,
            'extraction_mode': self.config.extraction_mode,
        }
    
    def run_full_pipeline(
        self,
        custom_query: str = None,
        dry_run: bool = False
    ) -> Dict:
        """
        Run the complete pipeline
        
        Args:
            custom_query: Optional custom PubMed query
            dry_run: If True, don't make LLM calls or write to DB
            
        Returns:
            Pipeline statistics
        """
        logger.info("="*60)
        logger.info("Starting Amyloid Aggregation Pipeline")
        logger.info("="*60)
        
        # Step 1: Search PubMed
        logger.info("\n[STEP 1] Searching PubMed...")
        articles = list(self._search_pubmed(custom_query))
        self.stats['papers_searched'] = len(articles)
        logger.info(f"Found {len(articles)} articles")
        
        if dry_run:
            logger.info("DRY RUN - stopping before LLM extraction")
            return self.stats
        
        # Step 2: Pre-filter (optional)
        if self.config.use_prefilter:
            logger.info("\n[STEP 2] Pre-filtering articles...")
            articles = list(self._prefilter_articles(articles))
            self.stats['papers_filtered'] = len(articles)
            logger.info(f"After filtering: {len(articles)} articles")
        
        # Step 3: LLM extraction
        logger.info("\n[STEP 3] Extracting protein information...")
        all_proteins = []
        skipped_count = 0
        
        for i, article in enumerate(articles):
            if self.config.skip_processed and self.db.is_paper_processed(article.pmid):
                skipped_count += 1
                continue
            
            if (i + 1) % 100 == 0 or i == 0:
                logger.info(f"Processing {i+1}/{len(articles)}: PMID {article.pmid}")
            
            # Extract based on mode
            if self.use_hybrid:
                result = self.extractor.extract(
                    title=article.title,
                    abstract=article.abstract,
                    pmid=article.pmid
                )
                proteins_data = result.proteins
                paper_relevance = result.paper_relevance
            else:
                result = self.extractor.extract_from_abstract(
                    title=article.title,
                    abstract=article.abstract,
                    pmid=article.pmid
                )
                proteins_data = result.get('proteins', [])
                paper_relevance = result.get('paper_relevance', 'unknown')
            
            # Store paper
            self.db.add_paper({
                'pmid': article.pmid,
                'title': article.title,
                'abstract': article.abstract,
                'authors': article.authors,
                'journal': article.journal,
                'year': article.year,
                'doi': article.doi,
                'relevance': paper_relevance,
            })
            
            # Normalize and enrich proteins before storing
            if proteins_data:
                proteins_data = normalize_and_enrich(proteins_data, use_uniprot=True)
            
            # Store proteins
            for protein_data in proteins_data:
                protein_id = self.db.add_protein(protein_data)
                self.db.link_protein_paper(protein_id, article.pmid)
                
                # Store methods
                if protein_data.get('experimental_methods'):
                    self.db.add_experimental_methods(
                        protein_id, 
                        article.pmid,
                        protein_data['experimental_methods']
                    )
                
                # Store cofactors from literature
                for cofactor in protein_data.get('assembly_cofactors', []):
                    self.db.add_cofactor(
                        protein_id, cofactor, 'assembly', 
                        'literature', 0.7, article.pmid
                    )
                    self.stats['cofactors_found'] += 1
                
                for cofactor in protein_data.get('disassembly_cofactors', []):
                    self.db.add_cofactor(
                        protein_id, cofactor, 'disassembly',
                        'literature', 0.7, article.pmid
                    )
                    self.stats['cofactors_found'] += 1
                
                all_proteins.append({
                    'id': protein_id,
                    **protein_data
                })
            
            self.stats['papers_processed'] += 1
            self.stats['proteins_extracted'] += len(proteins_data)
        
        if skipped_count > 0:
            logger.info(f"Skipped {skipped_count} already processed papers (use --force to reprocess)")
        
        # Step 4: STRING DB enrichment
        logger.info("\n[STEP 4] Enriching with STRING DB interactions...")
        self._enrich_with_string(all_proteins)
        
        # Final statistics
        logger.info("\n" + "="*60)
        logger.info("Pipeline Complete!")
        logger.info("="*60)
        self._print_statistics()
        
        return self.stats
    
    def _search_pubmed(self, custom_query: str = None) -> Generator[Article, None, None]:
        """Search PubMed for relevant articles"""
        if custom_query:
            query = custom_query
        else:
            query = self.pubmed.build_query(date_from=self.config.date_from)
        
        logger.debug(f"Query: {query[:200]}...")
        
        yield from self.pubmed.search_and_fetch(
            query=query,
            max_results=self.config.max_papers
        )
    
    def _prefilter_articles(
        self, 
        articles: List[Article]
    ) -> Generator[Article, None, None]:
        """Pre-filter articles using rule-based filtering"""
        for article in articles:
            should_process, confidence = RuleBasedFilter.should_process(
                article.title,
                article.abstract
            )
            if should_process:
                yield article
            else:
                logger.debug(f"Filtered out: {article.pmid}")
    
    def _enrich_with_string(self, proteins: List[Dict]):
        """Enrich proteins with STRING DB interactions"""
        for protein in proteins:
            protein_name = protein.get('protein_name')
            protein_id = protein.get('id')
            
            if not protein_name or not protein_id:
                continue
            
            # Resolve to STRING ID
            resolved = self.string_client.resolve_protein(protein_name)
            if not resolved:
                logger.debug(f"Could not resolve in STRING: {protein_name}")
                continue
            
            string_id = resolved[0].string_id
            
            # Get interactions
            interactions = self.string_client.get_interaction_partners(
                string_id,
                required_score=self.config.string_score_threshold,
                limit=30
            )
            
            if not interactions:
                continue
            
            # Store interactions
            interaction_dicts = [{
                'partner_name': i.protein_b if i.protein_a == resolved[0].preferred_name else i.protein_a,
                'partner_string_id': i.string_id_b if i.protein_a == resolved[0].preferred_name else i.string_id_a,
                'combined_score': i.combined_score,
                'experimental_score': i.experimental_score,
                'database_score': i.database_score,
                'textmining_score': i.textmining_score,
            } for i in interactions]
            
            self.db.add_string_interactions(protein_id, interaction_dicts)
            
            # Identify potential cofactors
            cofactors = identify_potential_cofactors(interactions)
            
            for chaperone in cofactors['chaperones']:
                self.db.add_cofactor(
                    protein_id, chaperone, 'disassembly',
                    'string_db', 0.8
                )
                self.stats['cofactors_found'] += 1
            
            for disaggregase in cofactors['disaggregases']:
                self.db.add_cofactor(
                    protein_id, disaggregase, 'disassembly',
                    'string_db', 0.9
                )
                self.stats['cofactors_found'] += 1
            
            self.stats['string_enriched'] += 1
            logger.debug(f"STRING enriched: {protein_name} ({len(interactions)} interactions)")
    
    def _print_statistics(self):
        """Print pipeline statistics"""
        logger.info(f"\nPipeline Statistics:")
        logger.info(f"  Extraction mode: {self.stats.get('extraction_mode', 'unknown')}")
        logger.info(f"  Papers searched: {self.stats['papers_searched']}")
        logger.info(f"  Papers after filter: {self.stats['papers_filtered']}")
        logger.info(f"  Papers processed: {self.stats['papers_processed']}")
        logger.info(f"  Proteins extracted: {self.stats['proteins_extracted']}")
        logger.info(f"  Cofactors found: {self.stats['cofactors_found']}")
        logger.info(f"  STRING enriched: {self.stats['string_enriched']}")
        
        # Hybrid-specific stats
        if self.use_hybrid and hasattr(self.extractor, 'get_stats'):
            hybrid_stats = self.extractor.get_stats()
            logger.info(f"\nHybrid Extraction Stats:")
            logger.info(f"  NER extractions: {hybrid_stats.get('ner_extractions', 0)}")
            logger.info(f"  LLM calls: {hybrid_stats.get('llm_calls', 0)}")
            logger.info(f"  LLM skipped: {hybrid_stats.get('llm_skipped', 0)}")
            if 'llm_call_rate' in hybrid_stats:
                logger.info(f"  LLM call rate: {hybrid_stats['llm_call_rate']:.1%}")
        
        # Database stats
        db_stats = self.db.get_statistics()
        logger.info(f"\nDatabase totals:")
        logger.info(f"  Total proteins: {db_stats['total_proteins']}")
        logger.info(f"  Total papers: {db_stats['total_papers']}")
        logger.info(f"  By class: {db_stats['by_class']}")
    
    # ==========================================================================
    # Specialized pipeline methods
    # ==========================================================================
    
    def run_cofactor_analysis(self, min_shared: int = 2) -> Dict[str, List[str]]:
        """
        Analyze cofactors across functional vs pathological proteins
        
        Returns:
            Dict with shared and unique cofactors
        """
        logger.info("Analyzing cofactors...")
        
        # Get proteins by class
        functional = self.db.get_all_proteins(functional_class='functional')
        pathological = self.db.get_all_proteins(functional_class='pathological')
        
        # Get cofactors
        func_cofactors = set()
        path_cofactors = set()
        
        for p in functional:
            cofactors = self.db.get_cofactors(protein_id=p['id'])
            func_cofactors.update(c['cofactor_name'] for c in cofactors)
        
        for p in pathological:
            cofactors = self.db.get_cofactors(protein_id=p['id'])
            path_cofactors.update(c['cofactor_name'] for c in cofactors)
        
        return {
            'shared': list(func_cofactors & path_cofactors),
            'functional_only': list(func_cofactors - path_cofactors),
            'pathological_only': list(path_cofactors - func_cofactors),
        }
    
    def analyze_partner_proteins(self, protein_ids: List[int] = None):
        """
        Secondary LLM analysis: check if interaction partners are also aggregating
        """
        logger.info("Analyzing partner proteins...")
        
        # Get all unique partners from STRING
        with self.db._connection() as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT DISTINCT partner_name FROM string_interactions
                WHERE combined_score >= ?
            ''', (self.config.string_score_threshold / 1000,))
            partners = [row[0] for row in cursor.fetchall()]
        
        logger.info(f"Found {len(partners)} unique interaction partners")
        
        # For each partner, search PubMed to check if it's also aggregating
        # This is a simplified version - in practice you'd want more sophisticated logic
        for partner in partners[:10]:  # Limit for testing
            query = f'"{partner}"[Title/Abstract] AND (amyloid OR aggregat*)'
            pmids = self.pubmed.search(query, max_results=5)
            
            if pmids:
                logger.info(f"  {partner}: {len(pmids)} papers found - potential aggregator")


def main():
    """CLI entry point"""
    parser = argparse.ArgumentParser(
        description='Amyloid/Aggregation Protein Database Pipeline'
    )
    
    parser.add_argument('--email', required=True, help='PubMed email (required)')
    parser.add_argument('--pubmed-key', help='PubMed API key (optional)')
    
    # Extraction mode
    parser.add_argument('--mode', default='hybrid', 
                        choices=['hybrid', 'llm_only', 'ner_only'],
                        help='Extraction mode (hybrid recommended for CPU)')
    parser.add_argument('--profile', default='balanced',
                        choices=['fast', 'balanced', 'quality'],
                        help='Hybrid profile (fast=less LLM, quality=more LLM)')
    parser.add_argument('--ner-model', default='scispacy',
                        choices=['scispacy', 'pubmedbert'],
                        help='NER model: scispacy (fast) or pubmedbert (quality)')
    
    # LLM settings
    parser.add_argument('--llm-provider', default='ollama', 
                        choices=['anthropic', 'openai', 'ollama'],
                        help='LLM provider (ollama = free local)')
    parser.add_argument('--llm-key', help='LLM API key (not needed for ollama)')
    parser.add_argument('--llm-model', help='Model name (e.g., qwen2.5:7b for ollama)')
    parser.add_argument('--ollama-url', default='http://localhost:11434',
                        help='Ollama server URL')
    parser.add_argument('--max-papers', type=int, default=100)
    parser.add_argument('--date-from', default='2020/01/01')
    parser.add_argument('--db-path', default='data/amyloid_aggregation.db')
    parser.add_argument('--dry-run', action='store_true', help='Search only, no extraction')
    parser.add_argument('--force', action='store_true', help='Reprocess already processed papers')
    parser.add_argument('--no-prefilter', action='store_true')
    parser.add_argument('--query', help='Custom PubMed query')
    parser.add_argument('-v', '--verbose', action='store_true')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Create config
    config = PipelineConfig(
        pubmed_email=args.email,
        pubmed_api_key=args.pubmed_key,
        extraction_mode=args.mode,
        ner_model=args.ner_model,
        hybrid_profile=args.profile,
        llm_provider=args.llm_provider,
        llm_model=args.llm_model,
        llm_api_key=args.llm_key,
        ollama_url=args.ollama_url,
        max_papers=args.max_papers,
        date_from=args.date_from,
        db_path=args.db_path,
        use_prefilter=not args.no_prefilter,
        skip_processed=not args.force,
    )
    
    # Ensure data directory exists
    os.makedirs(os.path.dirname(config.db_path), exist_ok=True)
    
    # Run pipeline
    pipeline = AmyloidPipeline(config)
    stats = pipeline.run_full_pipeline(
        custom_query=args.query,
        dry_run=args.dry_run
    )
    
    # Save stats
    with open('data/pipeline_stats.json', 'w') as f:
        json.dump(stats, f, indent=2)


if __name__ == "__main__":
    main()
