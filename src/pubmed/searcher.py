"""
PubMed Search Module
Handles searching and fetching articles from PubMed/NCBI
"""

import requests
import xml.etree.ElementTree as ET
import time
import re
from typing import List, Dict, Optional, Generator
from dataclasses import dataclass, field
from datetime import datetime
import logging
import os

import sys
project_dir = os.path.dirname(os.path.abspath("amyloid_pipeline"))
sys.path.append(project_dir)
from config.settings import PUBMED_API, AGGREGATION_TERMS, EXCLUSION_TERMS

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class Article:
    """Represents a PubMed article"""
    pmid: str
    title: str
    abstract: str
    authors: List[str]
    journal: str
    year: int
    doi: Optional[str] = None
    keywords: List[str] = field(default_factory=list)
    mesh_terms: List[str] = field(default_factory=list)
    full_text_available: bool = False


class PubMedSearcher:
    """Handles PubMed API interactions"""
    
    def __init__(self, email: str = None, api_key: str = None):
        self.base_url = PUBMED_API["base_url"]
        self.email = email or PUBMED_API["email"]
        self.api_key = api_key or PUBMED_API.get("api_key")
        self.rate_limit = PUBMED_API["rate_limit_delay"]
        self.last_request_time = 0
        
    def _rate_limit(self):
        """Enforce rate limiting"""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.rate_limit:
            time.sleep(self.rate_limit - elapsed)
        self.last_request_time = time.time()
    
    def _build_params(self, **kwargs) -> Dict:
        """Build request parameters with authentication"""
        params = {"email": self.email}
        if self.api_key:
            params["api_key"] = self.api_key
        params.update(kwargs)
        return params
    
    def build_query(
        self, 
        include_terms: List[str] = None,
        exclude_terms: List[str] = None,
        date_from: str = None,
        date_to: str = None,
        require_abstract: bool = True
    ) -> str:
        """
        Build PubMed search query
        
        Args:
            include_terms: Terms to include (OR'd together)
            exclude_terms: Terms to exclude (NOT'd)
            date_from: Start date (YYYY/MM/DD)
            date_to: End date (YYYY/MM/DD)
            require_abstract: Only return articles with abstracts
        """
        include_terms = include_terms or AGGREGATION_TERMS
        exclude_terms = exclude_terms or EXCLUSION_TERMS
        
        # Build inclusion part
        include_parts = [f'"{term}"[Title/Abstract]' for term in include_terms]
        include_query = " OR ".join(include_parts)
        
        # Build exclusion part
        exclude_parts = [f'"{term}"[Title/Abstract]' for term in exclude_terms]
        exclude_query = " OR ".join(exclude_parts)
        
        # Combine
        query = f"({include_query}) NOT ({exclude_query})"
        
        # Add date filter
        if date_from or date_to:
            date_from = date_from or "1900/01/01"
            date_to = date_to or datetime.now().strftime("%Y/%m/%d")
            query += f' AND ("{date_from}"[Date - Publication] : "{date_to}"[Date - Publication])'
        
        # Require abstract
        if require_abstract:
            query += " AND hasabstract"
        
        return query
    
    def search(
        self, 
        query: str, 
        max_results: int = 1000,
        sort: str = "relevance"
    ) -> List[str]:
        """
        Search PubMed and return list of PMIDs
        Handles pagination for >10000 results (NCBI limit)
        """
        all_pmids = []
        batch_size = min(max_results, 10000)  # NCBI max is 10000
        retstart = 0
        
        while len(all_pmids) < max_results:
            self._rate_limit()
            
            current_batch = min(batch_size, max_results - len(all_pmids))
            
            params = self._build_params(
                db="pubmed",
                term=query,
                retmax=current_batch,
                retstart=retstart,
                retmode="json",
                sort=sort,
                usehistory="y"
            )
            
            url = f"{self.base_url}esearch.fcgi"
            
            try:
                response = requests.get(url, params=params)
                response.raise_for_status()
                data = response.json()
                
                result = data.get("esearchresult", {})
                pmids = result.get("idlist", [])
                total_count = int(result.get("count", 0))
                
                if not pmids:
                    break
                
                all_pmids.extend(pmids)
                
                if len(all_pmids) == 0:
                    logger.info(f"Found {total_count} articles total")
                
                logger.info(f"Fetched {len(all_pmids)}/{min(total_count, max_results)} PMIDs...")
                
                # Check if we've got all available
                if len(pmids) < current_batch or len(all_pmids) >= total_count:
                    break
                
                retstart += len(pmids)
                
            except requests.RequestException as e:
                logger.error(f"Search failed: {e}")
                break
        
        logger.info(f"Total PMIDs retrieved: {len(all_pmids)}")
        return all_pmids[:max_results]
    
    def fetch_articles(
        self, 
        pmids: List[str], 
        batch_size: int = 100
    ) -> Generator[Article, None, None]:
        """
        Fetch full article details for given PMIDs
        
        Args:
            pmids: List of PMIDs to fetch
            batch_size: Number of articles per request
            
        Yields:
            Article objects
        """
        for i in range(0, len(pmids), batch_size):
            batch = pmids[i:i + batch_size]
            self._rate_limit()
            
            params = self._build_params(
                db="pubmed",
                id=",".join(batch),
                retmode="xml"
            )
            
            url = f"{self.base_url}efetch.fcgi"
            
            try:
                response = requests.get(url, params=params)
                response.raise_for_status()
                
                articles = self._parse_pubmed_xml(response.text)
                for article in articles:
                    yield article
                    
            except requests.RequestException as e:
                logger.error(f"Fetch failed for batch starting at {i}: {e}")
                continue
    
    def _parse_pubmed_xml(self, xml_text: str) -> List[Article]:
        """Parse PubMed XML response into Article objects"""
        articles = []
        
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError as e:
            logger.error(f"XML parse error: {e}")
            return articles
        
        for article_elem in root.findall(".//PubmedArticle"):
            try:
                article = self._parse_article_element(article_elem)
                if article:
                    articles.append(article)
            except Exception as e:
                logger.warning(f"Failed to parse article: {e}")
                continue
        
        return articles
    
    def _parse_article_element(self, elem: ET.Element) -> Optional[Article]:
        """Parse single PubmedArticle XML element"""
        
        # PMID
        pmid_elem = elem.find(".//PMID")
        if pmid_elem is None:
            return None
        pmid = pmid_elem.text
        
        # Title
        title_elem = elem.find(".//ArticleTitle")
        title = title_elem.text if title_elem is not None else ""
        
        # Abstract
        abstract_parts = []
        for abs_text in elem.findall(".//AbstractText"):
            label = abs_text.get("Label", "")
            text = abs_text.text or ""
            if label:
                abstract_parts.append(f"{label}: {text}")
            else:
                abstract_parts.append(text)
        abstract = " ".join(abstract_parts)
        
        # Authors
        authors = []
        for author in elem.findall(".//Author"):
            last = author.findtext("LastName", "")
            first = author.findtext("ForeName", "")
            if last:
                authors.append(f"{last} {first}".strip())
        
        # Journal
        journal = elem.findtext(".//Journal/Title", "")
        
        # Year
        year_elem = elem.find(".//PubDate/Year")
        if year_elem is None:
            year_elem = elem.find(".//PubDate/MedlineDate")
        year = 0
        if year_elem is not None and year_elem.text:
            # Extract first 4 digits
            match = re.search(r"(\d{4})", year_elem.text)
            if match:
                year = int(match.group(1))
        
        # DOI
        doi = None
        for article_id in elem.findall(".//ArticleId"):
            if article_id.get("IdType") == "doi":
                doi = article_id.text
                break
        
        # Keywords
        keywords = [kw.text for kw in elem.findall(".//Keyword") if kw.text]
        
        # MeSH terms
        mesh_terms = []
        for mesh in elem.findall(".//MeshHeading/DescriptorName"):
            if mesh.text:
                mesh_terms.append(mesh.text)
        
        return Article(
            pmid=pmid,
            title=title,
            abstract=abstract,
            authors=authors,
            journal=journal,
            year=year,
            doi=doi,
            keywords=keywords,
            mesh_terms=mesh_terms
        )
    
    def search_and_fetch(
        self,
        query: str = None,
        max_results: int = 1000,
        **query_kwargs
    ) -> Generator[Article, None, None]:
        """
        Convenience method: search and fetch in one call
        
        Args:
            query: Pre-built query string (if None, builds from settings)
            max_results: Maximum articles to return
            **query_kwargs: Arguments for build_query() if query not provided
            
        Yields:
            Article objects
        """
        if query is None:
            query = self.build_query(**query_kwargs)
        
        logger.info(f"Searching with query: {query[:200]}...")
        pmids = self.search(query, max_results=max_results)
        
        logger.info(f"Fetching {len(pmids)} articles...")
        yield from self.fetch_articles(pmids)


# =============================================================================
# Specialized search functions
# =============================================================================

def search_aggregating_proteins(
    searcher: PubMedSearcher = None,
    max_results: int = 1000,
    date_from: str = "2000/01/01"
) -> Generator[Article, None, None]:
    """
    Search for papers about aggregating proteins (excluding prions)
    """
    if searcher is None:
        searcher = PubMedSearcher()
    
    yield from searcher.search_and_fetch(
        max_results=max_results,
        date_from=date_from
    )


def search_functional_amyloids(
    searcher: PubMedSearcher = None,
    max_results: int = 500
) -> Generator[Article, None, None]:
    """
    Search specifically for functional amyloid papers
    """
    if searcher is None:
        searcher = PubMedSearcher()
    
    functional_terms = [
        "functional amyloid",
        "beneficial amyloid",
        "physiological amyloid",
        "biofilm amyloid",
        "secretory granule amyloid",
        "hormone storage amyloid",
    ]
    
    yield from searcher.search_and_fetch(
        include_terms=functional_terms,
        max_results=max_results
    )


def search_pathological_amyloids(
    searcher: PubMedSearcher = None,
    max_results: int = 500
) -> Generator[Article, None, None]:
    """
    Search specifically for pathological amyloid papers
    """
    if searcher is None:
        searcher = PubMedSearcher()
    
    pathological_terms = [
        "amyloidosis",
        "amyloid disease",
        "neurodegenerative amyloid",
        "toxic amyloid",
        "pathological aggregation",
    ]
    
    yield from searcher.search_and_fetch(
        include_terms=pathological_terms,
        max_results=max_results
    )


if __name__ == "__main__":
    # Test the module
    searcher = PubMedSearcher(email="test@example.com")
    
    # Build and print query
    query = searcher.build_query(date_from="2020/01/01")
    print(f"Query: {query}\n")
    
    # Test search (limited)
    print("Testing search...")
    pmids = searcher.search(query, max_results=5)
    print(f"Found PMIDs: {pmids}\n")
    
    # Test fetch
    if pmids:
        print("Testing fetch...")
        for article in searcher.fetch_articles(pmids):
            print(f"  - {article.pmid}: {article.title[:60]}...")
            print(f"    Year: {article.year}, Authors: {len(article.authors)}")
            print(f"    Abstract length: {len(article.abstract)} chars")
            print()
