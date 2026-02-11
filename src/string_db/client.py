"""
STRING DB Integration Module
Retrieves protein-protein interactions to identify assembly/disassembly cofactors
"""

import requests
import time
from typing import List, Dict, Optional, Set
from dataclasses import dataclass, field
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# STRING API base URL
STRING_API_URL = "https://string-db.org/api"
STRING_VERSION = "12.0"

# Common species NCBI Taxonomy IDs
SPECIES_IDS = {
    "human": 9606,
    "mouse": 10090,
    "rat": 10116,
    "yeast": 4932,
    "ecoli": 511145,
    "drosophila": 7227,
    "celegans": 6239,
    "zebrafish": 7955,
}

# Interaction channels in STRING
INTERACTION_CHANNELS = [
    "neighborhood",
    "neighborhood_transferred",
    "fusion",
    "cooccurence",
    "homology",
    "coexpression",
    "coexpression_transferred",
    "experiments",
    "experiments_transferred",
    "database",
    "database_transferred",
    "textmining",
    "textmining_transferred",
]


@dataclass
class ProteinInteraction:
    """Represents an interaction between two proteins"""
    protein_a: str
    protein_b: str
    combined_score: float
    experimental_score: float = 0.0
    database_score: float = 0.0
    textmining_score: float = 0.0
    coexpression_score: float = 0.0
    
    # STRING identifiers
    string_id_a: str = ""
    string_id_b: str = ""
    
    # Annotations
    protein_a_annotation: str = ""
    protein_b_annotation: str = ""


@dataclass
class ProteinInfo:
    """STRING protein information"""
    string_id: str
    preferred_name: str
    protein_size: int = 0
    annotation: str = ""
    
    # Additional identifiers
    uniprot_ids: List[str] = field(default_factory=list)
    gene_names: List[str] = field(default_factory=list)


class StringDBClient:
    """Client for STRING database API"""
    
    def __init__(self, species: int = 9606):
        """
        Initialize STRING client
        
        Args:
            species: NCBI Taxonomy ID (default: human)
        """
        self.base_url = f"{STRING_API_URL}/json"
        self.species = species
        self.rate_limit_delay = 0.5  # Be nice to the API
        self.last_request = 0
    
    def _rate_limit(self):
        """Enforce rate limiting"""
        elapsed = time.time() - self.last_request
        if elapsed < self.rate_limit_delay:
            time.sleep(self.rate_limit_delay - elapsed)
        self.last_request = time.time()
    
    def _make_request(self, endpoint: str, params: Dict) -> Optional[List[Dict]]:
        """Make API request with error handling"""
        self._rate_limit()
        
        url = f"{self.base_url}/{endpoint}"
        params["species"] = self.species
        
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            logger.error(f"STRING API error: {e}")
            return None
    
    def resolve_protein(self, identifier: str, limit: int = 1) -> List[ProteinInfo]:
        """
        Resolve protein identifier to STRING ID
        
        Args:
            identifier: Protein name, gene name, or UniProt ID
            limit: Max results to return
            
        Returns:
            List of ProteinInfo objects
        """
        data = self._make_request("get_string_ids", {
            "identifiers": identifier,
            "limit": limit,
        })
        
        if not data:
            return []
        
        results = []
        for item in data:
            info = ProteinInfo(
                string_id=item.get("stringId", ""),
                preferred_name=item.get("preferredName", ""),
                protein_size=item.get("proteinSize", 0),
                annotation=item.get("annotation", ""),
            )
            results.append(info)
        
        return results
    
    def resolve_proteins_batch(
        self, 
        identifiers: List[str],
        limit_per_id: int = 1
    ) -> Dict[str, List[ProteinInfo]]:
        """
        Resolve multiple protein identifiers
        
        Args:
            identifiers: List of protein names/IDs
            limit_per_id: Max results per identifier
            
        Returns:
            Dict mapping identifier to ProteinInfo list
        """
        results = {}
        
        # STRING API accepts multiple identifiers separated by newline
        # But let's process in batches of 100 to be safe
        batch_size = 100
        
        for i in range(0, len(identifiers), batch_size):
            batch = identifiers[i:i + batch_size]
            
            data = self._make_request("get_string_ids", {
                "identifiers": "\r".join(batch),
                "limit": limit_per_id,
            })
            
            if data:
                for item in data:
                    query_id = item.get("queryItem", "")
                    info = ProteinInfo(
                        string_id=item.get("stringId", ""),
                        preferred_name=item.get("preferredName", ""),
                        protein_size=item.get("proteinSize", 0),
                        annotation=item.get("annotation", ""),
                    )
                    
                    if query_id not in results:
                        results[query_id] = []
                    results[query_id].append(info)
        
        return results
    
    def get_interactions(
        self,
        proteins: List[str],
        required_score: int = 700,
        network_type: str = "physical",
        add_nodes: int = 0
    ) -> List[ProteinInteraction]:
        """
        Get protein-protein interactions for given proteins
        
        Args:
            proteins: List of STRING IDs or protein names
            required_score: Minimum combined score (0-1000)
            network_type: "physical" (physical subnetwork) or "functional"
            add_nodes: Number of additional interactors to add
            
        Returns:
            List of ProteinInteraction objects
        """
        data = self._make_request("network", {
            "identifiers": "\r".join(proteins),
            "required_score": required_score,
            "network_type": network_type,
            "add_nodes": add_nodes,
        })
        
        if not data:
            return []
        
        interactions = []
        for item in data:
            interaction = ProteinInteraction(
                protein_a=item.get("preferredName_A", ""),
                protein_b=item.get("preferredName_B", ""),
                combined_score=item.get("score", 0),
                string_id_a=item.get("stringId_A", ""),
                string_id_b=item.get("stringId_B", ""),
            )
            
            # Get detailed scores if available
            if "nscore" in item:
                interaction.experimental_score = item.get("escore", 0)
                interaction.database_score = item.get("dscore", 0)
                interaction.textmining_score = item.get("tscore", 0)
                interaction.coexpression_score = item.get("ascore", 0)
            
            interactions.append(interaction)
        
        return interactions
    
    def get_interaction_partners(
        self,
        protein: str,
        required_score: int = 700,
        limit: int = 50
    ) -> List[ProteinInteraction]:
        """
        Get all interaction partners for a single protein
        
        Args:
            protein: STRING ID or protein name
            required_score: Minimum combined score
            limit: Maximum partners to return
            
        Returns:
            List of ProteinInteraction objects
        """
        data = self._make_request("interaction_partners", {
            "identifiers": protein,
            "required_score": required_score,
            "limit": limit,
        })
        
        if not data:
            return []
        
        interactions = []
        for item in data:
            interaction = ProteinInteraction(
                protein_a=item.get("preferredName_A", ""),
                protein_b=item.get("preferredName_B", ""),
                combined_score=item.get("score", 0),
                string_id_a=item.get("stringId_A", ""),
                string_id_b=item.get("stringId_B", ""),
            )
            interactions.append(interaction)
        
        return interactions
    
    def get_functional_enrichment(
        self,
        proteins: List[str],
    ) -> Dict[str, List[Dict]]:
        """
        Get functional enrichment for a set of proteins
        
        Args:
            proteins: List of STRING IDs or protein names
            
        Returns:
            Dict with GO terms, KEGG pathways, etc.
        """
        data = self._make_request("enrichment", {
            "identifiers": "\r".join(proteins),
        })
        
        if not data:
            return {}
        
        # Group by category
        enrichment = {}
        for item in data:
            category = item.get("category", "unknown")
            if category not in enrichment:
                enrichment[category] = []
            enrichment[category].append({
                "term": item.get("term", ""),
                "description": item.get("description", ""),
                "p_value": item.get("p_value", 1.0),
                "fdr": item.get("fdr", 1.0),
                "proteins": item.get("inputGenes", "").split(","),
            })
        
        return enrichment


# =============================================================================
# Helper functions for aggregation research
# =============================================================================

# Known chaperones and proteostasis factors
KNOWN_CHAPERONES = {
    "HSP70", "HSPA1A", "HSPA8", "HSC70",
    "HSP90", "HSP90AA1", "HSP90AB1",
    "HSP40", "DNAJA1", "DNAJB1",
    "HSP60", "HSPD1",
    "HSP27", "HSPB1",
    "CHIP", "STUB1",
    "BAG1", "BAG3",
    "CRYAB", "HSPB5",
    "clusterin", "CLU",
}

KNOWN_DISAGGREGASES = {
    "HSP104",  # yeast
    "HSP70", "HSPA1A",  # human (with co-chaperones)
    "DNAJB1",
    "HSP110", "HSPH1",
    "p97", "VCP",
}


def identify_potential_cofactors(
    interactions: List[ProteinInteraction],
    min_score: float = 0.7
) -> Dict[str, List[str]]:
    """
    Identify potential assembly/disassembly cofactors from interactions
    
    Args:
        interactions: List of protein interactions
        min_score: Minimum interaction score
        
    Returns:
        Dict with "chaperones", "disaggregases", "other" lists
    """
    result = {
        "chaperones": [],
        "disaggregases": [],
        "other_interactors": [],
    }
    
    seen = set()
    
    for interaction in interactions:
        if interaction.combined_score < min_score:
            continue
        
        for partner in [interaction.protein_a, interaction.protein_b]:
            if partner in seen:
                continue
            seen.add(partner)
            
            partner_upper = partner.upper()
            
            if partner_upper in KNOWN_CHAPERONES:
                result["chaperones"].append(partner)
            elif partner_upper in KNOWN_DISAGGREGASES:
                result["disaggregases"].append(partner)
            else:
                result["other_interactors"].append(partner)
    
    return result


def find_common_interactors(
    protein_list: List[str],
    client: StringDBClient,
    min_shared: int = 2
) -> Dict[str, Set[str]]:
    """
    Find proteins that interact with multiple proteins in the list
    
    Args:
        protein_list: List of protein names
        client: StringDBClient instance
        min_shared: Minimum number of proteins an interactor must connect
        
    Returns:
        Dict mapping interactor to set of connected proteins
    """
    interactor_to_proteins = {}
    
    for protein in protein_list:
        interactions = client.get_interaction_partners(protein)
        
        for interaction in interactions:
            partner = interaction.protein_b
            if partner == protein:
                partner = interaction.protein_a
            
            if partner not in interactor_to_proteins:
                interactor_to_proteins[partner] = set()
            interactor_to_proteins[partner].add(protein)
    
    # Filter to shared interactors
    shared = {
        k: v for k, v in interactor_to_proteins.items()
        if len(v) >= min_shared
    }
    
    return shared


if __name__ == "__main__":
    # Test the module
    print("Testing STRING DB Client...")
    
    client = StringDBClient(species=9606)  # Human
    
    # Test protein resolution
    print("\n1. Resolving 'alpha-synuclein'...")
    proteins = client.resolve_protein("alpha-synuclein")
    for p in proteins:
        print(f"   {p.preferred_name} ({p.string_id}): {p.annotation[:50]}...")
    
    # Test getting interactions
    if proteins:
        print(f"\n2. Getting interactions for {proteins[0].preferred_name}...")
        interactions = client.get_interaction_partners(
            proteins[0].string_id,
            required_score=900,
            limit=10
        )
        print(f"   Found {len(interactions)} high-confidence interactions:")
        for i in interactions[:5]:
            print(f"   - {i.protein_b}: score={i.combined_score}")
    
    # Test cofactor identification
    print("\n3. Identifying potential cofactors...")
    cofactors = identify_potential_cofactors(interactions)
    print(f"   Chaperones: {cofactors['chaperones']}")
    print(f"   Disaggregases: {cofactors['disaggregases']}")
    print(f"   Other: {cofactors['other_interactors'][:5]}...")
