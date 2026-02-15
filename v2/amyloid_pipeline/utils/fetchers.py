"""UniProt and NCBI sequence fetching utilities."""

import json
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from urllib.request import urlopen, Request

from ..config import UNIPROT_API, REQUEST_DELAY, REQUEST_TIMEOUT


class SequenceCache:
    """Manages sequence caching to/from JSON file."""
    
    def __init__(self, filepath: str = None):
        self.filepath = filepath
        self.uniprot: Dict[str, Optional[str]] = {}
        self.protein_info: Dict[str, Dict] = {}
        self.failed: Set[str] = set()
        if filepath and Path(filepath).exists():
            self.load()
    
    def load(self):
        try:
            with open(self.filepath, 'r') as f:
                data = json.load(f)
                self.uniprot = data.get('uniprot', {})
                self.protein_info = data.get('protein_info', {})
                self.failed = set(data.get('failed', []))
        except Exception:
            pass
    
    def save(self):
        if not self.filepath:
            return
        try:
            data = {
                'uniprot': {k: v for k, v in self.uniprot.items() if v},
                'protein_info': self.protein_info,
                'failed': list(self.failed)
            }
            with open(self.filepath, 'w') as f:
                json.dump(data, f, indent=2)
        except Exception:
            pass


def fetch_uniprot_info(uniprot_id: str) -> Optional[Dict]:
    """Fetch protein info from UniProt including family."""
    if not uniprot_id or not re.match(r'^[A-Z][A-Z0-9]{5,9}$', uniprot_id.upper()):
        return None
    
    try:
        url = f"{UNIPROT_API}/{uniprot_id}.json"
        req = Request(url, headers={'User-Agent': 'Python/AmyloidDB'})
        with urlopen(req, timeout=REQUEST_TIMEOUT) as response:
            data = json.loads(response.read().decode('utf-8'))
            
            info = {
                'sequence': data.get('sequence', {}).get('value', ''),
                'organism': data.get('organism', {}).get('scientificName', ''),
                'protein_name': '',
                'protein_family': '',
            }
            
            if 'proteinDescription' in data:
                rec_name = data['proteinDescription'].get('recommendedName', {})
                if 'fullName' in rec_name:
                    info['protein_name'] = rec_name['fullName'].get('value', '')
            
            for comment in data.get('comments', []):
                if comment.get('commentType') == 'SIMILARITY':
                    texts = comment.get('texts', [])
                    if texts:
                        info['protein_family'] = texts[0].get('value', '')
                        break
            
            for xref in data.get('uniProtKBCrossReferences', []):
                if xref.get('database') in ['InterPro', 'Pfam']:
                    props = xref.get('properties', [])
                    for prop in props:
                        if prop.get('key') == 'EntryName' and not info['protein_family']:
                            info['protein_family'] = prop.get('value', '')
            
            return info
    except Exception:
        pass
    return None


def fetch_sequences_and_info(entries: List, cache: SequenceCache, logger=None):
    """Fetch sequences and protein family info for entries."""
    import logging
    if logger is None:
        logger = logging.getLogger(__name__)
    
    to_fetch = set()
    for entry in entries:
        if entry.uniprot_id and entry.uniprot_id not in cache.uniprot:
            to_fetch.add(entry.uniprot_id)
    
    if not to_fetch:
        return
    
    logger.info(f"Fetching info for {len(to_fetch)} UniProt IDs...")
    
    for i, uid in enumerate(to_fetch):
        if (i + 1) % 50 == 0:
            logger.info(f"  Progress: {i+1}/{len(to_fetch)}")
        
        info = fetch_uniprot_info(uid)
        if info:
            cache.uniprot[uid] = info.get('sequence')
            cache.protein_info[uid] = info
        else:
            cache.failed.add(uid)
        
        time.sleep(REQUEST_DELAY)
    
    cache.save()
    
    # Update entries
    for entry in entries:
        if entry.uniprot_id and entry.uniprot_id in cache.protein_info:
            info = cache.protein_info[entry.uniprot_id]
            
            if not entry.organism and info.get('organism'):
                entry.organism = info['organism']
            
            if not entry.protein_family and info.get('protein_family'):
                entry.protein_family = info['protein_family']
            
            if not entry.sequence and info.get('sequence'):
                entry.full_protein_sequence = info['sequence']
                if entry.region_start and entry.region_end:
                    seq = info['sequence']
                    start_idx = entry.region_start - 1
                    end_idx = entry.region_end
                    if 0 <= start_idx < len(seq) and start_idx < end_idx <= len(seq):
                        entry.sequence = seq[start_idx:end_idx]
