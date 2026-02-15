"""Dataset unification and deduplication."""

import logging
from collections import defaultdict
from typing import Dict, List, Tuple

from .models import AmyloidEntry


class DatasetUnifier:
    """Unifies entries from multiple databases with deduplication and conflict resolution."""
    
    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        self.all_entries: List[AmyloidEntry] = []
        self.dedup_map: Dict[Tuple, List[AmyloidEntry]] = defaultdict(list)
    
    def add_entries(self, entries: List[AmyloidEntry], db_name: str):
        self.logger.info(f"Adding {len(entries)} entries from {db_name}")
        self.all_entries.extend(entries)
        for entry in entries:
            self.dedup_map[entry.get_dedup_key()].append(entry)
    
    def get_stats(self) -> Dict[str, int]:
        return {
            "total_entries": len(self.all_entries),
            "unique_keys": len(self.dedup_map),
            "duplicates": len(self.all_entries) - len(self.dedup_map)
        }
    
    def detect_conflicts(self) -> List[Dict]:
        """Find entries with conflicting amyloid labels."""
        conflicts = []
        for key, entries in self.dedup_map.items():
            if len(entries) > 1:
                labels = set(e.experimental_label for e in entries if e.experimental_label)
                if len(labels) > 1:
                    conflicts.append({
                        "key": key,
                        "entries": entries,
                        "labels": labels,
                        "sources": [e.source_db for e in entries]
                    })
        self.logger.info(f"Detected {len(conflicts)} conflicts")
        return conflicts
    
    def build_consensus(self) -> List[AmyloidEntry]:
        """Build consensus using evidence-weighted voting."""
        consensus = []
        for key, entries in self.dedup_map.items():
            if len(entries) == 1:
                consensus.append(entries[0])
            else:
                labels = set(e.experimental_label for e in entries if e.experimental_label)
                if len(labels) <= 1:
                    best = max(entries, key=lambda e: e.evidence_weight)
                    consensus.append(best)
                else:
                    label_weights = defaultdict(float)
                    label_entries = defaultdict(list)
                    for entry in entries:
                        if entry.experimental_label:
                            label_weights[entry.experimental_label] += entry.evidence_weight
                            label_entries[entry.experimental_label].append(entry)
                    
                    sorted_labels = sorted(label_weights.items(), key=lambda x: x[1], reverse=True)
                    if len(sorted_labels) >= 2 and sorted_labels[0][1] > sorted_labels[1][1] * 1.5:
                        best = max(label_entries[sorted_labels[0][0]], key=lambda e: e.evidence_weight)
                        consensus.append(best)
        
        self.logger.info(f"Built consensus with {len(consensus)} entries")
        return consensus
