"""Database parsers."""

from .base import BaseParser
from .waltzdb import WaltzDBParser
from .crossbeta import CrossBetaDBParser
from .amyload import AmyLoadParser
from .amyloid_explorer import AmyloidExplorerParser
from .amyloid_atlas import AmyloidAtlasParser
from .amylobase import AmylobaseParser
from .amypro import AmyProParser
from .amylograph import AmyloGraphParser
from .cpad import CPADPeptideParser, CPADStructureParser

# Parser registry
PARSERS = {
    'waltzdb': WaltzDBParser,
    'crossbeta': CrossBetaDBParser,
    'amyload': AmyLoadParser,
    'amyloid_explorer': AmyloidExplorerParser,
    'amyloid_atlas': AmyloidAtlasParser,
    'amylobase': AmylobaseParser,
    'amypro': AmyProParser,
    'amylograph': AmyloGraphParser,
    'cpad_peptides': CPADPeptideParser,
    'cpad_structures': CPADStructureParser,
}

__all__ = [
    'BaseParser',
    'PARSERS',
    'WaltzDBParser',
    'CrossBetaDBParser',
    'AmyLoadParser',
    'AmyloidExplorerParser',
    'AmyloidAtlasParser',
    'AmylobaseParser',
    'AmyProParser',
    'AmyloGraphParser',
    'CPADPeptideParser',
    'CPADStructureParser',
]
