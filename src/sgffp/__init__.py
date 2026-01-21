"""
SnapGene File Format (SGFF) parser and writer
"""

from .reader import SgffReader
from .writer import SgffWriter
from .internal import SgffObject, Cookie, BlockList

__all__ = ["SgffReader", "SgffWriter", "SgffObject", "Cookie", "BlockList"]
__version__ = "0.6.0"

