"""
SnapGene File Format (SGFF) parser and writer
"""

from .reader import SgffReader
from .writer import SgffWriter
from .internal import SgffObject

__all__ = ["SgffReader", "SgffWriter", "SgffObject"]
__version__ = "0.1.0"
