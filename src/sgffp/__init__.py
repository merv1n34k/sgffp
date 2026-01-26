"""
SnapGene File Format (SGFF) parser and writer
"""

from .reader import SgffReader
from .writer import SgffWriter
from .internal import SgffObject, Cookie, BlockList
from .models import (
    SgffSequence,
    SgffFeature,
    SgffFeatureList,
    SgffSegment,
    SgffHistory,
    SgffHistoryNode,
    SgffPrimer,
    SgffPrimerList,
    SgffNotes,
    SgffProperties,
    SgffAlignment,
    SgffAlignmentList,
)

__all__ = [
    "SgffReader",
    "SgffWriter",
    "SgffObject",
    "Cookie",
    "BlockList",
    "SgffSequence",
    "SgffFeature",
    "SgffFeatureList",
    "SgffSegment",
    "SgffHistory",
    "SgffHistoryNode",
    "SgffPrimer",
    "SgffPrimerList",
    "SgffNotes",
    "SgffProperties",
    "SgffAlignment",
    "SgffAlignmentList",
]

