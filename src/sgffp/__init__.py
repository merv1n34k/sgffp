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
    SgffHistoryNodeContent,
    SgffHistoryTree,
    SgffHistoryTreeNode,
    SgffHistoryOligo,
    SgffInputSummary,
    HistoryOperation,
    SgffPrimer,
    SgffPrimerList,
    SgffNotes,
    SgffProperties,
    SgffAlignment,
    SgffAlignmentList,
    SgffTrace,
    SgffTraceList,
    SgffTraceClip,
    SgffTraceSamples,
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
    "SgffHistoryNodeContent",
    "SgffHistoryTree",
    "SgffHistoryTreeNode",
    "SgffHistoryOligo",
    "SgffInputSummary",
    "HistoryOperation",
    "SgffPrimer",
    "SgffPrimerList",
    "SgffNotes",
    "SgffProperties",
    "SgffAlignment",
    "SgffAlignmentList",
    "SgffTrace",
    "SgffTraceList",
    "SgffTraceClip",
    "SgffTraceSamples",
]
