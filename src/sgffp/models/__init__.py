"""
Data models for SnapGene file format
"""

from .base import SgffModel, SgffListModel
from .sequence import SgffSequence
from .feature import SgffFeature, SgffFeatureList, SgffSegment
from .history import (
    SgffHistory,
    SgffHistoryNode,
    SgffHistoryNodeContent,
    SgffHistoryTree,
    SgffHistoryTreeNode,
    SgffHistoryOligo,
    SgffInputSummary,
    HistoryOperation,
)
from .primer import SgffPrimer, SgffPrimerList
from .notes import SgffNotes
from .properties import SgffProperties
from .alignment import SgffAlignment, SgffAlignmentList

__all__ = [
    "SgffModel",
    "SgffListModel",
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
]
