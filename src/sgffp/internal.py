"""
Internal data structures for SGFF representation
"""

from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Iterator

from .models import (
    SgffSequence,
    SgffFeature,
    SgffFeatureList,
    SgffSegment,
    SgffHistory,
    SgffPrimer,
    SgffPrimerList,
    SgffNotes,
    SgffProperties,
    SgffAlignmentList,
    SgffTraceList,
)


@dataclass
class Cookie:
    """File header metadata"""

    type_of_sequence: int = 1
    export_version: int = 15
    import_version: int = 19


class BlockList:
    """Wrapper for accessing blocks of a single type"""

    def __init__(self, block_type: int, items: List[Any]):
        self._type = block_type
        self._items = items

    @property
    def type(self) -> int:
        return self._type

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self) -> Iterator[Any]:
        return iter(self._items)

    def __getitem__(self, idx: int) -> Any:
        return self._items[idx]

    def get(self, idx: int = 0) -> Any:
        """Get item at index, default first"""
        if not self._items:
            return None
        return self._items[idx]

    @property
    def first(self) -> Any:
        """First item"""
        return self.get(0)

    @property
    def last(self) -> Any:
        """Last item"""
        return self.get(-1) if self._items else None


@dataclass
class SgffObject:
    """Container for SnapGene file data."""

    cookie: Cookie
    blocks: Dict[int, List[Any]] = field(default_factory=dict)

    TYPE_MAP = {"dna": 1, "protein": 2, "rna": 7}
    BLOCK_MAP = {"dna": 0, "protein": 21, "rna": 32}

    @classmethod
    def new(
        cls,
        sequence: str = "",
        *,
        topology: str = "linear",
        strandedness: str = "double",
        sequence_type: str = "dna",
    ) -> "SgffObject":
        """Create a new SnapGene file with sensible defaults."""
        if sequence_type not in cls.TYPE_MAP:
            raise ValueError(
                f"Invalid sequence_type: {sequence_type!r} "
                f"(expected one of {list(cls.TYPE_MAP)})"
            )

        sgff = cls(cookie=Cookie(type_of_sequence=cls.TYPE_MAP[sequence_type]))
        if sequence:
            sgff.sequence.data.value = sequence
            sgff.sequence.data.topology = topology
            sgff.sequence.data.strandedness = strandedness
            sgff.sequence._block_id = cls.BLOCK_MAP[sequence_type]
            sgff.sequence._sync()
        return sgff

    def __post_init__(self):
        self._sequence: Optional[SgffSequence] = None
        self._features: Optional[SgffFeatureList] = None
        self._history: Optional[SgffHistory] = None
        self._primers: Optional[SgffPrimerList] = None
        self._notes: Optional[SgffNotes] = None
        self._properties: Optional[SgffProperties] = None
        self._alignments: Optional[SgffAlignmentList] = None
        self._traces: Optional[SgffTraceList] = None

    def invalidate(self) -> None:
        """Clear cached model instances."""
        self._sequence = None
        self._features = None
        self._history = None
        self._primers = None
        self._notes = None
        self._properties = None
        self._alignments = None
        self._traces = None

    @property
    def types(self) -> List[int]:
        """List of all block type IDs present"""
        return list(self.blocks.keys())

    def type(self, block_id: int) -> BlockList:
        """Get blocks by type ID."""
        items = self.blocks.get(block_id, [])
        return BlockList(block_id, items)

    def set(self, block_id: int, value: Any) -> None:
        """Add value to block type."""
        if block_id not in self.blocks:
            self.blocks[block_id] = []
        self.blocks[block_id].append(value)

    def remove(self, block_id: int, idx: int = 0) -> bool:
        """Remove single item from block type"""
        if block_id not in self.blocks:
            return False
        items = self.blocks[block_id]
        if idx >= len(items):
            return False
        items.pop(idx)
        if not items:
            del self.blocks[block_id]
        return True

    def bset(self, block_id: int, value: Any) -> None:
        """Set entire block."""
        if not isinstance(value, list):
            value = [value]
        self.blocks[block_id] = value

    def bremove(self, block_id: int) -> bool:
        """Remove entire block"""
        if block_id in self.blocks:
            del self.blocks[block_id]
            return True
        return False

    def block(self, block_id: int) -> Optional[Any]:
        """Get first block value by ID."""
        items = self.blocks.get(block_id, [])
        return items[0] if items else None

    @property
    def sequence(self) -> SgffSequence:
        """DNA/RNA/Protein sequence."""
        if self._sequence is None:
            self._sequence = SgffSequence(self.blocks)
        return self._sequence

    @property
    def features(self) -> SgffFeatureList:
        """Annotation features."""
        if self._features is None:
            self._features = SgffFeatureList(self.blocks)
        return self._features

    @property
    def history(self) -> SgffHistory:
        """Edit history."""
        if self._history is None:
            self._history = SgffHistory(self.blocks)
        return self._history

    @property
    def primers(self) -> SgffPrimerList:
        """Primers."""
        if self._primers is None:
            self._primers = SgffPrimerList(self.blocks)
        return self._primers

    @property
    def notes(self) -> SgffNotes:
        """File notes."""
        if self._notes is None:
            self._notes = SgffNotes(self.blocks)
        return self._notes

    @property
    def properties(self) -> SgffProperties:
        """Sequence properties."""
        if self._properties is None:
            self._properties = SgffProperties(self.blocks)
        return self._properties

    @property
    def alignments(self) -> SgffAlignmentList:
        """Alignable sequences."""
        if self._alignments is None:
            self._alignments = SgffAlignmentList(self.blocks)
        return self._alignments

    @property
    def traces(self) -> SgffTraceList:
        """Sequence traces / chromatograms."""
        if self._traces is None:
            self._traces = SgffTraceList(self.blocks)
        return self._traces

    @property
    def has_notes(self) -> bool:
        return 6 in self.blocks

    @property
    def has_properties(self) -> bool:
        return 8 in self.blocks

    @property
    def has_history(self) -> bool:
        return any(bid in self.blocks for bid in (7, 11, 29, 30))

    @property
    def has_features(self) -> bool:
        return 10 in self.blocks

    @property
    def has_primers(self) -> bool:
        return 5 in self.blocks

    @property
    def has_alignments(self) -> bool:
        return 17 in self.blocks

    @property
    def has_traces(self) -> bool:
        return 16 in self.blocks

    def add_feature(
        self,
        name: str,
        type: str,
        start: int,
        end: int,
        strand: str = "+",
        **kwargs,
    ) -> "SgffObject":
        """Add an annotation feature. Returns self for chaining."""
        feature = SgffFeature(
            name=name,
            type=type,
            strand=strand,
            segments=[SgffSegment(start=start, end=end)],
            **kwargs,
        )
        self.features.add(feature)
        return self

    def add_primer(
        self,
        name: str,
        sequence: str,
        bind_position: Optional[int] = None,
        bind_strand: str = "+",
    ) -> "SgffObject":
        """Add a primer. Returns self for chaining."""
        primer = SgffPrimer(
            name=name,
            sequence=sequence,
            bind_position=bind_position,
            bind_strand=bind_strand,
        )
        self.primers.add(primer)
        return self
