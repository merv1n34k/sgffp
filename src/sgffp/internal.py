"""
Internal data structures for SGFF representation
"""

from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Iterator, TYPE_CHECKING

if TYPE_CHECKING:
    from .models import (
        SgffSequence,
        SgffFeatureList,
        SgffHistory,
        SgffPrimerList,
        SgffNotes,
        SgffProperties,
        SgffAlignmentList,
    )


@dataclass
class Cookie:
    """File header metadata"""

    type_of_sequence: int
    export_version: int
    import_version: int


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
    """
    Container for SnapGene file data

    Blocks stored as {int: [values]} for clean access
    """

    cookie: Cookie
    blocks: Dict[int, List[Any]] = field(default_factory=dict)

    @property
    def types(self) -> List[int]:
        """List of all block type IDs present"""
        return list(self.blocks.keys())

    def type(self, block_id: int) -> BlockList:
        """Get blocks by type ID, returns BlockList for chaining"""
        items = self.blocks.get(block_id, [])
        return BlockList(block_id, items)

    def set(self, block_id: int, value: Any) -> None:
        """Add value to block type (creates or appends)"""
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
        """Set entire block (replaces if exists)"""
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
        """Raw block access by ID"""
        items = self.blocks.get(block_id, [])
        return items[0] if items else None

    @property
    def sequence(self) -> "SgffSequence":
        """DNA/RNA/Protein sequence (blocks 0, 1, 21, 32)"""
        from .models import SgffSequence
        return SgffSequence(self.blocks)

    @property
    def features(self) -> "SgffFeatureList":
        """Annotation features (block 10)"""
        from .models import SgffFeatureList
        return SgffFeatureList(self.blocks)

    @property
    def history(self) -> "SgffHistory":
        """Edit history (blocks 7, 11, 29, 30)"""
        from .models import SgffHistory
        return SgffHistory(self.blocks)

    @property
    def primers(self) -> "SgffPrimerList":
        """Primers (block 5)"""
        from .models import SgffPrimerList
        return SgffPrimerList(self.blocks)

    @property
    def notes(self) -> "SgffNotes":
        """File notes (block 6)"""
        from .models import SgffNotes
        return SgffNotes(self.blocks)

    @property
    def properties(self) -> "SgffProperties":
        """Sequence properties (block 8)"""
        from .models import SgffProperties
        return SgffProperties(self.blocks)

    @property
    def alignments(self) -> "SgffAlignmentList":
        """Alignable sequences (block 17)"""
        from .models import SgffAlignmentList
        return SgffAlignmentList(self.blocks)

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

