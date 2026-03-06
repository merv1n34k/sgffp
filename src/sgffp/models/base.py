"""
Base classes for SGFF data models
"""

from typing import Dict, List, Any, Optional, Iterator, TypeVar, Generic

T = TypeVar("T")


class SgffModel:
    """Base wrapper for block data."""

    BLOCK_IDS: tuple = ()

    def __init__(self, blocks: Dict[int, List[Any]]):
        self._blocks = blocks

    def _get_block(self, block_id: int) -> Optional[Any]:
        """Get first item from block"""
        items = self._blocks.get(block_id, [])
        return items[0] if items else None

    def _set_block(self, block_id: int, value: Any) -> None:
        """Set block value."""
        self._blocks[block_id] = [value]

    def _get_blocks(self, block_id: int) -> List[Any]:
        """Get all items from block"""
        return self._blocks.get(block_id, [])

    def _set_blocks(self, block_id: int, values: List[Any]) -> None:
        """Set all block values"""
        if values:
            self._blocks[block_id] = values
        elif block_id in self._blocks:
            del self._blocks[block_id]

    def _remove_block(self, block_id: int) -> bool:
        """Remove block entirely"""
        if block_id in self._blocks:
            del self._blocks[block_id]
            return True
        return False

    @property
    def exists(self) -> bool:
        """Check if any relevant blocks exist."""
        return any(bid in self._blocks for bid in self.BLOCK_IDS)


class SgffListModel(SgffModel, Generic[T]):
    """Base for list-based block data."""

    def __init__(self, blocks: Dict[int, List[Any]]):
        super().__init__(blocks)
        self._items: Optional[List[T]] = None

    def _load(self) -> List[T]:
        """Load and parse items from block storage."""
        raise NotImplementedError

    def _sync(self) -> None:
        """Write items to block storage."""
        raise NotImplementedError

    @property
    def items(self) -> List[T]:
        """Get parsed items."""
        if self._items is None:
            self._items = self._load()
        return self._items

    def __iter__(self) -> Iterator[T]:
        return iter(self.items)

    def __len__(self) -> int:
        return len(self.items)

    def __getitem__(self, idx: int) -> T:
        return self.items[idx]

    def add(self, item: T) -> None:
        """Add item and sync."""
        self.items.append(item)
        self._sync()

    def remove(self, idx: int) -> bool:
        """Remove item by index and sync."""
        if 0 <= idx < len(self.items):
            self.items.pop(idx)
            self._sync()
            return True
        return False

    def clear(self) -> None:
        """Remove all items"""
        self._items = []
        self._sync()
