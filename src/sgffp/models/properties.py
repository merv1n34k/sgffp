"""
Sequence properties model (block 8)
"""

from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional

from .base import SgffModel


@dataclass
class SgffProperties:
    """Sequence properties from block 8"""

    BLOCK_IDS = (8,)

    _blocks: Dict[int, List[Any]] = field(repr=False)
    _data: Optional[Dict] = field(default=None, repr=False)

    def __init__(self, blocks: Dict[int, List[Any]]):
        self._blocks = blocks
        self._data = None

    def _load(self) -> Dict:
        items = self._blocks.get(8, [])
        return items[0] if items else {}

    @property
    def data(self) -> Dict:
        if self._data is None:
            self._data = self._load()
        return self._data

    def get(self, key: str, default: Any = None) -> Any:
        return self.data.get(key, default)

    def set(self, key: str, value: Any) -> None:
        self.data[key] = value
        self._sync()

    @property
    def exists(self) -> bool:
        return 8 in self._blocks

    def _sync(self) -> None:
        if self._data:
            self._blocks[8] = [self._data]
        elif 8 in self._blocks:
            del self._blocks[8]

    def __repr__(self) -> str:
        return f"SgffProperties(keys={list(self.data.keys())})"
