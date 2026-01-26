"""
Notes model (block 6)
"""

from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional

from .base import SgffModel


@dataclass
class SgffNotes:
    """File notes and metadata from block 6"""

    BLOCK_IDS = (6,)

    _blocks: Dict[int, List[Any]] = field(repr=False)
    _data: Optional[Dict] = field(default=None, repr=False)

    def __init__(self, blocks: Dict[int, List[Any]]):
        self._blocks = blocks
        self._data = None

    def _load(self) -> Dict:
        items = self._blocks.get(6, [])
        if items:
            data = items[0]
            return data.get("Notes", {}) if data else {}
        return {}

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

    def remove(self, key: str) -> bool:
        if key in self.data:
            del self.data[key]
            self._sync()
            return True
        return False

    @property
    def description(self) -> str:
        return self.get("Description", "")

    @description.setter
    def description(self, value: str) -> None:
        self.set("Description", value)

    @property
    def created(self) -> Optional[str]:
        return self.get("Created")

    @property
    def last_modified(self) -> Optional[str]:
        return self.get("LastModified")

    @property
    def exists(self) -> bool:
        return 6 in self._blocks

    def _sync(self) -> None:
        if self._data:
            self._blocks[6] = [{"Notes": self._data}]
        elif 6 in self._blocks:
            del self._blocks[6]

    def __repr__(self) -> str:
        return f"SgffNotes(keys={list(self.data.keys())})"
