"""
Notes model (block 6)
"""

from typing import Dict, List, Any, Optional

from .base import SgffModel


class SgffNotes(SgffModel):
    """File notes and metadata from block 6"""

    BLOCK_IDS = (6,)

    def __init__(self, blocks: Dict[int, List[Any]]):
        super().__init__(blocks)
        self._data: Optional[Dict] = None

    def _load(self) -> Dict:
        data = self._get_block(6)
        if data:
            return data.get("Notes", {})
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

    def _sync(self) -> None:
        if self._data:
            self._set_block(6, {"Notes": self._data})
        else:
            self._remove_block(6)

    def __repr__(self) -> str:
        return f"SgffNotes(keys={list(self.data.keys())})"
