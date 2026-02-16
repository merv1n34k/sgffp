"""
Sequence properties model (block 8)
"""

from typing import Dict, List, Any, Optional

from .base import SgffModel


class SgffProperties(SgffModel):
    """Sequence properties from block 8"""

    BLOCK_IDS = (8,)

    def __init__(self, blocks: Dict[int, List[Any]]):
        super().__init__(blocks)
        self._data: Optional[Dict] = None

    def _load(self) -> Dict:
        return self._get_block(8) or {}

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

    def _sync(self) -> None:
        if self._data:
            self._set_block(8, self._data)
        else:
            self._remove_block(8)

    def __repr__(self) -> str:
        return f"SgffProperties(keys={list(self.data.keys())})"
