"""
Primer model (block 5)
"""

from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional, Iterator

from .base import SgffModel


@dataclass
class SgffPrimer:
    """Single primer definition"""

    name: str = ""
    sequence: str = ""
    bind_position: Optional[int] = None
    bind_strand: str = "+"
    _raw: Dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffPrimer":
        return cls(
            name=data.get("name", ""),
            sequence=data.get("sequence", ""),
            bind_position=data.get("bindingSite"),
            bind_strand=data.get("strand", "+"),
            _raw=data,
        )

    def to_dict(self) -> Dict:
        result = dict(self._raw)
        result.update({
            "name": self.name,
            "sequence": self.sequence,
        })
        if self.bind_position is not None:
            result["bindingSite"] = self.bind_position
        if self.bind_strand:
            result["strand"] = self.bind_strand
        return result


class SgffPrimerList(SgffModel):
    """Primer list wrapper for block 5"""

    BLOCK_IDS = (5,)

    def __init__(self, blocks: Dict[int, List[Any]]):
        super().__init__(blocks)
        self._items: Optional[List[SgffPrimer]] = None
        self._raw: Optional[Dict] = None

    def _load(self) -> List[SgffPrimer]:
        data = self._get_block(5)
        if not data:
            return []

        self._raw = data
        primers_data = data.get("Primers", {})
        if primers_data is None:
            return []

        primers = primers_data.get("Primer", [])
        if not isinstance(primers, list):
            primers = [primers] if primers else []

        return [SgffPrimer.from_dict(p) for p in primers]

    @property
    def items(self) -> List[SgffPrimer]:
        if self._items is None:
            self._items = self._load()
        return self._items

    def add(self, primer: SgffPrimer) -> None:
        self.items.append(primer)
        self._sync()

    def remove(self, idx: int) -> bool:
        if 0 <= idx < len(self.items):
            self.items.pop(idx)
            self._sync()
            return True
        return False

    def clear(self) -> None:
        self._items = []
        self._sync()

    def _sync(self) -> None:
        if self._items is None:
            return

        if self._items:
            data = {"Primers": {"Primer": [p.to_dict() for p in self._items]}}
            self._set_block(5, data)
        else:
            self._remove_block(5)

    def __iter__(self) -> Iterator[SgffPrimer]:
        return iter(self.items)

    def __len__(self) -> int:
        return len(self.items)

    def __getitem__(self, idx: int) -> SgffPrimer:
        return self.items[idx]

    def __repr__(self) -> str:
        return f"SgffPrimerList(count={len(self)})"
