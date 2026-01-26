"""
Alignable sequences model (block 17)
"""

from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional, Iterator

from .base import SgffModel


@dataclass
class SgffAlignment:
    """Single alignable sequence"""

    name: str = ""
    sequence: str = ""
    _raw: Dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffAlignment":
        return cls(
            name=data.get("name", ""),
            sequence=data.get("sequence", ""),
            _raw=data,
        )

    def to_dict(self) -> Dict:
        result = dict(self._raw)
        result["name"] = self.name
        result["sequence"] = self.sequence
        return result


class SgffAlignmentList(SgffModel):
    """Alignable sequences wrapper for block 17"""

    BLOCK_IDS = (17,)

    def __init__(self, blocks: Dict[int, List[Any]]):
        super().__init__(blocks)
        self._items: Optional[List[SgffAlignment]] = None
        self._raw: Optional[Dict] = None

    def _load(self) -> List[SgffAlignment]:
        data = self._get_block(17)
        if not data:
            return []

        self._raw = data
        seqs_data = data.get("AlignableSequences", {})
        if seqs_data is None:
            return []

        seqs = seqs_data.get("Sequence", [])
        if not isinstance(seqs, list):
            seqs = [seqs] if seqs else []

        return [SgffAlignment.from_dict(s) for s in seqs]

    @property
    def items(self) -> List[SgffAlignment]:
        if self._items is None:
            self._items = self._load()
        return self._items

    def add(self, alignment: SgffAlignment) -> None:
        self.items.append(alignment)
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
            data = {"AlignableSequences": {"Sequence": [a.to_dict() for a in self._items]}}
            self._set_block(17, data)
        else:
            self._remove_block(17)

    def __iter__(self) -> Iterator[SgffAlignment]:
        return iter(self.items)

    def __len__(self) -> int:
        return len(self.items)

    def __getitem__(self, idx: int) -> SgffAlignment:
        return self.items[idx]

    def __repr__(self) -> str:
        return f"SgffAlignmentList(count={len(self)})"
