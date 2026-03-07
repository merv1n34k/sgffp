"""
Alignable sequences model
"""

from dataclasses import dataclass, field
from typing import Dict, List

from .base import SgffListModel

_ALIGNMENT_KNOWN_KEYS = frozenset({"name", "sequence"})


@dataclass
class SgffAlignment:
    """Single alignable sequence"""

    name: str = ""
    sequence: str = ""
    extras: Dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffAlignment":
        extras = {k: v for k, v in data.items() if k not in _ALIGNMENT_KNOWN_KEYS}
        return cls(
            name=data.get("name", ""),
            sequence=data.get("sequence", ""),
            extras=extras,
        )

    def to_dict(self) -> Dict:
        result = dict(self.extras)
        result["name"] = self.name
        result["sequence"] = self.sequence
        return result


class SgffAlignmentList(SgffListModel[SgffAlignment]):
    """Alignable sequences list with block management."""

    BLOCK_IDS = (17,)

    def __init__(self, blocks, **kwargs):
        super().__init__(blocks, **kwargs)
        self._wrapper_extras: Dict = {}

    def _load(self) -> List[SgffAlignment]:
        data = self._get_block(17)
        if not data:
            return []

        seqs_data = data.get("AlignableSequences", {})
        if seqs_data is None:
            return []

        # Capture wrapper-level extras (trimStringency, etc.)
        self._wrapper_extras = {
            k: v for k, v in seqs_data.items() if k != "Sequence"
        }

        seqs = seqs_data.get("Sequence", [])
        if not isinstance(seqs, list):
            seqs = [seqs] if seqs else []

        return [SgffAlignment.from_dict(s) for s in seqs]

    def _sync(self) -> None:
        if self._items is None:
            return
        if self._items:
            wrapper = dict(self._wrapper_extras)
            wrapper["Sequence"] = [a.to_dict() for a in self._items]
            self._set_block(17, {"AlignableSequences": wrapper})
        else:
            self._remove_block(17)

    def __repr__(self) -> str:
        return f"SgffAlignmentList(count={len(self)})"
