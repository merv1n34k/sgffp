"""
Alignable sequences model (block 17)
"""

from dataclasses import dataclass, field
from typing import Dict, List

from .base import SgffListModel


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


class SgffAlignmentList(SgffListModel[SgffAlignment]):
    """Alignable sequences wrapper for block 17"""

    BLOCK_IDS = (17,)

    def _load(self) -> List[SgffAlignment]:
        data = self._get_block(17)
        if not data:
            return []

        seqs_data = data.get("AlignableSequences", {})
        if seqs_data is None:
            return []

        seqs = seqs_data.get("Sequence", [])
        if not isinstance(seqs, list):
            seqs = [seqs] if seqs else []

        return [SgffAlignment.from_dict(s) for s in seqs]

    def _sync(self) -> None:
        if self._items is None:
            return
        if self._items:
            self._set_block(
                17,
                {
                    "AlignableSequences": {
                        "Sequence": [a.to_dict() for a in self._items]
                    }
                },
            )
        else:
            self._remove_block(17)

    def __repr__(self) -> str:
        return f"SgffAlignmentList(count={len(self)})"
