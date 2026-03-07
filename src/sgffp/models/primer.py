"""
Primer model for oligonucleotide definitions
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .base import SgffListModel

_PRIMER_KNOWN_KEYS = frozenset({"name", "sequence", "bindingSite", "strand"})


@dataclass
class SgffPrimer:
    """Single primer definition"""

    name: str = ""
    sequence: str = ""
    bind_position: Optional[int] = None
    bind_strand: str = "+"
    extras: Dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffPrimer":
        extras = {k: v for k, v in data.items() if k not in _PRIMER_KNOWN_KEYS}
        return cls(
            name=data.get("name", ""),
            sequence=data.get("sequence", ""),
            bind_position=data.get("bindingSite"),
            bind_strand=data.get("strand", "+"),
            extras=extras,
        )

    def clear_binding_sites(self) -> None:
        """Remove cached binding site data so SnapGene recalculates it."""
        self.extras.pop("BindingSite", None)

    def to_dict(self) -> Dict:
        result = dict(self.extras)
        result["name"] = self.name
        result["sequence"] = self.sequence
        if self.bind_position is not None:
            result["bindingSite"] = self.bind_position
        if self.bind_strand:
            result["strand"] = self.bind_strand
        return result


class SgffPrimerList(SgffListModel[SgffPrimer]):
    """Primer list with block management."""

    BLOCK_IDS = (5,)

    def __init__(self, blocks, **kwargs):
        super().__init__(blocks, **kwargs)
        self._wrapper_extras: Dict = {}

    def _load(self) -> List[SgffPrimer]:
        data = self._get_block(5)
        if not data:
            return []

        primers_data = data.get("Primers", {})
        if primers_data is None:
            return []

        # Capture wrapper-level extras (HybridizationParams, nextValidID, etc.)
        self._wrapper_extras = {
            k: v for k, v in primers_data.items() if k != "Primer"
        }

        primers = primers_data.get("Primer", [])
        if not isinstance(primers, list):
            primers = [primers] if primers else []

        return [SgffPrimer.from_dict(p) for p in primers]

    def _sync(self) -> None:
        if self._items is None:
            return
        if self._items:
            wrapper = dict(self._wrapper_extras)
            wrapper["Primer"] = [p.to_dict() for p in self._items]
            self._set_block(5, {"Primers": wrapper})
        else:
            self._remove_block(5)

    def __repr__(self) -> str:
        return f"SgffPrimerList(count={len(self)})"
