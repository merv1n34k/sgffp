"""
Primer model for oligonucleotide definitions
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .base import SgffListModel

_PRIMER_KNOWN_KEYS = frozenset({"name", "sequence", "bindingSite", "strand", "BindingSite"})


@dataclass
class SgffBindingSite:
    """A single primer binding site on the target sequence.

    SnapGene stores detailed and simplified versions of each binding site.
    The ``simplified`` flag distinguishes them.
    """

    start: int = 0
    end: int = 0
    bound_strand: str = "+"
    annealed_bases: str = ""
    melting_temperature: Optional[float] = None
    simplified: bool = False
    extras: Dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffBindingSite":
        extras = {}
        start, end = 0, 0
        location = data.get("location", "")
        if location and "-" in location:
            parts = location.split("-", 1)
            start = int(parts[0]) - 1  # convert to 0-based
            end = int(parts[1])

        raw_strand = data.get("boundStrand", "0")
        bound_strand = "-" if raw_strand == "1" else "+"

        mt = data.get("meltingTemperature")
        melting_temperature = float(mt) if mt is not None else None

        simplified = data.get("simplified", "0") == "1"

        # Everything not explicitly modeled goes to extras
        _known = {"location", "boundStrand", "meltingTemperature", "simplified", "annealedBases"}
        extras = {k: v for k, v in data.items() if k not in _known}

        return cls(
            start=start,
            end=end,
            bound_strand=bound_strand,
            annealed_bases=data.get("annealedBases", ""),
            melting_temperature=melting_temperature,
            simplified=simplified,
            extras=extras,
        )

    def to_dict(self) -> Dict:
        result = dict(self.extras)
        if self.simplified:
            result["simplified"] = "1"
        result["location"] = f"{self.start + 1}-{self.end}"  # back to 1-based
        result["boundStrand"] = "1" if self.bound_strand == "-" else "0"
        if self.annealed_bases:
            result["annealedBases"] = self.annealed_bases
        if self.melting_temperature is not None:
            result["meltingTemperature"] = str(int(self.melting_temperature)) if self.melting_temperature == int(self.melting_temperature) else str(self.melting_temperature)
        return result


@dataclass
class SgffPrimer:
    """Single primer definition"""

    name: str = ""
    sequence: str = ""
    binding_sites: List[SgffBindingSite] = field(default_factory=list)
    extras: Dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffPrimer":
        extras = {k: v for k, v in data.items() if k not in _PRIMER_KNOWN_KEYS}

        # Parse BindingSite child elements
        binding_sites = []
        raw_bs = data.get("BindingSite")
        if raw_bs is not None:
            if not isinstance(raw_bs, list):
                raw_bs = [raw_bs]
            binding_sites = [SgffBindingSite.from_dict(bs) for bs in raw_bs]

        return cls(
            name=data.get("name", ""),
            sequence=data.get("sequence", ""),
            binding_sites=binding_sites,
            extras=extras,
        )

    @property
    def bind_position(self) -> Optional[int]:
        """Position of the first non-simplified binding site (backward compat)."""
        for bs in self.binding_sites:
            if not bs.simplified:
                return bs.start
        return None

    @property
    def bind_strand(self) -> str:
        """Strand of the first non-simplified binding site (backward compat)."""
        for bs in self.binding_sites:
            if not bs.simplified:
                return bs.bound_strand
        return "+"

    @property
    def melting_temperature(self) -> Optional[float]:
        """Melting temperature of the first non-simplified binding site."""
        for bs in self.binding_sites:
            if not bs.simplified:
                return bs.melting_temperature
        return None

    def clear_binding_sites(self) -> None:
        """Remove all binding site data so SnapGene recalculates it."""
        self.binding_sites.clear()

    def to_dict(self) -> Dict:
        result = dict(self.extras)
        result["name"] = self.name
        result["sequence"] = self.sequence
        if self.binding_sites:
            sites = [bs.to_dict() for bs in self.binding_sites]
            result["BindingSite"] = sites if len(sites) > 1 else sites[0]
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
