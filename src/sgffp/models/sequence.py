"""
Sequence model for DNA, RNA, and Protein data
"""

from dataclasses import dataclass
from typing import Dict, List, Any, Optional

from .base import SgffModel


@dataclass
class SgffSequenceData:
    """Sequence data container"""

    value: str = ""
    topology: str = "linear"
    strandedness: str = "single"
    dam_methylated: bool = False
    dcm_methylated: bool = False
    ecoki_methylated: bool = False

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffSequenceData":
        return cls(
            value=data.get("sequence", ""),
            topology=data.get("topology", "linear"),
            strandedness=data.get("strandedness", "single"),
            dam_methylated=data.get("dam_methylated", False),
            dcm_methylated=data.get("dcm_methylated", False),
            ecoki_methylated=data.get("ecoki_methylated", False),
        )

    def to_dict(self) -> Dict:
        return {
            "sequence": self.value,
            "topology": self.topology,
            "strandedness": self.strandedness,
            "dam_methylated": self.dam_methylated,
            "dcm_methylated": self.dcm_methylated,
            "ecoki_methylated": self.ecoki_methylated,
        }

    @property
    def length(self) -> int:
        return len(self.value)


_RNA_BLOCK_ID = 32

# SnapGene stores RNA on disk as DNA bytes (T/t), but the user-facing alphabet
# is U/u. Translate on the parse/sync boundary so the model is always correct.
_T_TO_U = str.maketrans("Tt", "Uu")
_U_TO_T = str.maketrans("Uu", "Tt")


class SgffSequence(SgffModel):
    """Sequence wrapper for DNA, RNA, and Protein blocks."""

    BLOCK_IDS = (0, 1, 21, 32)

    def __init__(self, blocks: Dict[int, List[Any]]):
        super().__init__(blocks)
        self._data: Optional[SgffSequenceData] = None
        self._block_id: Optional[int] = None

    def _detect_block(self) -> Optional[int]:
        """Find which sequence block exists"""
        for bid in self.BLOCK_IDS:
            if bid in self._blocks:
                return bid
        return None

    def _load(self) -> SgffSequenceData:
        """Load sequence from blocks"""
        self._block_id = self._detect_block()
        if self._block_id is None:
            return SgffSequenceData()

        data = self._get_block(self._block_id)
        if not data:
            return SgffSequenceData()

        loaded = SgffSequenceData.from_dict(data)
        if self._block_id == _RNA_BLOCK_ID:
            loaded.value = loaded.value.translate(_T_TO_U)
        return loaded

    @property
    def data(self) -> SgffSequenceData:
        if self._data is None:
            self._data = self._load()
        return self._data

    @property
    def value(self) -> str:
        return self.data.value

    @value.setter
    def value(self, seq: str) -> None:
        self.data.value = seq
        self._sync()

    @property
    def length(self) -> int:
        return self.data.length

    @property
    def topology(self) -> str:
        return self.data.topology

    @topology.setter
    def topology(self, value: str) -> None:
        self.data.topology = value
        self._sync()

    @property
    def strandedness(self) -> str:
        return self.data.strandedness

    @strandedness.setter
    def strandedness(self, value: str) -> None:
        self.data.strandedness = value
        self._sync()

    @property
    def is_circular(self) -> bool:
        return self.topology == "circular"

    @property
    def is_double_stranded(self) -> bool:
        return self.strandedness == "double"

    @property
    def block_id(self) -> Optional[int]:
        """Active block type for this sequence."""
        if self._block_id is None:
            self._block_id = self._detect_block()
        return self._block_id

    @block_id.setter
    def block_id(self, value: int) -> None:
        """Set the sequence block type."""
        if value not in self.BLOCK_IDS:
            raise ValueError(f"Invalid block id: {value}")

        old_id = self._block_id
        if old_id and old_id != value:
            self._remove_block(old_id)

        self._block_id = value
        self._sync()

    def _sync(self) -> None:
        """Write data to block storage."""
        if self._data is None:
            return

        bid = self._block_id or 0
        out = self._data.to_dict()
        if bid == _RNA_BLOCK_ID:
            # On-disk RNA storage uses T/t; the model holds U/u.
            out["sequence"] = out["sequence"].translate(_U_TO_T)
        self._set_block(bid, out)

    def __repr__(self) -> str:
        return f"SgffSequence(length={self.length}, topology={self.topology})"
