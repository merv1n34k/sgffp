"""
Sequence model for DNA/RNA/Protein (blocks 0, 1, 21, 32)
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
    _mystery: Optional[bytes] = None

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffSequenceData":
        return cls(
            value=data.get("sequence", ""),
            topology=data.get("topology", "linear"),
            strandedness=data.get("strandedness", "single"),
            dam_methylated=data.get("dam_methylated", False),
            dcm_methylated=data.get("dcm_methylated", False),
            ecoki_methylated=data.get("ecoki_methylated", False),
            _mystery=data.get("mystery"),
        )

    def to_dict(self) -> Dict:
        result = {
            "sequence": self.value,
            "topology": self.topology,
            "strandedness": self.strandedness,
            "dam_methylated": self.dam_methylated,
            "dcm_methylated": self.dcm_methylated,
            "ecoki_methylated": self.ecoki_methylated,
        }
        if self._mystery:
            result["mystery"] = self._mystery
        return result

    @property
    def length(self) -> int:
        return len(self.value)


class SgffSequence(SgffModel):
    """
    Sequence wrapper for blocks 0 (DNA), 1 (compressed), 21 (protein), 32 (RNA).
    """

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
        if data:
            return SgffSequenceData.from_dict(data)
        return SgffSequenceData()

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
        """Which block type is used (0, 1, 21, or 32)"""
        if self._block_id is None:
            self._block_id = self._detect_block()
        return self._block_id

    @block_id.setter
    def block_id(self, value: int) -> None:
        """Change block type (e.g., switch from compressed to uncompressed)"""
        if value not in self.BLOCK_IDS:
            raise ValueError(f"Invalid block id: {value}")

        old_id = self._block_id
        if old_id and old_id != value:
            self._remove_block(old_id)

        self._block_id = value
        self._sync()

    def _sync(self) -> None:
        """Write data back to blocks"""
        if self._data is None:
            return

        bid = self._block_id or 0
        self._set_block(bid, self._data.to_dict())

    def __repr__(self) -> str:
        return f"SgffSequence(length={self.length}, topology={self.topology})"
