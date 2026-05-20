"""
Trace alignment model (block 27) — BGZF-compressed BAM data
"""

from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional

from .base import SgffModel


@dataclass
class SgffBamReference:
    """BAM reference sequence entry."""

    name: str = ""
    length: int = 0

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SgffBamReference":
        return cls(name=data.get("name", ""), length=data.get("length", 0))

    def to_dict(self) -> Dict[str, Any]:
        return {"name": self.name, "length": self.length}


@dataclass
class SgffBamRecord:
    """Single BAM alignment record."""

    read_name: str = ""
    flag: int = 0
    ref_id: int = 0
    pos: int = 0
    mapq: int = 255
    cigar: str = ""
    sequence: str = ""
    quality: List[int] = field(default_factory=list)
    next_ref_id: int = 0
    next_pos: int = 0
    tlen: int = 0
    bin: int = 0
    aux: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_unmapped(self) -> bool:
        return bool(self.flag & 0x4)

    @property
    def is_reverse(self) -> bool:
        return bool(self.flag & 0x10)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SgffBamRecord":
        return cls(
            read_name=data.get("read_name", ""),
            flag=data.get("flag", 0),
            ref_id=data.get("ref_id", 0),
            pos=data.get("pos", 0),
            mapq=data.get("mapq", 255),
            cigar=data.get("cigar", ""),
            sequence=data.get("sequence", ""),
            quality=data.get("quality", []),
            next_ref_id=data.get("next_ref_id", 0),
            next_pos=data.get("next_pos", 0),
            tlen=data.get("tlen", 0),
            bin=data.get("bin", 0),
            aux=data.get("aux", {}),
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "read_name": self.read_name,
            "flag": self.flag,
            "ref_id": self.ref_id,
            "pos": self.pos,
            "mapq": self.mapq,
            "cigar": self.cigar,
            "sequence": self.sequence,
            "quality": self.quality,
            "next_ref_id": self.next_ref_id,
            "next_pos": self.next_pos,
            "tlen": self.tlen,
            "bin": self.bin,
            "aux": self.aux,
        }


class SgffTraceAlignment(SgffModel):
    """Trace alignment data (block 27) — BAM alignments of traces against reference."""

    BLOCK_IDS = (27,)

    def __init__(self, blocks: Dict[int, List[Any]]):
        super().__init__(blocks)
        self._data: Optional[Dict] = None

    def _load(self) -> Dict:
        data = self._get_block(27)
        if data:
            return data
        return {"header": "", "references": [], "records": []}

    @property
    def data(self) -> Dict:
        if self._data is None:
            self._data = self._load()
        return self._data

    @property
    def header(self) -> str:
        return self.data.get("header", "")

    @property
    def references(self) -> List[SgffBamReference]:
        return [SgffBamReference.from_dict(r) for r in self.data.get("references", [])]

    @property
    def records(self) -> List[SgffBamRecord]:
        return [SgffBamRecord.from_dict(r) for r in self.data.get("records", [])]

    @property
    def record_count(self) -> int:
        return len(self.data.get("records", []))

    @property
    def reference_count(self) -> int:
        return len(self.data.get("references", []))

    def _sync(self) -> None:
        if self._data:
            self._set_block(27, self._data)
        else:
            self._remove_block(27)

    def __repr__(self) -> str:
        return f"SgffTraceAlignment(references={self.reference_count}, records={self.record_count})"
