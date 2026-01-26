"""
Feature model for annotations (block 10)
"""

from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional, Iterator

from .base import SgffModel


@dataclass
class SgffSegment:
    """Feature segment (location range)"""

    start: int
    end: int
    color: Optional[str] = None

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffSegment":
        range_str = data.get("range", "1-1")
        parts = sorted(int(x) for x in range_str.split("-"))
        return cls(
            start=parts[0] - 1,
            end=parts[1] if len(parts) > 1 else parts[0],
            color=data.get("color"),
        )

    def to_dict(self) -> Dict:
        result = {"range": f"{self.start + 1}-{self.end}"}
        if self.color:
            result["color"] = self.color
        return result


@dataclass
class SgffFeature:
    """Single annotation feature"""

    name: str
    type: str
    strand: str = "+"
    segments: List[SgffSegment] = field(default_factory=list)
    qualifiers: Dict[str, Any] = field(default_factory=dict)
    color: Optional[str] = None

    @property
    def start(self) -> int:
        if not self.segments:
            return 0
        return min(s.start for s in self.segments)

    @property
    def end(self) -> int:
        if not self.segments:
            return 0
        return max(s.end for s in self.segments)

    @property
    def length(self) -> int:
        return self.end - self.start

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffFeature":
        segments = []
        for seg in data.get("segments", []):
            segments.append(SgffSegment.from_dict(seg))

        return cls(
            name=data.get("name", ""),
            type=data.get("type", ""),
            strand=data.get("strand", "+"),
            segments=segments,
            qualifiers=data.get("qualifiers", {}),
            color=data.get("color"),
        )

    def to_dict(self) -> Dict:
        return {
            "name": self.name,
            "type": self.type,
            "strand": self.strand,
            "segments": [s.to_dict() for s in self.segments],
            "qualifiers": self.qualifiers,
            "color": self.color,
        }


class SgffFeatureList(SgffModel):
    """Feature list wrapper for block 10"""

    BLOCK_IDS = (10,)

    def __init__(self, blocks: Dict[int, List[Any]]):
        super().__init__(blocks)
        self._items: Optional[List[SgffFeature]] = None

    def _load(self) -> List[SgffFeature]:
        data = self._get_block(10)
        if not data:
            return []

        features = data.get("features", [])
        return [SgffFeature.from_dict(f) for f in features]

    @property
    def items(self) -> List[SgffFeature]:
        if self._items is None:
            self._items = self._load()
        return self._items

    def add(self, feature: SgffFeature) -> None:
        self.items.append(feature)
        self._sync()

    def remove(self, idx: int) -> bool:
        if 0 <= idx < len(self.items):
            self.items.pop(idx)
            self._sync()
            return True
        return False

    def find_by_name(self, name: str) -> Optional[SgffFeature]:
        for f in self.items:
            if f.name == name:
                return f
        return None

    def find_by_type(self, type_: str) -> List[SgffFeature]:
        return [f for f in self.items if f.type == type_]

    def clear(self) -> None:
        self._items = []
        self._sync()

    def _sync(self) -> None:
        if self._items is None:
            return

        if self._items:
            data = {"features": [f.to_dict() for f in self._items]}
            self._set_block(10, data)
        else:
            self._remove_block(10)

    def __iter__(self) -> Iterator[SgffFeature]:
        return iter(self.items)

    def __len__(self) -> int:
        return len(self.items)

    def __getitem__(self, idx: int) -> SgffFeature:
        return self.items[idx]

    def __repr__(self) -> str:
        return f"SgffFeatureList(count={len(self)})"
