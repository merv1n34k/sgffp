"""
Trace model for SnapGene sequence trace data (block 18)

ZTR format contains chromatogram data from Sanger sequencing.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional

from .base import SgffListModel


@dataclass
class SgffTraceClip:
    """Quality clip boundaries for trace data"""

    left: int = 0  # Left clip position (0-based)
    right: int = 0  # Right clip position

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffTraceClip":
        return cls(
            left=data.get("left", 0),
            right=data.get("right", 0),
        )

    def to_dict(self) -> Dict:
        return {"left": self.left, "right": self.right}


@dataclass
class SgffTraceSamples:
    """Trace sample intensities for each channel"""

    a: List[int] = field(default_factory=list)  # Adenine channel
    c: List[int] = field(default_factory=list)  # Cytosine channel
    g: List[int] = field(default_factory=list)  # Guanine channel
    t: List[int] = field(default_factory=list)  # Thymine channel

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffTraceSamples":
        return cls(
            a=data.get("A", []),
            c=data.get("C", []),
            g=data.get("G", []),
            t=data.get("T", []),
        )

    def to_dict(self) -> Dict:
        result = {}
        if self.a:
            result["A"] = self.a
        if self.c:
            result["C"] = self.c
        if self.g:
            result["G"] = self.g
        if self.t:
            result["T"] = self.t
        return result

    @property
    def length(self) -> int:
        """Number of sample points"""
        return max(len(self.a), len(self.c), len(self.g), len(self.t), 0)

    def __len__(self) -> int:
        return self.length


@dataclass
class SgffTrace:
    """Single sequence trace (chromatogram) data"""

    bases: str = ""
    positions: List[int] = field(default_factory=list)
    confidence: List[int] = field(default_factory=list)
    samples: Optional[SgffTraceSamples] = None
    clip: Optional[SgffTraceClip] = None
    text: Dict[str, str] = field(default_factory=dict)
    comments: List[str] = field(default_factory=list)

    @property
    def sequence(self) -> str:
        """Alias for bases"""
        return self.bases

    @property
    def length(self) -> int:
        """Number of bases in trace"""
        return len(self.bases)

    @property
    def sample_count(self) -> int:
        """Number of sample points in trace"""
        return len(self.samples) if self.samples else 0

    def get_metadata(self, key: str, default: str = "") -> str:
        """Get metadata value by key"""
        return self.text.get(key, default)

    def get_confidence_at(self, index: int) -> Optional[int]:
        """Get confidence score at base index"""
        if self.confidence and 0 <= index < len(self.confidence):
            return self.confidence[index]
        return None

    def get_position_at(self, index: int) -> Optional[int]:
        """Get sample position at base index"""
        if self.positions and 0 <= index < len(self.positions):
            return self.positions[index]
        return None

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffTrace":
        samples = None
        if "samples" in data:
            samples = SgffTraceSamples.from_dict(data["samples"])

        clip = None
        if "clip" in data:
            clip = SgffTraceClip.from_dict(data["clip"])

        return cls(
            bases=data.get("bases", ""),
            positions=data.get("positions", []),
            confidence=data.get("confidence", []),
            samples=samples,
            clip=clip,
            text=data.get("text", {}),
            comments=data.get("comments", []),
        )

    def to_dict(self) -> Dict:
        result: Dict[str, Any] = {}

        if self.bases:
            result["bases"] = self.bases
        if self.positions:
            result["positions"] = self.positions
        if self.confidence:
            result["confidence"] = self.confidence
        if self.samples:
            result["samples"] = self.samples.to_dict()
        if self.clip:
            result["clip"] = self.clip.to_dict()
        if self.text:
            result["text"] = self.text
        if self.comments:
            result["comments"] = self.comments

        return result

    def __len__(self) -> int:
        return self.length

    def __repr__(self) -> str:
        return f"SgffTrace(bases={self.length}, samples={self.sample_count})"


class SgffTraceList(SgffListModel[SgffTrace]):
    """Trace list wrapper for block 16 (trace container)"""

    BLOCK_IDS = (16,)

    def _load(self) -> List[SgffTrace]:
        """Load traces from block 16 containers"""
        traces = []
        for container in self._get_blocks(16):
            nested = container.get("blocks", {})
            for trace_data in nested.get(18, []):
                traces.append(SgffTrace.from_dict(trace_data))
        return traces

    def _sync(self) -> None:
        """Write traces back to block 16 containers"""
        if self._items is None:
            return
        if self._items:
            containers = [
                {"flags": 0, "blocks": {18: [t.to_dict()]}} for t in self._items
            ]
            self._set_blocks(16, containers)
        else:
            self._remove_block(16)

    def __repr__(self) -> str:
        return f"SgffTraceList(count={len(self)})"
