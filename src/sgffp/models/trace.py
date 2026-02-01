"""
Trace model for SnapGene sequence trace data (block 18)

ZTR format contains chromatogram data from Sanger sequencing.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional

from .base import SgffModel


@dataclass
class SgffTraceClip:
    """Quality clip boundaries for trace data"""

    left: int  # Left clip position (0-based)
    right: int  # Right clip position

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


class SgffTrace(SgffModel):
    """
    SnapGene sequence trace data (block 18).

    Contains chromatogram data from Sanger sequencing in ZTR format.
    """

    BLOCK_IDS = (18,)

    def __init__(self, blocks: Dict[int, List[Any]]):
        super().__init__(blocks)
        self._bases: Optional[str] = None
        self._positions: Optional[List[int]] = None
        self._confidence: Optional[List[int]] = None
        self._samples: Optional[SgffTraceSamples] = None
        self._clip: Optional[SgffTraceClip] = None
        self._text: Optional[Dict[str, str]] = None
        self._comments: Optional[List[str]] = None
        self._loaded: bool = False

    def _load(self) -> None:
        """Load and parse trace data from block 18"""
        if self._loaded:
            return

        data = self._get_block(18)
        if not data:
            self._loaded = True
            return

        self._bases = data.get("bases")
        self._positions = data.get("positions")
        self._confidence = data.get("confidence")
        self._text = data.get("text")
        self._comments = data.get("comments")

        if "samples" in data:
            self._samples = SgffTraceSamples.from_dict(data["samples"])

        if "clip" in data:
            self._clip = SgffTraceClip.from_dict(data["clip"])

        self._loaded = True

    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------

    @property
    def bases(self) -> Optional[str]:
        """Base calls (sequence) from trace"""
        self._load()
        return self._bases

    @bases.setter
    def bases(self, value: str) -> None:
        self._load()
        self._bases = value
        self._sync()

    @property
    def positions(self) -> Optional[List[int]]:
        """Base-to-sample position mapping (one per base)"""
        self._load()
        return self._positions

    @positions.setter
    def positions(self, value: List[int]) -> None:
        self._load()
        self._positions = value
        self._sync()

    @property
    def confidence(self) -> Optional[List[int]]:
        """Confidence/quality scores (one per base)"""
        self._load()
        return self._confidence

    @confidence.setter
    def confidence(self, value: List[int]) -> None:
        self._load()
        self._confidence = value
        self._sync()

    @property
    def samples(self) -> Optional[SgffTraceSamples]:
        """Trace sample intensities for ACGT channels"""
        self._load()
        return self._samples

    @samples.setter
    def samples(self, value: SgffTraceSamples) -> None:
        self._load()
        self._samples = value
        self._sync()

    @property
    def clip(self) -> Optional[SgffTraceClip]:
        """Quality clip boundaries"""
        self._load()
        return self._clip

    @clip.setter
    def clip(self, value: SgffTraceClip) -> None:
        self._load()
        self._clip = value
        self._sync()

    @property
    def text(self) -> Optional[Dict[str, str]]:
        """Metadata key-value pairs"""
        self._load()
        return self._text

    @text.setter
    def text(self, value: Dict[str, str]) -> None:
        self._load()
        self._text = value
        self._sync()

    @property
    def comments(self) -> Optional[List[str]]:
        """Free-text comments"""
        self._load()
        return self._comments

    @comments.setter
    def comments(self, value: List[str]) -> None:
        self._load()
        self._comments = value
        self._sync()

    # -------------------------------------------------------------------------
    # Convenience properties
    # -------------------------------------------------------------------------

    @property
    def sequence(self) -> Optional[str]:
        """Alias for bases"""
        return self.bases

    @property
    def length(self) -> int:
        """Number of bases in trace"""
        return len(self.bases) if self.bases else 0

    @property
    def sample_count(self) -> int:
        """Number of sample points in trace"""
        return len(self.samples) if self.samples else 0

    # -------------------------------------------------------------------------
    # Methods
    # -------------------------------------------------------------------------

    def get_metadata(self, key: str, default: str = "") -> str:
        """Get metadata value by key"""
        if self.text:
            return self.text.get(key, default)
        return default

    def set_metadata(self, key: str, value: str) -> None:
        """Set metadata value"""
        self._load()
        if self._text is None:
            self._text = {}
        self._text[key] = value
        self._sync()

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

    # -------------------------------------------------------------------------
    # Sync
    # -------------------------------------------------------------------------

    def _sync(self) -> None:
        """Write trace data back to block 18"""
        data: Dict[str, Any] = {}

        if self._bases:
            data["bases"] = self._bases
        if self._positions:
            data["positions"] = self._positions
        if self._confidence:
            data["confidence"] = self._confidence
        if self._samples:
            data["samples"] = self._samples.to_dict()
        if self._clip:
            data["clip"] = self._clip.to_dict()
        if self._text:
            data["text"] = self._text
        if self._comments:
            data["comments"] = self._comments

        if data:
            self._set_block(18, data)
        else:
            self._remove_block(18)

    # -------------------------------------------------------------------------
    # Dunder methods
    # -------------------------------------------------------------------------

    def __len__(self) -> int:
        return self.length

    def __repr__(self) -> str:
        return f"SgffTrace(bases={self.length}, samples={self.sample_count})"
