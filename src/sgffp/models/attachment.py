"""
Attachment model for embedded file attachments (block 23)
"""

from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional

from .base import SgffListModel


@dataclass
class SgffAttachment:
    """Single file attachment embedded in a SnapGene file."""

    id: int = 0
    name: str = ""
    data: bytes = b""
    size: int = 0
    mtime: str = ""
    compressible: str = "0"
    extras: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_blocks(
        cls,
        file_block: Dict[str, Any],
        manifest_entry: Optional[Dict[str, Any]] = None,
    ) -> "SgffAttachment":
        """Create attachment from parsed file block + manifest entry."""
        file_id = file_block["id"]
        file_data = file_block["data"]

        name = ""
        size = len(file_data)
        mtime = ""
        compressible = "0"
        extras = {}

        if manifest_entry:
            name = manifest_entry.get("name", "")
            size = int(manifest_entry.get("size", size))
            mtime = manifest_entry.get("mtime", "")
            compressible = manifest_entry.get("compressible", "0")
            extras = {
                k: v
                for k, v in manifest_entry.items()
                if k not in ("id", "name", "size", "mtime", "compressible")
            }

        return cls(
            id=file_id,
            name=name,
            data=file_data,
            size=size,
            mtime=mtime,
            compressible=compressible,
            extras=extras,
        )

    def to_manifest_dict(self) -> Dict[str, Any]:
        """Convert to manifest XML entry dict."""
        d: Dict[str, Any] = {
            "id": str(self.id),
            "name": self.name,
            "size": str(self.size),
        }
        if self.mtime:
            d["mtime"] = self.mtime
        d["compressible"] = self.compressible
        d.update(self.extras)
        return d

    def to_file_block(self) -> Dict[str, Any]:
        """Convert to parsed file block dict."""
        return {"_type": "file", "id": self.id, "data": self.data}

    def __repr__(self) -> str:
        return f"SgffAttachment(id={self.id}, name={self.name!r}, size={self.size})"


class SgffAttachmentList(SgffListModel[SgffAttachment]):
    """Attachment list with block 23 management."""

    BLOCK_IDS = (23,)

    def _load(self) -> List[SgffAttachment]:
        """Load attachments from block 23 entries."""
        blocks = self._get_blocks(23)
        if not blocks:
            return []

        # Separate manifest from file data
        manifest_entries: Dict[int, Dict[str, Any]] = {}
        file_blocks: List[Dict[str, Any]] = []

        for block in blocks:
            if block.get("_type") == "manifest":
                manifest = block.get("manifest", {})
                files_wrapper = manifest.get("Files", {})
                if files_wrapper is None:
                    continue
                file_list = files_wrapper.get("File", [])
                if not isinstance(file_list, list):
                    file_list = [file_list]
                for entry in file_list:
                    fid = int(entry.get("id", 0))
                    manifest_entries[fid] = entry
            elif block.get("_type") == "file":
                file_blocks.append(block)

        # Correlate file data with manifest metadata
        attachments = []
        for fb in file_blocks:
            entry = manifest_entries.get(fb["id"])
            attachments.append(SgffAttachment.from_blocks(fb, entry))

        return attachments

    def _sync(self) -> None:
        """Write attachments back to block 23 storage."""
        if self._items is None:
            return
        if not self._items:
            self._remove_block(23)
            return

        result = []

        # File data blocks (one per attachment)
        for att in self._items:
            result.append(att.to_file_block())

        # Manifest block
        file_entries = [att.to_manifest_dict() for att in self._items]
        manifest_dict = {
            "_type": "manifest",
            "manifest": {"Files": {"File": file_entries}},
        }
        result.append(manifest_dict)

        self._set_blocks(23, result)

    def _next_id(self) -> int:
        """Get next available file ID."""
        if not self.items:
            return 1
        return max(a.id for a in self.items) + 1

    def add(self, item: SgffAttachment) -> None:
        """Add attachment with auto-assigned ID if needed."""
        if item.id == 0:
            item.id = self._next_id()
        if item.size == 0 and item.data:
            item.size = len(item.data)
        self.items.append(item)
        self._sync()

    def get_by_name(self, name: str) -> Optional[SgffAttachment]:
        """Find attachment by filename."""
        for att in self.items:
            if att.name == name:
                return att
        return None

    def get_by_id(self, file_id: int) -> Optional[SgffAttachment]:
        """Find attachment by file ID."""
        for att in self.items:
            if att.id == file_id:
                return att
        return None

    def __repr__(self) -> str:
        return f"SgffAttachmentList(count={len(self)})"
