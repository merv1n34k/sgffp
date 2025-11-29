"""
SnapGene file writer
"""

import struct
from typing import Union, BinaryIO
from pathlib import Path

from .internal import SgffObject


class SgffWriter:
    """
    Write SgffObject to SnapGene file format

    Can write to filepath or file-like object (stream)
    """

    def __init__(self, target: Union[str, Path, BinaryIO]):
        """
        Initialize writer with file path or stream

        Args:
            target: File path (str/Path) or file-like object with write()
        """
        if isinstance(target, (str, Path)):
            self.stream = open(target, "wb")
            self.should_close = True
        else:
            self.stream = target
            self.should_close = False

    def write(self, sgff: SgffObject) -> None:
        """Write SgffObject to file"""
        try:
            self._write_file(sgff)
        finally:
            if self.should_close:
                self.stream.close()

    def _write_file(self, sgff: SgffObject) -> None:
        """Internal writing logic"""
        # Write header
        self.stream.write(b"\t")
        self.stream.write(struct.pack(">I", 14))
        self.stream.write(b"SnapGene")

        # Write cookie
        self.stream.write(struct.pack(">H", sgff.cookie.type_of_sequence))
        self.stream.write(struct.pack(">H", sgff.cookie.export_version))
        self.stream.write(struct.pack(">H", sgff.cookie.import_version))

        # Write blocks in order
        for block_key in sorted(sgff.blocks.keys(), key=self._sort_key):
            block_type = int(block_key.split(".")[0])
            block_data = self._serialize_block(block_type, sgff.blocks[block_key])

            self.stream.write(bytes([block_type]))
            self.stream.write(struct.pack(">I", len(block_data)))
            self.stream.write(block_data)

    def _sort_key(self, key: str) -> tuple:
        """Sort blocks by type then index"""
        parts = key.split(".")
        block_type = int(parts[0])
        index = int(parts[1]) if len(parts) > 1 else 0
        return (block_type, index)

    def _serialize_block(self, block_type: int, data) -> bytes:
        """
        Convert parsed data back to binary format

        TODO: Implement full serialization for all block types
        Currently supports simple string/bytes blocks
        """
        if isinstance(data, str):
            return data.encode("utf-8")
        elif isinstance(data, bytes):
            return data
        elif isinstance(data, dict):
            # For complex blocks, would need type-specific serialization
            raise NotImplementedError(
                f"Serialization for block type {block_type} not yet implemented"
            )
        else:
            raise ValueError(f"Unknown data type for block {block_type}")

    @classmethod
    def to_file(cls, sgff: SgffObject, filepath: Union[str, Path]) -> None:
        """Convenience method to write to file path"""
        cls(filepath).write(sgff)

    @classmethod
    def to_bytes(cls, sgff: SgffObject) -> bytes:
        """Convenience method to write to bytes"""
        from io import BytesIO

        stream = BytesIO()
        cls(stream).write(sgff)
        return stream.getvalue()
