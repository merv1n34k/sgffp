"""
SnapGene file reader
"""

import struct
from io import BytesIO
from typing import Union, BinaryIO
from pathlib import Path

from .internal import SgffObject, Cookie
from .parsers import parse_blocks


class SgffReader:
    """Read and parse SnapGene files into SgffObject"""

    def __init__(self, source: Union[str, Path, BinaryIO]):
        if isinstance(source, (str, Path)):
            self.stream = open(source, "rb")
            self.should_close = True
        else:
            self.stream = source
            self.should_close = False

    def read(self) -> SgffObject:
        """Parse file and return SgffObject"""
        try:
            return self._parse()
        finally:
            if self.should_close:
                self.stream.close()

    def _parse(self) -> SgffObject:
        """Internal parsing logic"""
        if self.stream.read(1) != b"\t":
            raise ValueError("Invalid SnapGene file: wrong magic byte")

        length = struct.unpack(">I", self.stream.read(4))[0]
        title = self.stream.read(8)

        if length != 14 or title != b"SnapGene":
            raise ValueError("Invalid SnapGene file: wrong header")

        cookie = Cookie(
            type_of_sequence=struct.unpack(">H", self.stream.read(2))[0],
            export_version=struct.unpack(">H", self.stream.read(2))[0],
            import_version=struct.unpack(">H", self.stream.read(2))[0],
        )

        blocks = parse_blocks(self.stream)

        return SgffObject(cookie=cookie, blocks=blocks)

    @classmethod
    def from_file(cls, filepath: Union[str, Path]) -> SgffObject:
        """Read from file path"""
        return cls(filepath).read()

    @classmethod
    def from_bytes(cls, data: bytes) -> SgffObject:
        """Read from bytes"""
        return cls(BytesIO(data)).read()

