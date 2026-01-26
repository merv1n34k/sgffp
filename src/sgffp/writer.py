"""
SnapGene file writer
"""

import struct
import lzma
from typing import Union, BinaryIO, Any, Dict
from pathlib import Path
from io import BytesIO

import xmltodict

from .internal import SgffObject

# Keys that should be XML attributes (need @ prefix for xmltodict)
XML_ATTR_KEYS = {
    "name", "type", "directionality", "range", "color", "text", "int",
    "nextValidID", "minContinuousMatchLen", "allowMismatch",
    "minMeltingTemperature", "showAdditionalFivePrimeMatches",
    "minimumFivePrimeAnnealing", "UTC"
}


def _to_xmltodict(obj: Any) -> Any:
    """Convert clean JSON back to xmltodict format (add @ prefix for attributes)"""
    if isinstance(obj, dict):
        result = {}
        for key, value in obj.items():
            # Add @ prefix for known attribute keys
            if key in XML_ATTR_KEYS:
                new_key = f"@{key}"
            elif key == "_text":
                new_key = "#text"
            else:
                new_key = key
            result[new_key] = _to_xmltodict(value)
        return result
    elif isinstance(obj, list):
        return [_to_xmltodict(item) for item in obj]
    else:
        return obj


class SgffWriter:
    """Write SgffObject to SnapGene file format"""

    def __init__(self, target: Union[str, Path, BinaryIO]):
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
        # Header
        self.stream.write(b"\t")
        self.stream.write(struct.pack(">I", 14))
        self.stream.write(b"SnapGene")

        # Cookie
        self.stream.write(struct.pack(">H", sgff.cookie.type_of_sequence))
        self.stream.write(struct.pack(">H", sgff.cookie.export_version))
        self.stream.write(struct.pack(">H", sgff.cookie.import_version))

        # Blocks sorted by type
        for block_type in sorted(sgff.blocks.keys()):
            for item in sgff.blocks[block_type]:
                block_data = self._serialize(block_type, item)
                if block_data is not None:
                    self.stream.write(bytes([block_type]))
                    self.stream.write(struct.pack(">I", len(block_data)))
                    self.stream.write(block_data)

    def _serialize(self, block_type: int, data: Any) -> bytes:
        """Serialize block data to bytes"""
        # Already bytes
        if isinstance(data, bytes):
            return data

        # String (legacy)
        if isinstance(data, str):
            return data.encode("utf-8")

        # Dict - type-specific serialization
        if isinstance(data, dict):
            return self._serialize_dict(block_type, data)

        raise ValueError(f"Cannot serialize {type(data)} for block {block_type}")

    def _serialize_dict(self, block_type: int, data: Dict) -> bytes:
        """Serialize dict data based on block type"""
        # Sequence blocks (0, 21, 32)
        if block_type in (0, 21, 32):
            return self._serialize_sequence(data)

        # Compressed DNA (1)
        if block_type == 1:
            return self._serialize_compressed_dna(data)

        # Features (10)
        if block_type == 10:
            return self._serialize_features(data)

        # History tree (7) - LZMA compressed XML
        if block_type == 7:
            return self._serialize_lzma_xml(data)

        # History node (11) - binary format
        if block_type == 11:
            return self._serialize_history_node(data)

        # History modifier (29) - LZMA compressed XML
        if block_type == 29:
            return self._serialize_lzma_xml(data)

        # History content (30) - LZMA compressed nested TLV
        if block_type == 30:
            return self._serialize_lzma_nested(data)

        # XML blocks
        if block_type in (5, 6, 8, 17):
            return self._serialize_xml(data)

        # Default: try XML conversion
        return self._serialize_xml(data)

    def _serialize_sequence(self, data: Dict) -> bytes:
        """Serialize uncompressed sequence with properties"""
        props = 0
        if data.get("topology") == "circular":
            props |= 0x01
        if data.get("strandedness") == "double":
            props |= 0x02
        if data.get("dam_methylated"):
            props |= 0x04
        if data.get("dcm_methylated"):
            props |= 0x08
        if data.get("ecoki_methylated"):
            props |= 0x10

        sequence = data.get("sequence", "")
        return bytes([props]) + sequence.encode("utf-8")

    def _serialize_compressed_dna(self, data: Dict) -> bytes:
        """Serialize compressed DNA with mystery bytes"""
        sequence = data.get("sequence", "")
        length = len(sequence)
        mystery = data.get("mystery", b"\x00" * 14)

        # Encode sequence to 2-bit
        encoded = self._dna_to_octet(sequence)

        # compressed_length = 4 (uncompressed_length) + 14 (mystery) + len(encoded)
        compressed_length = 4 + 14 + len(encoded)

        buf = BytesIO()
        buf.write(struct.pack(">I", compressed_length))
        buf.write(struct.pack(">I", length))
        buf.write(mystery)
        buf.write(encoded)

        return buf.getvalue()

    def _dna_to_octet(self, sequence: str) -> bytes:
        """Convert DNA sequence to 2-bit encoding"""
        base_map = {"G": 0, "A": 1, "T": 2, "C": 3}
        result = bytearray()

        for i in range(0, len(sequence), 4):
            byte = 0
            for j, shift in enumerate([6, 4, 2, 0]):
                if i + j < len(sequence):
                    base = sequence[i + j].upper()
                    byte |= base_map.get(base, 0) << shift
            result.append(byte)

        return bytes(result)

    def _serialize_features(self, data: Dict) -> bytes:
        """Serialize features back to XML"""
        features = data.get("features", [])
        if not features:
            return self._serialize_xml(data)

        # Convert back to XML structure
        xml_features = []
        for f in features:
            strand_rev = {".": "0", "+": "1", "-": "2", "=": "3"}

            xml_f = {
                "@name": f.get("name", ""),
                "@type": f.get("type", ""),
                "@directionality": strand_rev.get(f.get("strand", "."), "0"),
            }

            # Segments - convert back to xmltodict format
            if f.get("segments"):
                xml_f["Segment"] = _to_xmltodict(f["segments"])

            # Qualifiers
            quals = f.get("qualifiers", {})
            if quals:
                xml_f["Q"] = [
                    {
                        "@name": k,
                        "V": {"@text": str(v)}
                        if not isinstance(v, list)
                        else [{"@text": str(x)} for x in v],
                    }
                    for k, v in quals.items()
                ]

            xml_features.append(xml_f)

        xml_dict = {"Features": {"Feature": xml_features}}
        return xmltodict.unparse(xml_dict, full_document=False).encode("utf-8")

    def _serialize_xml(self, data: Dict) -> bytes:
        """Serialize dict to XML"""
        try:
            xml_data = _to_xmltodict(data)
            return xmltodict.unparse(xml_data, full_document=False).encode("utf-8")
        except:
            raise ValueError("Cannot serialize dict to XML")

    def _serialize_lzma_xml(self, data: Dict) -> bytes:
        """Serialize dict to LZMA-compressed XML (blocks 7, 29)"""
        xml_bytes = self._serialize_xml(data)
        return lzma.compress(xml_bytes)

    def _serialize_history_node(self, data: Dict) -> bytes:
        """Serialize history node (block 11) to binary format"""
        buf = BytesIO()

        node_index = data.get("node_index", 0)
        seq_type = data.get("sequence_type", 0)

        buf.write(struct.pack(">I", node_index))
        buf.write(bytes([seq_type]))

        if seq_type == 1:
            # Compressed DNA
            sequence = data.get("sequence", "")
            encoded = self._dna_to_octet(sequence)
            compressed_length = 4 + 14 + len(encoded)

            buf.write(struct.pack(">I", compressed_length))
            buf.write(struct.pack(">I", len(sequence)))

            mystery = data.get("mystery", b"\x00" * 14)
            if len(mystery) < 14:
                mystery = mystery + b"\x00" * (14 - len(mystery))
            buf.write(mystery[:14])
            buf.write(encoded)

        elif seq_type in (0, 21, 32):
            # Uncompressed sequence
            sequence = data.get("sequence", "")
            seq_bytes = sequence.encode("ascii", errors="ignore")
            buf.write(struct.pack(">I", len(seq_bytes)))
            buf.write(seq_bytes)

        # seq_type == 29: modifier only, no sequence data

        # Nested node_info (block 30 content)
        node_info = data.get("node_info")
        if node_info:
            lzma_data = self._serialize_lzma_nested(node_info)
            buf.write(bytes([0x1E]))  # Block type 30
            buf.write(struct.pack(">I", len(lzma_data)))
            buf.write(lzma_data)

        return buf.getvalue()

    def _serialize_lzma_nested(self, data: Dict) -> bytes:
        """Serialize nested TLV blocks with LZMA compression (block 30)"""
        buf = BytesIO()

        for block_type, items in data.items():
            if not isinstance(block_type, int):
                continue

            for item in items:
                block_data = self._serialize(block_type, item)
                if block_data:
                    buf.write(bytes([block_type]))
                    buf.write(struct.pack(">I", len(block_data)))
                    buf.write(block_data)

        return lzma.compress(buf.getvalue())

    @classmethod
    def to_file(cls, sgff: SgffObject, filepath: Union[str, Path]) -> None:
        """Write to file path"""
        cls(filepath).write(sgff)

    @classmethod
    def to_bytes(cls, sgff: SgffObject) -> bytes:
        """Write to bytes"""
        stream = BytesIO()
        cls(stream).write(sgff)
        return stream.getvalue()
