"""
SnapGene file writer
"""

import struct
import json
import lzma
from typing import Union, BinaryIO, Any, Dict
from pathlib import Path
from io import BytesIO

import xmltodict

from .internal import SgffObject

# Uppercase keys that are XML attributes, not child elements
_UPPERCASE_ATTRS = {"ID"}


def _qual_value_to_xml(v: object) -> dict:
    """Convert a qualifier value to xmltodict V element format."""
    if isinstance(v, int):
        return {"@int": str(v)}
    return {"@text": str(v)}


def _to_xmltodict(obj: Any) -> Any:
    """Convert clean JSON dict back to xmltodict format.

    SnapGene XML convention: element names start with uppercase,
    attribute names start with lowercase. When _text is present,
    all sibling keys are attributes regardless of case.
    """
    if isinstance(obj, dict):
        has_text = "_text" in obj
        result = {}
        for key, value in obj.items():
            if key == "_text":
                result["#text"] = _to_xmltodict(value)
            elif has_text:
                # Sibling of _text — always an attribute
                result[f"@{key}"] = _to_xmltodict(value)
            elif key[0:1].islower() or key in _UPPERCASE_ATTRS:
                result[f"@{key}"] = _to_xmltodict(value)
            else:
                result[key] = _to_xmltodict(value)
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

        # Trace container (16) - 4 byte flags + nested TLV blocks
        if block_type == 16:
            return self._serialize_trace_container(data)

        # Trace (18) - ZTR binary format
        if block_type == 18:
            return self._serialize_ztr(data)

        # RNA structure predictions (34) - LZMA JSON
        if block_type == 34:
            return self._serialize_lzma_json(data)

        # XML blocks with declaration (matches SnapGene behavior)
        if block_type in (5, 14, 28):
            return self._serialize_xml(data, xml_declaration=True)

        # XML blocks without declaration
        if block_type in (6, 8, 17, 20):
            return self._serialize_xml(data)

        # Default: try XML conversion
        return self._serialize_xml(data)

    def _serialize_sequence(self, data: Dict) -> bytes:
        """Serialize uncompressed sequence with property flags."""
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
        """Serialize compressed DNA block."""
        sequence = data.get("sequence", "")
        length = len(sequence)

        # Reconstruct 14-byte metadata header from decoded fields
        header = bytearray(14)
        header[0] = data.get("format_version", 30)
        header[4] = data.get("strandedness_flag", 1)
        struct.pack_into(">H", header, 8, data.get("property_flags", 1))
        struct.pack_into(">H", header, 12, data.get("header_seq_length", length))

        # Encode sequence to 2-bit
        encoded = self._dna_to_octet(sequence)

        # compressed_length = 4 (uncompressed_length) + 14 (header) + len(encoded)
        compressed_length = 4 + 14 + len(encoded)

        buf = BytesIO()
        buf.write(struct.pack(">I", compressed_length))
        buf.write(struct.pack(">I", length))
        buf.write(bytes(header))
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
        """Serialize features to XML."""
        features = data.get("features", [])
        if not features:
            xml_dict = {"Features": None}
            xml_str = xmltodict.unparse(
                xml_dict, full_document=False, short_empty_elements=True
            )
            return ('<?xml version="1.0"?>' + xml_str).encode("utf-8")

        strand_rev = {".": "0", "+": "1", "-": "2", "=": "3"}
        wrapper_extras = data.get("wrapper_extras", {})

        xml_features = []
        for f in features:
            # Start from feature-level extras (unmodeled XML attrs)
            xml_f = _to_xmltodict(f.get("extras", {}))

            # Named feature attrs
            xml_f["@name"] = f.get("name", "")
            xml_f["@type"] = f.get("type", "")
            xml_f["@directionality"] = strand_rev.get(f.get("strand", "."), "0")

            # Segments
            if f.get("segments"):
                xml_f["Segment"] = _to_xmltodict(f["segments"])

            # Qualifiers: prefer raw_qualifiers for lossless roundtrip
            raw_quals = f.get("raw_qualifiers")
            if raw_quals:
                xml_f["Q"] = _to_xmltodict(raw_quals)
            else:
                quals = f.get("qualifiers", {})
                if quals:
                    xml_f["Q"] = self._qualifiers_to_xml(quals)

            xml_features.append(xml_f)

        features_wrapper = _to_xmltodict(wrapper_extras)
        features_wrapper["Feature"] = xml_features

        xml_dict = {"Features": features_wrapper}
        xml_str = xmltodict.unparse(
            xml_dict, full_document=False, short_empty_elements=True
        )
        return ('<?xml version="1.0"?>' + xml_str).encode("utf-8")

    @staticmethod
    def _qualifiers_to_xml(quals: Dict) -> list:
        """Convert qualifiers dict to xmltodict Q list format."""
        result = []
        for k, v in quals.items():
            if isinstance(v, list):
                result.append({
                    "@name": k,
                    "V": [_qual_value_to_xml(x) for x in v],
                })
            else:
                result.append({
                    "@name": k,
                    "V": _qual_value_to_xml(v),
                })
        return result

    def _serialize_xml(
        self, data: Dict, *, xml_declaration: bool = False
    ) -> bytes:
        """Serialize dict to XML."""
        try:
            xml_data = _to_xmltodict(data)
            xml_str = xmltodict.unparse(
                xml_data, full_document=False, short_empty_elements=True
            )
            if xml_declaration:
                xml_str = '<?xml version="1.0"?>' + xml_str
            return xml_str.encode("utf-8")
        except Exception:
            raise ValueError("Cannot serialize dict to XML")

    def _serialize_lzma_xml(self, data: Dict) -> bytes:
        """Serialize dict to LZMA-compressed XML."""
        xml_data = _to_xmltodict(data)
        xml_str = xmltodict.unparse(
            xml_data, full_document=False, short_empty_elements=True
        )
        xml_str = '<?xml version="1.0" encoding="UTF-8"?>' + xml_str
        return lzma.compress(xml_str.encode("utf-8"))

    def _serialize_lzma_json(self, data: Any) -> bytes:
        """Serialize to LZMA-compressed JSON."""
        json_bytes = json.dumps(data, separators=(",", ":")).encode("utf-8")
        return lzma.compress(json_bytes)

    def _serialize_history_node(self, data: Dict) -> bytes:
        """Serialize history node to binary format."""
        buf = BytesIO()

        node_index = data.get("node_index", 0)
        seq_type = data.get("sequence_type", 0)

        buf.write(struct.pack(">I", node_index))
        buf.write(bytes([seq_type]))

        if seq_type == 1:
            # Compressed DNA with decoded metadata header
            sequence = data.get("sequence", "")
            encoded = self._dna_to_octet(sequence)
            compressed_length = 4 + 14 + len(encoded)

            buf.write(struct.pack(">I", compressed_length))
            buf.write(struct.pack(">I", len(sequence)))

            # Reconstruct 14-byte metadata header
            header = bytearray(14)
            header[0] = data.get("format_version", 30)
            header[4] = data.get("strandedness_flag", 1)
            struct.pack_into(">H", header, 8, data.get("property_flags", 1))
            struct.pack_into(
                ">H", header, 12, data.get("header_seq_length", len(sequence))
            )
            buf.write(bytes(header))
            buf.write(encoded)

        elif seq_type in (0, 21, 32):
            # Uncompressed sequence
            sequence = data.get("sequence", "")
            seq_bytes = sequence.encode("ascii", errors="ignore")
            buf.write(struct.pack(">I", len(seq_bytes)))
            buf.write(seq_bytes)

        # seq_type == 29: modifier only, no sequence data

        # Nested node_info — bare TLV blocks (not LZMA-wrapped)
        node_info = data.get("node_info")
        if node_info:
            for block_type, items in node_info.items():
                if not isinstance(block_type, int):
                    continue
                for item in items:
                    block_data = self._serialize(block_type, item)
                    if block_data:
                        buf.write(bytes([block_type]))
                        buf.write(struct.pack(">I", len(block_data)))
                        buf.write(block_data)

        return buf.getvalue()

    def _serialize_trace_container(self, data: Dict) -> bytes:
        """Serialize trace container with flags and nested TLV blocks."""
        buf = BytesIO()

        # 4-byte flags header
        flags = data.get("flags", 0)
        buf.write(struct.pack(">I", flags))

        # Nested blocks (typically block 18 trace + optional block 8 properties)
        nested = data.get("blocks", {})
        for block_type in sorted(nested.keys()):
            if not isinstance(block_type, int):
                continue
            for item in nested[block_type]:
                block_data = self._serialize(block_type, item)
                if block_data:
                    buf.write(bytes([block_type]))
                    buf.write(struct.pack(">I", len(block_data)))
                    buf.write(block_data)

        return buf.getvalue()

    def _serialize_lzma_nested(self, data: Dict) -> bytes:
        """Serialize nested TLV blocks with LZMA compression."""
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

    def _serialize_ztr(self, data: Dict) -> bytes:
        """Serialize trace data to ZTR format."""
        buf = BytesIO()

        # ZTR magic and version
        buf.write(b"\xaeZTR\r\n\x1a\n")  # Magic
        buf.write(struct.pack(">BB", 1, 2))  # Version 1.2

        # Helper to write a chunk
        def write_chunk(chunk_type: str, chunk_data: bytes, metadata: bytes = b""):
            # Type (4 bytes, space-padded)
            type_bytes = chunk_type.encode("ascii").ljust(4)[:4]
            buf.write(type_bytes)
            # Metadata length + metadata
            buf.write(struct.pack(">I", len(metadata)))
            if metadata:
                buf.write(metadata)
            # Data length + data
            buf.write(struct.pack(">I", len(chunk_data)))
            buf.write(chunk_data)

        # BASE chunk: format byte (0) + padding (1) + ASCII bases
        if data.get("bases"):
            bases = data["bases"].encode("ascii")
            chunk_data = b"\x00\x00" + bases
            write_chunk("BASE", chunk_data)

        # BPOS chunk: format byte (0) + 3 padding + 4-byte positions
        if data.get("positions"):
            positions = data["positions"]
            chunk_data = b"\x00\x00\x00\x00"
            for pos in positions:
                chunk_data += struct.pack(">I", pos)
            write_chunk("BPOS", chunk_data)

        # CNF4 chunk: format byte (0) + confidence values (1 byte per base)
        if data.get("confidence"):
            confidence = data["confidence"]
            chunk_data = b"\x00" + bytes(confidence)
            write_chunk("CNF4", chunk_data)

        # SMP4 chunk: format byte (0) + padding + channel data
        if data.get("samples"):
            samples = data["samples"]
            # Check if all channels present (use SMP4), otherwise use SAMP
            if all(ch in samples for ch in ["A", "C", "G", "T"]):
                chunk_data = b"\x00\x00"
                for channel in ["A", "C", "G", "T"]:
                    for val in samples.get(channel, []):
                        chunk_data += struct.pack(">H", val)
                write_chunk("SMP4", chunk_data)
            else:
                # Individual SAMP chunks per channel
                for channel in ["A", "C", "G", "T"]:
                    if samples.get(channel):
                        chunk_data = b"\x00\x00"
                        for val in samples[channel]:
                            chunk_data += struct.pack(">H", val)
                        # Metadata is 4-byte channel name
                        metadata = channel.encode("ascii") + b"\x00\x00\x00"
                        write_chunk("SAMP", chunk_data, metadata)

        # TEXT chunk: format byte (0) + padding + null-terminated key-value pairs
        if data.get("text"):
            text_data = b"\x00\x00"
            for key, val in data["text"].items():
                text_data += key.encode("ascii") + b"\x00"
                text_data += str(val).encode("ascii") + b"\x00"
            write_chunk("TEXT", text_data)

        # CLIP chunk: format byte (0) + left (4) + right (4)
        if data.get("clip"):
            clip = data["clip"]
            chunk_data = b"\x00"
            chunk_data += struct.pack(">I", clip.get("left", 0))
            chunk_data += struct.pack(">I", clip.get("right", 0))
            write_chunk("CLIP", chunk_data)

        # COMM chunks: format byte (0) + free text
        if data.get("comments"):
            for comment in data["comments"]:
                chunk_data = b"\x00" + comment.encode("ascii")
                write_chunk("COMM", chunk_data)

        return buf.getvalue()

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
