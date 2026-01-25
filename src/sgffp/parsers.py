"""
Block type parsers and parsing scheme
"""

import struct
import lzma
import zlib
from io import BytesIO
from typing import Dict, Tuple, Optional, Callable, Any, List

import xmltodict


def parse_blocks(stream) -> Dict[int, List[Any]]:
    """Parse multiple TLV blocks from stream into {type: [values]}"""
    result: Dict[int, List[Any]] = {}

    while True:
        block_type, block_length = read_header(stream)
        if block_type is None:
            break

        # Skip unknown blocks
        if block_type not in SCHEME:
            stream.read(block_length)
            continue

        length_override, parser = SCHEME[block_type]

        if length_override is not None:
            block_length = length_override

        data = stream.read(block_length)

        if parser is None:
            continue

        parsed = parser(data)
        if parsed is not None:
            if block_type not in result:
                result[block_type] = []
            result[block_type].append(parsed)

    return result


def read_header(stream) -> Tuple[Optional[int], Optional[int]]:
    """Read TLV header: 1 byte type + 4 bytes length"""
    type_byte = stream.read(1)
    if not type_byte:
        return None, None
    return type_byte[0], struct.unpack(">I", stream.read(4))[0]


def octet_to_dna(raw_data: bytes, base_count: int) -> bytes:
    """Convert 2-bit GATC encoding to ASCII"""
    bases = b"GATC"
    result = bytearray()
    for byte in raw_data:
        for shift in [6, 4, 2, 0]:
            if len(result) < base_count:
                result.append(bases[(byte >> shift) & 3])
    return bytes(result[:base_count])


# =============================================================================
# SEQUENCE PARSERS
# =============================================================================


def parse_sequence(data: bytes) -> Dict[str, Any]:
    """Type 0, 21, 32: Uncompressed sequence with properties"""
    props = data[0]
    sequence = data[1:].decode("utf-8", errors="ignore")

    return {
        "sequence": sequence,
        "length": len(sequence),
        "topology": "circular" if props & 0x01 else "linear",
        "strandedness": "double" if props & 0x02 else "single",
        "dam_methylated": bool(props & 0x04),
        "dcm_methylated": bool(props & 0x08),
        "ecoki_methylated": bool(props & 0x10),
    }


def parse_compressed_dna(data: bytes) -> Dict[str, Any]:
    """Type 1: Compressed DNA sequence with mystery bytes preserved"""
    offset = 0

    compressed_length = struct.unpack(">I", data[offset : offset + 4])[0]
    offset += 4

    uncompressed_length = struct.unpack(">I", data[offset : offset + 4])[0]
    offset += 4

    # Mystery bytes (14 bytes) - preserve for round-trip
    mystery = data[offset : offset + 10]
    offset += 14

    total_bytes = (uncompressed_length * 2 + 7) // 8
    seq_data = data[offset : offset + total_bytes]

    return {
        "sequence": octet_to_dna(seq_data, uncompressed_length).decode("ascii"),
        "length": uncompressed_length,
        "mystery": mystery,
    }


# =============================================================================
# XML PARSERS
# =============================================================================


def _clean_xml_dict(obj: Any) -> Any:
    """Convert xmltodict format to pure JSON (remove @ prefixes, handle #text)"""
    if isinstance(obj, dict):
        result = {}
        for key, value in obj.items():
            # Remove @ prefix from attribute keys
            clean_key = key[1:] if key.startswith("@") else key
            # Convert #text to _text
            if clean_key == "#text":
                clean_key = "_text"
            result[clean_key] = _clean_xml_dict(value)
        return result
    elif isinstance(obj, list):
        return [_clean_xml_dict(item) for item in obj]
    else:
        return obj


def parse_xml(data: bytes) -> Optional[Dict]:
    """Parse XML blocks into dict with clean JSON format"""
    try:
        xml_str = data.decode("utf-8", errors="ignore")
        parsed = xmltodict.parse(xml_str)
        return _clean_xml_dict(parsed)
    except:
        return None


def parse_lzma_xml(data: bytes) -> Optional[Dict]:
    """Parse LZMA-compressed XML into dict with clean JSON format"""
    try:
        decompressed = lzma.decompress(data)
        parsed = xmltodict.parse(decompressed.decode("utf-8", errors="ignore"))
        return _clean_xml_dict(parsed)
    except:
        return None


def parse_lzma_nested(data: bytes) -> Optional[Dict[int, List[Any]]]:
    """Type 30: LZMA with nested TLV blocks"""
    try:
        decompressed = lzma.decompress(data)
        return parse_blocks(BytesIO(decompressed))
    except:
        return None


# =============================================================================
# FEATURE PARSER
# =============================================================================

STRAND_MAP = {"0": ".", "1": "+", "2": "-", "3": "="}


def parse_features(data: bytes) -> Optional[Dict]:
    """Type 10: Features with full qualifier extraction"""
    parsed = parse_xml(data)
    if not parsed or "Features" not in parsed:
        return parsed

    # Handle empty <Features></Features> element (xmltodict returns None)
    if parsed["Features"] is None:
        return {"features": []}

    features_data = parsed["Features"].get("Feature", [])
    if not isinstance(features_data, list):
        features_data = [features_data]

    features = []
    for feature in features_data:
        segments = feature.get("Segment", [])
        if not isinstance(segments, list):
            segments = [segments]

        # Parse segment ranges
        ranges = []
        for seg in segments:
            if "range" in seg:
                r = sorted(int(x) for x in seg["range"].split("-"))
                ranges.append(r)

        # Parse qualifiers
        qualifiers = _parse_qualifiers(feature.get("Q", []))

        # Defaults
        if "label" not in qualifiers:
            qualifiers["label"] = feature.get("name", "")

        color = segments[0].get("color", "") if segments else ""

        features.append(
            {
                "name": feature.get("name", ""),
                "type": feature.get("type", ""),
                "strand": STRAND_MAP.get(feature.get("directionality", "0"), "."),
                "start": min(r[0] - 1 for r in ranges) if ranges else 0,
                "end": max(r[1] for r in ranges) if ranges else 0,
                "color": color,
                "segments": segments,
                "qualifiers": qualifiers,
            }
        )

    return {"features": features}


def _parse_qualifiers(quals: Any) -> Dict[str, Any]:
    """Extract qualifiers from feature"""
    if not quals:
        return {}
    if not isinstance(quals, list):
        quals = [quals]

    result = {}
    for q in quals:
        name = q.get("name", "")
        val = q.get("V")

        if val is None:
            continue
        elif isinstance(val, dict):
            result[name] = _extract_value(val)
        elif isinstance(val, list):
            result[name] = [_extract_value(v) for v in val]
        else:
            result[name] = val

    return result


def _extract_value(v: Dict) -> Any:
    """Extract typed value from qualifier"""
    if "text" in v:
        return v["text"]
    if "int" in v:
        return int(v["int"])
    # Return first value found
    for key, val in v.items():
        if key != "_text":
            return val
    return v


# =============================================================================
# ZTR PARSER
# =============================================================================

ZTR_MAGIC = b"\xaeZTR\r\n\x1a\n"


def parse_ztr(data: bytes) -> Optional[Dict[str, Any]]:
    """Type 18: Sequence trace (ZTR format)"""
    if len(data) < 10 or data[:8] != ZTR_MAGIC:
        return None

    result = {}
    offset = 10

    while offset + 12 <= len(data):
        chunk_type = data[offset : offset + 4].decode("ascii", errors="ignore").strip()
        meta_len = struct.unpack(">I", data[offset + 4 : offset + 8])[0]
        offset += 8 + meta_len

        if offset + 4 > len(data):
            break

        data_len = struct.unpack(">I", data[offset : offset + 4])[0]
        offset += 4

        if offset + data_len > len(data):
            break

        chunk_data = data[offset : offset + data_len]

        # Decompress if zlib compressed
        if chunk_data and chunk_data[0] == 2:
            try:
                chunk_data = b"\x00" + zlib.decompress(chunk_data[5:])
            except:
                pass

        # Parse chunks
        if chunk_type == "BASE" and chunk_data and chunk_data[0] == 0:
            result["bases"] = chunk_data[2:].decode("ascii", errors="ignore")

        elif chunk_type == "TEXT" and chunk_data and chunk_data[0] == 0:
            items = chunk_data[2:-2].split(b"\x00")
            text = {}
            for i in range(0, len(items) - 1, 2):
                key = items[i].decode("ascii", errors="ignore")
                val = items[i + 1].decode("ascii", errors="ignore")
                text[key] = val
            result["text"] = text

        elif chunk_type == "SMP4":
            trace_len = len(chunk_data) // 8
            samples = {}
            for i, base in enumerate(["A", "C", "G", "T"]):
                start = i * trace_len * 2
                trace = [
                    struct.unpack(">H", chunk_data[start + j : start + j + 2])[0]
                    for j in range(0, trace_len * 2, 2)
                ]
                samples[base] = trace
            result["samples"] = samples

        elif chunk_type == "CLIP" and len(chunk_data) >= 9:
            result["clip"] = {
                "left": struct.unpack(">I", chunk_data[1:5])[0],
                "right": struct.unpack(">I", chunk_data[5:9])[0],
            }

        offset += data_len

    return result if result else None


# =============================================================================
# HISTORY NODE PARSER
# =============================================================================


def parse_history_node(data: bytes) -> Dict[str, Any]:
    """Type 11: History node - delegates to other parsers"""
    node = {}
    offset = 0

    node["node_index"] = struct.unpack(">I", data[offset : offset + 4])[0]
    offset += 4

    seq_type = data[offset]
    node["sequence_type"] = seq_type
    offset += 1

    # Type 29: modifier only
    if seq_type == 29:
        if offset < len(data):
            nested = parse_blocks(BytesIO(data[offset:]))
            if nested:
                node["node_info"] = nested
        return node

    # Type 1: compressed DNA
    if seq_type == 1:
        compressed_length = struct.unpack(">I", data[offset : offset + 4])[0]
        compressed_start = offset + 4

        block_data = data[offset : compressed_start + compressed_length]
        result = parse_compressed_dna(block_data)
        if result:
            node.update(result)

        offset = compressed_start + compressed_length

    # Types 0, 21, 32: uncompressed
    elif seq_type in [0, 21, 32]:
        seq_length = struct.unpack(">I", data[offset : offset + 4])[0]
        offset += 4
        node["sequence"] = data[offset : offset + seq_length].decode(
            "ascii", errors="ignore"
        )
        node["length"] = seq_length
        offset += seq_length

    # Parse remaining nested blocks
    if offset < len(data):
        nested = parse_blocks(BytesIO(data[offset:]))
        if nested:
            node["node_info"] = nested

    return node


# =============================================================================
# PARSING SCHEME
# =============================================================================

# Format: block_type -> (length_override, parser_function)
# Only known blocks are included - unknown blocks are skipped
SCHEME: Dict[int, Tuple[Optional[int], Optional[Callable]]] = {
    0: (None, parse_sequence),  #       DNA
    1: (None, parse_compressed_dna),  # 2bit compressed DNA
    5: (None, parse_xml),  #            Primers
    6: (None, parse_xml),  #            Notes
    7: (None, parse_lzma_xml),  #       History tree
    8: (None, parse_xml),  #            Sequence properties
    10: (None, parse_features),  #      Features
    11: (None, parse_history_node),  #  History node container
    16: (4, None),  #                   Legacy trace - skip content
    17: (None, parse_xml),  #           Alignable sequences
    18: (None, parse_ztr),  #           Sequence trace
    21: (None, parse_sequence),  #      Protein
    29: (None, parse_lzma_xml),  #      History modifier
    30: (None, parse_lzma_nested),  #   History node content
    32: (None, parse_sequence),  #      RNA
}
