# src/parsers.py
"""
Block type parsers and parsing scheme
"""

import struct
import lzma
from io import BytesIO
from typing import Dict, Tuple, Optional, Callable, Any


# Forward declaration for recursive parsing
def parse_blocks(stream) -> Dict[str, Any]:
    """Parse multiple TLV blocks from stream"""
    result = {}
    block_counters = {}

    while True:
        block_type, block_length = read_header(stream)
        if block_type is None:
            break

        length_override, parser = SCHEME.get(block_type, (None, None))

        if length_override is not None:
            block_length = length_override

        data = stream.read(block_length)

        if parser is not None:
            parsed_data = parser(data)
            if parsed_data is not None:
                block_key = str(block_type)

                if block_key in result:
                    if block_key not in block_counters:
                        block_counters[block_key] = 1
                    else:
                        block_counters[block_key] += 1
                    block_key = f"{block_key}.{block_counters[block_key]}"

                result[block_key] = parsed_data

    return result


def read_header(stream):
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


# Block parsers (simplified implementations)
def parse_dna_sequence(data):
    """Type 0, 21, 32: Uncompressed sequence"""
    return data[1:].decode("utf-8", errors="ignore")


def parse_compressed_dna(data):
    """Type 1: Compressed DNA sequence"""
    offset = 0
    compressed_length = struct.unpack(">I", data[offset : offset + 4])[0]
    offset += 4
    uncompressed_length = struct.unpack(">I", data[offset : offset + 4])[0]
    offset += 4 + 4 + 10  # Skip mystery fields

    total_bytes = (uncompressed_length * 2 + 7) // 8
    seq_data = data[offset : offset + total_bytes]

    return {
        "sequence": octet_to_dna(seq_data, uncompressed_length).decode("ascii"),
        "length": uncompressed_length,
    }


def parse_xml(data):
    """Parse XML blocks"""
    return data.decode("utf-8", errors="ignore")


def parse_lzma_xml(data):
    """Parse LZMA-compressed XML"""
    try:
        return lzma.decompress(data).decode("utf-8", errors="ignore")
    except:
        return None


def parse_lzma_nested(data):
    """Type 30: LZMA with nested TLV blocks"""
    try:
        decompressed = lzma.decompress(data)
        return parse_blocks(BytesIO(decompressed))
    except:
        return None


def parse_history_node(data):
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


# Global parsing scheme: (length_override, parser_function)
SCHEME: Dict[int, Tuple[Optional[int], Optional[Callable]]] = {
    0: (None, parse_dna_sequence),
    1: (None, parse_compressed_dna),
    5: (None, parse_xml),
    6: (None, parse_xml),
    7: (None, parse_lzma_xml),
    8: (None, parse_xml),
    10: (None, parse_xml),
    11: (None, parse_history_node),
    14: (None, parse_xml),
    16: (4, None),  # Legacy trace - skip
    17: (None, parse_xml),
    21: (None, parse_dna_sequence),
    28: (None, parse_xml),
    29: (None, parse_lzma_xml),
    30: (None, parse_lzma_nested),
    32: (None, parse_dna_sequence),
}
