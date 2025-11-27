#!/usr/bin/env python3
"""
Parse SnapGeneFileFormat (SGFF) file to JSON using recursive TLV parsing
"""

import struct
import json
import lzma
import zlib
import sys
from io import BytesIO


# =============================================================================
# SECTION 1: COMMON FUNCTIONS
# =============================================================================


def read_header(stream):
    """Read TLV block header: 1 byte type + 4 bytes length"""
    type_byte = stream.read(1)
    if not type_byte:
        return None, None

    block_type = type_byte[0]
    block_length = struct.unpack(">I", stream.read(4))[0]

    return block_type, block_length


def parse_blocks(stream):
    """Generic TLV block parser - parses multiple blocks from stream using SCHEME"""
    result = {}
    block_counters = {}  # Track how many times each block type is repeated

    while True:
        block_type, block_length = read_header(stream)
        if block_type is None:
            break

        # Get scheme for this block type (default to (None, None) for unknown types)
        length_override, parser = SCHEME.get(block_type, (None, None))

        # Apply length override BEFORE reading data
        if length_override is not None:
            block_length = length_override

        # Read block data
        data = stream.read(block_length)

        # Parse and store if parser exists and returns non-None
        if parser is not None:
            parsed_data = parser(data)
            if parsed_data is not None:
                block_key = str(block_type)

                # Handle repeated blocks with .1, .2, .3 etc suffixes
                if block_key in result:
                    # Block already exists - add numbered suffix
                    if block_key not in block_counters:
                        block_counters[block_key] = 1  # First repeat
                    else:
                        block_counters[block_key] += 1

                    block_key = f"{block_key}.{block_counters[block_key]}"

                result[block_key] = parsed_data

    return result


def decode_xml(data):
    """Decode data as UTF-8 XML text"""
    try:
        return data.decode("utf-8", errors="ignore")
    except:
        return data.hex()


def decode_lzma(data):
    """Decompress LZMA data and decode as UTF-8"""
    try:
        decompressed = lzma.decompress(data)
        return decompressed.decode("utf-8", errors="ignore")
    except:
        return data.hex()


def octet_to_dna(raw_data, base_count):
    """Convert 2-bit encoded DNA (GATC) to ASCII"""
    bases = b"GATC"
    result = bytearray()

    for byte in raw_data:
        for shift in [6, 4, 2, 0]:
            if len(result) < base_count:
                result.append(bases[(byte >> shift) & 3])
            else:
                break
        if len(result) >= base_count:
            break

    return bytes(result[:base_count])


# =============================================================================
# SECTION 2: BLOCK TYPE PARSERS
# =============================================================================


def parse_block_0(data):
    """Type 0: DNA sequence (uncompressed ASCII)"""
    try:
        return data[1:].decode("utf-8", errors="strict")
    except:
        return data.hex()


def parse_block_1(data):
    """Type 1: Compressed DNA sequence"""
    node = {}
    offset = 0

    # Parse compression header (12 bytes total)
    compressed_length = struct.unpack(">I", data[offset : offset + 4])[0]
    offset += 4

    uncompressed_length = struct.unpack(">I", data[offset : offset + 4])[0]
    offset += 4

    # TODO: INVESTIGATE - 4 bytes unknown flag/property
    # Possible meanings: topology, edit flags, quality markers
    offset += 4

    node["compressed_length"] = compressed_length
    node["uncompressed_length"] = uncompressed_length

    # TODO: INVESTIGATE - 10 bytes unknown header data
    # Structure: 1B + 4B + 1B + 4B (always: 01 + flag + 01 + length_repeat)
    offset += 10

    # Extract and decode compressed DNA sequence
    total_bytes_needed = (uncompressed_length * 2 + 7) // 8
    seq_data = data[offset : offset + total_bytes_needed]
    node["sequence"] = octet_to_dna(seq_data, uncompressed_length).decode("ascii")

    return node


def parse_block_3(data):
    """Type 3: Enzyme cutters"""
    if len(data) < 5 or data[0] != 1:
        return data.hex()

    csv_length = struct.unpack(">I", data[1:5])[0]
    if len(data) < 5 + csv_length:
        return data.hex()

    csv_data = data[5 : 5 + csv_length]
    sequences = csv_data.decode("ascii", errors="ignore").split(",")

    # values = []
    # offset = 5 + csv_length
    # while offset + 4 <= len(data):
    #    val = struct.unpack(">I", data[offset : offset + 4])[0]
    #    values.append(val)
    #    offset += 4

    return {"sequences": sequences}  # , "values": values}


def parse_block_5(data):
    """Type 5: Primers XML"""
    return decode_xml(data)


def parse_block_6(data):
    """Type 6: Notes XML"""
    return decode_xml(data)


def parse_block_7(data):
    """Type 7: History tree (LZMA compressed XML)"""
    return decode_lzma(data)


def parse_block_8(data):
    """Type 8: Additional sequence properties XML"""
    return decode_xml(data)


def parse_block_10(data):
    """Type 10: Features XML"""
    return decode_xml(data)


def parse_block_11(data):
    """Type 11: History node - organizes data and calls other parsers"""
    node = {}
    offset = 0

    # Parse node index (4 bytes)
    node["node_index"] = struct.unpack(">I", data[offset : offset + 4])[0]
    offset += 4

    # Parse sequence type (1 byte)
    seq_type = data[offset]
    node["sequence_type"] = seq_type
    offset += 1

    # Type 29: Modifier only (no sequence data, just type 0x1E block)
    if seq_type == 29:
        if offset < len(data):
            # Create stream and parse remaining blocks
            stream = BytesIO(data[offset:])
            nested = parse_blocks(stream)
            if nested:
                node["node_info"] = nested
        return node

    # Type 1: Compressed DNA
    if seq_type == 1:
        # Read compressed_length to know where block 1 ends
        compressed_length = struct.unpack(">I", data[offset : offset + 4])[0]
        compressed_data_start = offset + 4

        # Extract block 1 data and parse it
        block_1_data = data[offset : compressed_data_start + compressed_length]
        compressed_result = parse_block_1(block_1_data)
        if compressed_result:
            node.update(compressed_result)

        # Move offset to after compressed section
        offset = compressed_data_start + compressed_length

    # Types 0, 21, 32: Uncompressed sequences
    elif seq_type in [0, 21, 32]:
        seq_length = struct.unpack(">I", data[offset : offset + 4])[0]
        offset += 4

        # Extract sequence data and parse using appropriate parser
        seq_data_with_type = bytes([seq_type]) + data[offset - 4 : offset + seq_length]

        if seq_type == 0:
            result = parse_block_0(seq_data_with_type)
        elif seq_type == 21:
            result = parse_block_21(seq_data_with_type)
        elif seq_type == 32:
            result = parse_block_32(seq_data_with_type)

        if result:
            node["sequence"] = result
            node["sequence_length"] = seq_length

        offset += seq_length

    # Parse remaining blocks (typically type 0x1E with LZMA data)
    if offset < len(data):
        stream = BytesIO(data[offset:])
        nested = parse_blocks(stream)
        if nested:
            node["node_info"] = nested

    return node


def parse_block_14(data):
    """Type 14: Enzyme custom XML"""
    return decode_xml(data)


def parse_block_17(data):
    """Type 17: Alignable sequences XML"""
    return decode_xml(data)


def parse_block_18(data):
    """Type 18: Sequence trace (ZTR format)"""
    if len(data) < 10 or data[:8] != b"\xaeZTR\r\n\x1a\n":
        return f"Invalid ZTR ({len(data)} bytes)"

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

        # Parse specific chunk types
        if chunk_type == "BASE" and chunk_data[0] == 0:
            result["BASE"] = chunk_data[2:].decode("ascii", errors="ignore")
        elif chunk_type == "TEXT" and chunk_data[0] == 0:
            items = chunk_data[2:-2].split(b"\x00")
            text = {}
            for i in range(0, len(items) - 1, 2):
                key = items[i].decode("ascii", errors="ignore")
                val = items[i + 1].decode("ascii", errors="ignore")
                text[key] = val
            result["TEXT"] = text
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
            result["SMP4"] = samples
        elif chunk_type == "CLIP" and len(chunk_data) >= 9:
            result["CLIP"] = {
                "left": struct.unpack(">I", chunk_data[1:5])[0],
                "right": struct.unpack(">I", chunk_data[5:9])[0],
            }

        offset += data_len

    return result


def parse_block_21(data):
    """Type 21: Protein sequence (uncompressed ASCII)"""
    try:
        return data[1:].decode("utf-8", errors="strict")
    except:
        return data.hex()


def parse_block_28(data):
    """Type 28: Enzyme visualization XML"""
    return decode_xml(data)


def parse_block_29(data):
    """Type 29: History modifier (LZMA compressed XML)"""
    return decode_lzma(data)


def parse_block_30(data):
    """Type 30: History content (LZMA compressed, contains nested TLV blocks)"""
    try:
        decompressed = lzma.decompress(data)
        stream = BytesIO(decompressed)
        return parse_blocks(stream)
    except:
        return data.hex()


def parse_block_32(data):
    """Type 32: RNA sequence (uncompressed ASCII)"""
    try:
        return data[1:].decode("utf-8", errors="strict")
    except:
        return data.hex()


# =============================================================================
# GLOBAL SCHEME
# =============================================================================

# Define parsing scheme: block_type -> (length_override, parser_function)
SCHEME = {
    0: (None, parse_block_0),
    1: (None, parse_block_1),
    2: (None, None),  # Unknown - skip
    3: (None, parse_block_3),
    4: (None, None),  # Unknown - skip
    5: (None, parse_block_5),
    6: (None, parse_block_6),
    7: (None, parse_block_7),
    8: (None, parse_block_8),
    9: (None, None),  # File description (legacy) - skip
    10: (None, parse_block_10),
    11: (None, parse_block_11),
    12: (None, None),  # Unknown - skip
    13: (None, None),  # Enzyme info - skip
    14: (None, parse_block_14),
    15: (None, None),  # Unknown - skip
    16: (4, None),  # Legacy trace: always 4 bytes - skip
    17: (None, parse_block_17),
    18: (None, parse_block_18),
    19: (None, None),  # Uracil positions - skip
    20: (None, None),  # Custom colors XML - skip
    21: (None, parse_block_21),
    22: (None, None),  # Unknown - skip
    23: (None, None),  # Unknown - skip
    24: (None, None),  # Unknown - skip
    25: (None, None),  # Unknown - skip
    26: (None, None),  # Unknown - skip
    27: (None, None),  # Unknown - skip
    28: (None, parse_block_28),
    29: (None, parse_block_29),
    30: (None, parse_block_30),
    31: (None, None),  # Unknown - skip
    32: (None, parse_block_32),
}


# =============================================================================
# SECTION 3: MAIN PARSING FUNCTIONS
# =============================================================================


def parse_sgff(filepath):
    """Parse SGFF file into JSON-serializable dictionary"""
    with open(filepath, "rb") as f:
        # Validate header
        if f.read(1) != b"\t":
            raise ValueError("Wrong format for a SnapGene file!")

        length = struct.unpack(">I", f.read(4))[0]
        title = f.read(8).decode("ascii")

        if length != 14 or title != "SnapGene":
            raise ValueError("Wrong format for a SnapGene file!")

        # Parse cookie
        result = {
            "cookie": {
                "magic": title,
                "type_of_sequence": struct.unpack(">H", f.read(2))[0],
                "export_version": struct.unpack(">H", f.read(2))[0],
                "import_version": struct.unpack(">H", f.read(2))[0],
            }
        }

        # Parse all blocks using generic parse_blocks function
        blocks = parse_blocks(f)
        result.update(blocks)

        return result


def main(filepath, output_file=None):
    """Parse SnapGene file and save as JSON"""
    result = parse_sgff(filepath)

    if output_file:
        with open(output_file, "w") as f:
            json.dump(result, f, indent=2)
        print(f"Saved to {output_file}")
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <snapgene_file> [output.json]")
        sys.exit(1)

    output = sys.argv[2] if len(sys.argv) > 2 else None
    main(sys.argv[1], output)
