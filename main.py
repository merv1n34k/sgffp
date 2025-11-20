#!/usr/bin/env python3
"""
Parse SnapGeneFileFormat (SGFF) file to JSON using recursive TLV parsing
"""

import struct
import json
import lzma
import zlib
import sys


def unpack(fileobject, size, mode):
    return struct.unpack(">" + mode, fileobject.read(size))[0]


def read_header(fileobject):
    """Read TLV block header: 1 byte type + 4 bytes length"""
    type_byte = fileobject.read(1)
    if not type_byte:
        return None, None

    block_type = type_byte[0]
    block_length = unpack(fileobject, 4, "I")

    return block_type, block_length


def parse_sgff(filepath):
    """Parse SGFF into JSON"""
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

        # file structure table
        # 0: sequence_dna:ascii
        # 1: compressed sequence:WHICH COMPRESSOR?!
        # 2:
        # 3: enzyme_library:mixed (enzyme sites + id?)
        # 4:
        # 5: primers:xml
        # 6: notes:xml
        # 7: history tree:xml
        # 8: additional sequence properties:xml
        # 9: file Description:?
        # 10: features:xml
        # 11: history_node:tlv container for ((0/1/32)+30)/29
        # 12:
        # 13: enzyme_info:mixed
        # 14: enzyme_custom:xml
        # 15:
        # 16: sequence trace(legacy): 4 empty bytes
        # 17: alignable sequences:xml
        # 18: sequence trace:zrt-trace format
        # 19: uracil_positions:?
        # 20: custom_colors:xml
        # 21: sequence_protein:utf-8
        # 22:
        # 23:
        # 24:
        # 25:
        # 26:
        # 27: unknown:binary
        # 28: enzyme_vizualisation:xml
        # 29: history_modifier:lzma
        # 30: history_content:lzma (content from file it was taken, except for sequence)
        # 31:
        # 32: sequence_rna:ascii

        # Define main parsing scheme
        scheme = {
            0: parse_sequence,  # DNA sequence
            1: None,  # Compressed sequence - skip
            2: None,  # Unknown binary - skip
            3: parse_enzyme_cutter,  # Enzyme cutters - skip
            4: None,  # Unknown binary - skip
            5: parse_xml,  # Primers XML
            6: parse_xml,  # Notes XML
            7: parse_lzma,  # History tree (LZMA compressed XML)
            8: parse_xml,  # Additional sequence properties XML
            9: None,  # File description - skip
            10: parse_xml,  # Features XML
            11: None,  # History node (nested blocks) - skip for now
            12: None,  # Unknown binary - skip
            13: None,  # Enzyme info - skip
            14: parse_xml,  # Enzyme custom XML
            15: None,  # Unknown binary - skip
            16: (4, None),  # Legacy trace - read 4 bytes, skip
            17: parse_xml,  # Alignable sequences XML
            18: parse_ztr,  # Sequence trace ZTR
            19: None,  # Uracil positions - skip
            20: None,  # parse_xml,  # Custom colors XML
            21: parse_sequence,  # Protein sequence
            22: None,  # Unknown binary - skip
            23: None,  # Unknown binary - skip
            24: None,  # Unknown binary - skip
            25: None,  # Unknown binary - skip
            26: None,  # Unknown binary - skip
            27: None,  # Unknown binary - skip
            28: parse_xml,  # Enzyme visualization XML
            29: None,  # parse_lzma,  # History modifier (LZMA)
            30: None,  # parse_lzma,  # History content (LZMA)
            31: None,  # Unknown binary - skip
            32: parse_sequence,  # RNA sequence
        }

        # Parse all blocks in file using scheme
        while True:
            block_type, block_length = read_header(f)
            if block_type is None:
                break

            # Get parser from scheme (default to skip if not defined)
            block_scheme = scheme.get(block_type)

            # Handle different scheme formats
            if block_scheme is None:
                # Block type not in scheme - skip it
                f.read(block_length)
                continue
            elif type(block_scheme) is tuple:
                # Tuple format: (length_override, parser)
                block_length = block_scheme[0]
                block_parser = block_scheme[1]
            else:
                # Direct parser function
                block_parser = block_scheme

            # Read block data
            data = f.read(block_length)

            # If parser is None, skip this block
            if block_parser is None:
                continue

            # Parse and store result
            result[str(block_type)] = block_parser(data)

        return result


# === SPECIFIC PARSERS ===


def parse_xml(data):
    """Parse XML as UTF-8 text"""
    try:
        return data.decode("utf-8", errors="ignore")
    except:
        return data.hex()


def parse_sequence(data):
    """Parse as UTF-8 text, fallback to hex"""
    try:
        return data[1:].decode("utf-8", errors="strict")
    except:
        return data.hex()


def parse_lzma(data):
    """Parse LZMA compressed data as text"""
    try:
        return lzma.decompress(data).decode("utf-8", errors="ignore")
    except:
        return data.hex()


def parse_ztr(data):
    """Parse ZTR chromatogram format"""
    if len(data) < 10 or data[:8] != b"\xaeZTR\r\n\x1a\n":
        return f"Invalid ZTR ({len(data)} bytes)"

    result = {}
    offset = 10

    while offset + 12 <= len(data):
        # Read chunk header
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


def parse_enzyme_cutter(data):
    """Parse enzyme cutter sequences (type 3)"""
    if len(data) < 5 or data[0] != 1:
        return data.hex()

    # Read CSV length
    csv_length = struct.unpack(">I", data[1:5])[0]
    if len(data) < 5 + csv_length:
        return data.hex()

    # Parse CSV sequences
    csv_data = data[5 : 5 + csv_length]
    sequences = csv_data.decode("ascii", errors="ignore").split(",")

    # Parse remaining 4-byte values
    values = []
    offset = 5 + csv_length
    while offset + 4 <= len(data):
        val = struct.unpack(">I", data[offset : offset + 4])[0]
        values.append(val)
        offset += 4

    return {"sequences": sequences, "values": values}


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
