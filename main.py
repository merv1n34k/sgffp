#!/usr/bin/env python3
"""
Parse SnapGeneFileFormat file to JSON
"""

import struct
import json
import lzma
import zlib
import sys


def parse_sgff(filepath):
    """Parse SnapGene file format to JSON structure"""

    with open(filepath, "rb") as f:

        def unpack(size, mode):
            return struct.unpack(">" + mode, f.read(size))[0]

        # Validate header
        if f.read(1) != b"\t":
            raise ValueError("Wrong format for a SnapGene file!")

        length = unpack(4, "I")
        title = f.read(8).decode("ascii")

        if length != 14 or title != "SnapGene":
            raise ValueError("Wrong format for a SnapGene file!")

        result = {
            "cookie": {
                "magic": title,
                "type_of_sequence": unpack(2, "H"),
                "export_version": unpack(2, "H"),
                "import_version": unpack(2, "H"),
            }
        }

        # Parse blocks
        while True:
            # file structure table
            # 0: sequence_dna:utf-8
            # 1: compressed sequence:?
            # 2:
            # 3: enzyme_library:mixed
            # 4:
            # 5: primers:xml
            # 6: notes:xml
            # 7: history tree:xml
            # 8: additional sequence properties:xml
            # 9: file Description:?
            # 10: features:xml
            # 11: history node:lzma compressed:sgff content
            # 12:
            # 13: enzyme_info:mixed
            # 14: enzyme_custom:xml
            # 15:
            # 16: sequence trace(legacy): 4 empty bytes
            # 17: alignable sequences:xml
            # 18: sequence trace:zrt-trace format
            # 19: uracil Positions:?
            # 20: custom colors:xml
            # 21: sequence_protein:utf-8
            # 22:
            # 23:
            # 24:
            # 25:
            # 26:
            # 27: unknown:binary
            # 28: enzyme_vizualisation:xml
            # 29
            # 30:
            # 31:
            # 32: sequence_rna:utf-8
            type_byte = f.read(1)
            if not type_byte:
                break

            block_type = type_byte[0]
            block_size = unpack(4, "I")

            # History tree - decompress LZMA
            if block_type == 0x07:
                try:
                    result[str(block_type)] = lzma.decompress(
                        f.read(block_size)
                    ).decode("utf-8", errors="ignore")
                except:
                    pass

            # History node - extract and decompress
            elif block_type == 0x0B:
                node_index = unpack(4, "I")
                node_data = f.read(block_size - 4)

                if "11" not in result.keys():
                    result["11"] = {}

                # Find XZ signature and decompress
                xz_pos = node_data.find(b"\xfd7zXZ\x00")
                if xz_pos != -1:
                    try:
                        decompressed = lzma.decompress(node_data[xz_pos:])
                        result["11"][f"{node_index}"] = parse_node(decompressed)
                    except:
                        pass

            # Other blocks
            else:
                result[str(block_type)] = decode_block(f.read(block_size))

    return result


# THIS DOES NOT PARSE ZRT FILE BUT SHOULD: see ZRT spec
# https://staden.sourceforge.net/ztr.html
def parse_ztr(data):
    """Parse ZTR chromatogram format"""

    # Check ZTR magic
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

        # Decompress if needed (format byte != 0)
        if chunk_data and chunk_data[0] == 2:  # zlib
            try:
                chunk_data = b"\x00" + zlib.decompress(chunk_data[5:])
            except:
                pass

        # Parse chunk based on type
        if chunk_type == "BASE" and chunk_data[0] == 0:
            result["BASE"] = chunk_data[2:].decode("ascii", errors="ignore")
        elif chunk_type == "TEXT" and chunk_data[0] == 0:
            items = chunk_data[2:-2].split(b"\x00")
            text = {}
            for i in range(0, len(items) - 1, 2):
                text[items[i].decode("ascii", errors="ignore")] = items[i + 1].decode(
                    "ascii", errors="ignore"
                )
            result["TEXT"] = text
        elif chunk_type == "SMP4":
            trace_len = len(data) // 8

            samples = {}
            for i, base in enumerate(["A", "C", "G", "T"]):
                start = i * trace_len * 2
                trace = [
                    struct.unpack(">H", data[start + j : start + j + 2])[0]
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


def parse_node(data):
    """Parse decompressed node content - either XML or TLV blocks"""

    # If starts with XML, return directly
    if data[:5] == b"<?xml":
        return data.decode("utf-8", errors="ignore")

    # Parse as TLV blocks
    result = {}
    offset = 0
    block_counter = {}  # Track count for each block type

    while offset + 5 <= len(data):
        block_type = data[offset]
        block_size = struct.unpack(">I", data[offset + 1 : offset + 5])[0]

        if block_size > len(data) - offset - 5:
            break

        # Track block occurrence
        if block_type not in block_counter:
            block_counter[block_type] = 0
        block_counter[block_type] += 1

        # Create two-byte key: counter + type
        key = f"{block_type}.{block_counter[block_type]}"

        # Special handling for different block types
        if block_type == 16:  # 0x10 - legacy type - skip
            offset += 5 + 4
            continue

        block_data = data[offset + 5 : offset + 5 + block_size]

        if block_type == 18:  # 0x12 - ZTR trace data
            result[key] = parse_ztr(block_data)

        else:
            # Try to decode as string
            try:
                result[key] = block_data.decode("utf-8", errors="strict")
            except:
                if b"<?xml" in block_data[:100] or (
                    block_data and block_data[0:1] == b"<"
                ):
                    result[key] = block_data.decode("utf-8", errors="ignore")
                else:
                    result[key] = block_data.hex()

        offset += 5 + block_size

    return result if result else data.hex()


def decode_block(data):
    """Decode block data to string or hex"""

    try:
        return data.decode("utf-8", errors="strict")
    except:
        if b"<?xml" in data[:100]:
            return data.decode("utf-8", errors="ignore")
        return data.hex()


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
