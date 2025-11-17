#!/usr/bin/env python3
"""
Parse SnapGene file to JSON
"""

import struct
import json
import lzma
import sys


def parse_snapgene(filepath):
    """Parse SnapGene file to JSON structure"""

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
            # 3:
            # 4:
            # 5: primers:xml
            # 6: notes:xml
            # 7: history tree:xml
            # 8: additional sequence properties:xml
            # 9: file Description:?
            # 10: features:xml
            # 11: history node:lzma compressed:1.gzip zrt-trace file + 2.xml
            # 12:
            # 13: legacy?:binary
            # 14: enzyme_custom:xml
            # 15:
            # 16: alignable sequence:binary
            # 17: alignable sequence:xml
            # 18: sequence trace:?
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
            block_data = f.read(block_size)

            # History tree - decompress LZMA
            if block_type == 0x07:
                try:
                    result[str(block_type)] = lzma.decompress(block_data).decode(
                        "utf-8", errors="ignore"
                    )
                except:
                    pass

            # History node - extract and decompress
            elif block_type == 0x0B:
                node_index = struct.unpack(">I", block_data[0:4])[0]

                if "11" not in result.keys():
                    result["11"] = {}

                # Find XZ signature and decompress
                xz_pos = block_data.find(b"\xfd7zXZ\x00")
                if xz_pos != -1:
                    try:
                        decompressed = lzma.decompress(block_data[xz_pos:])
                        result["11"][f"{node_index}"] = process_node(decompressed)
                    except:
                        pass

            # Other blocks
            else:
                result[str(block_type)] = decode_block(block_data)

    return result


# THIS DOES NOT PARSE ZRT FILE BUT SHOULD: see ZRT spec
# https://staden.sourceforge.net/ztr.html
def process_node(data):
    """Extract XML content from node"""

    # If has CLIP marker, extract XML after it
    if b"CLIP" in data:
        clip_pos = data.find(b"CLIP")
        # Find first '<' after CLIP
        for i in range(clip_pos, min(clip_pos + 100, len(data))):
            if data[i : i + 1] == b"<":
                return data[i:].decode("utf-8", errors="ignore")

    # If starts with XML, return as is
    if data[:5] == b"<?xml":
        return data.decode("utf-8", errors="ignore")

    # Try TLV parsing
    if len(data) > 5:
        tlv = parse_tlv(data)
        if tlv:
            return tlv

    return data.hex()


def parse_tlv(data):
    """Parse TLV structure"""

    result = {}
    offset = 0

    while offset + 5 <= len(data):
        block_type = data[offset]
        block_size = struct.unpack(">I", data[offset + 1 : offset + 5])[0]

        if block_size > len(data) - offset - 5 or block_size > 1000000:
            break

        block_data = data[offset + 5 : offset + 5 + block_size]
        result[str(block_type)] = decode_block(block_data)
        offset += 5 + block_size

    return result if result else None


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

    result = parse_snapgene(filepath)

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
