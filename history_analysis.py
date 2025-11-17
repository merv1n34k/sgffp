#!/usr/bin/env python3
"""
Analyze SnapGene file history blocks
"""

import struct
import sys
import lzma


def parse_blocks(filepath):
    """Parse TestView file and return all blocks"""
    blocks = []

    with open(filepath, "rb") as f:
        # Skip header
        f.seek(0x13)

        block_num = 0
        while True:
            pos = f.tell()
            type_byte = f.read(1)

            if not type_byte:
                break

            block_type = type_byte[0]
            block_size = struct.unpack(">I", f.read(4))[0]
            block_data = f.read(block_size)

            blocks.append(
                {
                    "num": block_num,
                    "type": block_type,
                    "size": block_size,
                    "data": block_data,
                    "offset": pos,
                }
            )

            block_num += 1

    return blocks


def list_all_blocks(filepath):
    """List all blocks in the file"""
    blocks = parse_blocks(filepath)

    for block in blocks:
        print(
            f"Block {block['num']:2d}: Type 0x{block['type']:02X}, {block['size']:6d} bytes at 0x{block['offset']:04X}"
        )


def decompress(filepath):
    """Search entire node for LZMA/XZ and decompress"""
    blocks = parse_blocks(filepath)

    for block in blocks:
        if block["type"] == 0x0B:
            data = block["data"]
            node_index = struct.unpack(">I", data[0:4])[0]

            print(f"Node {node_index}: {block['size']} bytes")

            # Search entire node for XZ signature
            found = False
            for offset in range(len(data) - 6):
                if data[offset : offset + 6] == b"\xfd7zXZ\x00":
                    print(f"  XZ found at offset {offset} (0x{offset:X})")

                    try:
                        decompressed = lzma.decompress(data[offset:])
                        print(f"  Decompressed: {len(decompressed)} bytes")

                        # Save decompressed
                        filename = f"node_{node_index}_decompressed.bin"
                        with open(filename, "wb") as f:
                            f.write(decompressed)
                        print(f"  Saved to: {filename}")

                        found = True
                        break
                    except Exception as e:
                        print(f"  Decompression failed at {offset}: {e}")
                        continue

            if not found:
                print(f"  No XZ signature found")

                # Save raw data for manual inspection
                filename = f"node_{node_index}_raw.bin"
                with open(filename, "wb") as f:
                    f.write(data)
                print(f"  Saved raw to: {filename}")

            print()


def main(filepath):
    """Main function - uncomment the analysis you want to run"""

    # list_all_blocks(filepath)
    decompress(filepath)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <snapgene_file>")
        sys.exit(1)

    main(sys.argv[1])
