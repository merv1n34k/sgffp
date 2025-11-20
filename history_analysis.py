#!/usr/bin/env python3
"""
Analyze SnapGene file history blocks
"""

import struct
import sys
import lzma
import zlib


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


# possible structure for nodes can be
# <1b-type-18><4b-node-length><4b-node-id>
# <1b-type-0/1><4b-sequence-length><1b-padding><nb-sequence>
# <1b-type-30><4b-compressed-info-length><nb-compressed-info>


def parse_node_structure(filepath):
    """Parse history nodes with the discovered structure"""
    blocks = parse_blocks(filepath)

    for block in blocks:
        if block["type"] == 0x0B:  # History nodes
            data = block["data"]
            node_index = struct.unpack(">I", data[0:4])[0]

            print(f"node:{node_index} - length:{len(data)}")

            offset = 4  # After node index

            # Sequence type (0=dna, 1=compressed, 32=rna, 21=protein)
            seq_type = data[offset]
            offset += 1
            # print(f"type  : {seq_type}")
            if seq_type == 29:
                continue

            # Sequence length
            seq_length = struct.unpack(">I", data[offset : offset + 4])[0]
            offset += 4
            print(f"type:{seq_type} - length:{seq_length}")
            seq_type2 = data[offset]
            offset += 0
            seq_length2 = struct.unpack(">I", data[offset : offset + 4])[0]
            offset += 0
            print(f"type:{seq_type2} - length:{seq_length2}")

            # Extract sequence
            sequence_data = data[offset : offset + seq_length]
            print(
                f"seq:{' '.join(sequence_data[i : i + 4].hex() for i in range(0, 30, 4))}"
            )
            offset += seq_length + 4

            if seq_type == 0:
                print(
                    f"  Raw DNA: {sequence_data[:50].decode('ascii', errors='ignore')}..."
                )
            else:
                try:
                    decompressed = zlib.decompress(sequence_data)
                    print(f"  Decompressed DNA: {len(decompressed)} bytes")
                    sequence_data = decompressed
                except:
                    pass
                    # print(f"  Failed to decompress")

            # Parse type 0x30 block
            if data[offset] == 0x1E:
                offset += 1
                info_length = struct.unpack(">I", data[offset : offset + 4])[0]
                offset += 4
                # print(f"  Compressed info length: {info_length}")

                # XZ data follows
                # print(f"XZ offset {offset}")
                if data[offset : offset + 6] == b"\xfd7zXZ\x00":
                    try:
                        xz_data = lzma.decompress(data[offset:])
                        # print(f"  XZ decompressed: {len(xz_data)} bytes")
                    except:
                        pass
                        # print(f"  XZ decompression failed")


def main(filepath):
    """Main function - uncomment the analysis you want to run"""

    # list_all_blocks(filepath)
    # decompress(filepath)
    parse_node_structure(filepath)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <snapgene_file>")
        sys.exit(1)

    main(sys.argv[1])
