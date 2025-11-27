#!/usr/bin/env python3
"""
Analyze SnapGene file history blocks
"""

import struct
import sys
import lzma
import json


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

    seq_dict = {"node": [], "data": []}

    for block in blocks:
        if block["type"] == 0x0B:  # History nodes
            data = block["data"]
            node_index = struct.unpack(">I", data[0:4])[0]

            offset = 4  # After node index

            # Sequence type (0=dna, 1=compressed, 32=rna, 21=protein)
            seq_type = data[offset]
            offset += 1
            if seq_type == 29 or seq_type == 32:
                continue

            # Sequence length
            seq_length = struct.unpack(">I", data[offset : offset + 4])[0]
            offset += 4
            seq_data = data[offset - 5 : offset - 5 + seq_length]
            offset += seq_length
            seq_dict["node"].append(node_index)
            seq_dict["data"].append(seq_data)

            # switch to 1 if need to parse xz part
            parse_xz = 0

            # Parse type 0x30 block
            if data[offset] == 0x1E and parse_xz:
                offset += 1
                info_length = struct.unpack(">I", data[offset : offset + 4])[0]
                offset += 4
                print(f"length {info_length}")

                # XZ data follows
                print(f"XZ offset {offset}")
                if data[offset : offset + 6] == b"\xfd7zXZ\x00":
                    try:
                        pass
                        xz_data = lzma.decompress(data[offset:])
                        print(f"  XZ decompressed: {len(xz_data)} bytes")
                    except:
                        print(f"  XZ decompression failed")

    for idx, node in enumerate(seq_dict["node"]):
        data = seq_dict["data"][idx]

        def pad(val, n=4):
            diff = abs(n - len(str(val)))
            return f"{val}{' ' * diff}"

        def fmt(d, start, end, show_len=False):
            chunk = d[start:end]
            if not chunk:
                return ""
            result = chunk.hex()
            if show_len and len(chunk) == 4:
                result += f"({pad(int.from_bytes(chunk, 'big'), 5)})"
            return result

        seq_hex = " ".join(
            filter(
                None,
                [
                    fmt(data, 0, 1),  # 1b
                    fmt(data, 1, 5, True),  # 4b(length)
                    fmt(data, 5, 9, True),  # 4b(length)
                    fmt(data, 9, 13),  # 4b(flag?)
                    fmt(data, 13, 14),  # 1b
                    fmt(data, 14, 18),  # 4b
                    fmt(data, 18, 19),  # 1b
                    fmt(data, 19, 23, True),  # 4b(length)
                    # fmt(data, 23, 30),  # rest
                ],
            )
        )
        print(f"node:{pad(node, 3)}:seq:{seq_hex}")


def parse_history_nodes_to_json(filepath):
    """Parse history nodes and output as JSON with properly decoded DNA sequences"""

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

    blocks = parse_blocks(filepath)
    nodes = []

    for block in blocks:
        if block["type"] != 0x0B:  # History nodes only
            continue

        data = block["data"]
        node = {}

        # Parse node index (4 bytes)
        node["node_index"] = struct.unpack(">I", data[0:4])[0]
        offset = 4

        # Parse sequence type (1 byte)
        seq_type = data[offset]
        node["sequence_type"] = seq_type
        offset += 1

        # Type 29: Modifier only (no sequence data)
        if seq_type == 29:
            # Read LZMA block: type (0x1E) + length (4B) + LZMA data
            if offset < len(data) and data[offset] == 0x1E:
                offset += 1  # Skip type byte
                lzma_length = struct.unpack(">I", data[offset : offset + 4])[0]
                offset += 4

                try:
                    lzma_data = data[offset : offset + lzma_length]
                    decompressed = lzma.decompress(lzma_data)
                    node["node_info"] = decompressed.decode("utf-8", errors="ignore")
                except:
                    node["node_info"] = None

            nodes.append(node)
            continue

        # Type 1: Compressed DNA (2-bit GATC encoding)
        if seq_type == 1:
            # Read compression header (12 bytes total)
            compressed_length = struct.unpack(">I", data[offset : offset + 4])[0]
            offset += 4
            compressed_data_start = offset  # Mark for type 0x1E calculation

            uncompressed_length = struct.unpack(">I", data[offset : offset + 4])[0]
            offset += 4

            # TODO: INVESTIGATE - 4 bytes unknown flag/property
            # Possible meanings: topology, edit flags, quality markers
            offset += 4

            node["compressed_length"] = compressed_length
            node["uncompressed_length"] = uncompressed_length

            # TODO: INVESTIGATE - 10 bytes unknown header data
            # Structure: 1B + 4B + 1B + 4B (always: 01 + flag + 01 + length_repeat)
            # Likely: sequence properties, validation flags, or compression metadata
            offset += 10

            # Extract and decode compressed DNA sequence
            total_bytes_needed = (uncompressed_length * 2 + 7) // 8
            seq_data = data[offset : offset + total_bytes_needed]
            node["sequence"] = octet_to_dna(seq_data, uncompressed_length).decode(
                "ascii"
            )

            # Jump to type 0x1E block using compressed_length
            offset = compressed_data_start + compressed_length

        # Types 0, 21, 32: Uncompressed sequences (bytewise ASCII)
        elif seq_type in [0, 21, 32]:
            seq_length = struct.unpack(">I", data[offset : offset + 4])[0]
            offset += 4
            node["sequence_length"] = seq_length

            seq_data = data[offset : offset + seq_length]
            node["sequence"] = seq_data.decode("ascii", errors="ignore")
            offset += seq_length

        # Parse type 0x1E block: LZMA-compressed node metadata
        if offset < len(data) and data[offset] == 0x1E:
            offset += 1
            info_length = struct.unpack(">I", data[offset : offset + 4])[0]
            offset += 4

            try:
                lzma_data = data[offset : offset + info_length]
                decompressed = lzma.decompress(lzma_data)
                node["node_info"] = decompressed.decode("utf-8", errors="ignore")
                node["node_info_compressed_size"] = info_length
            except:
                node["node_info"] = None

        nodes.append(node)

    return json.dumps(nodes, indent=2)


def show_compressed_dna_padding(filepath):
    """Analyze unused bits at the end of compressed DNA sequences"""

    blocks = parse_blocks(filepath)

    for block in blocks:
        if block["type"] != 0x0B:
            continue

        data = block["data"]
        node_index = struct.unpack(">I", data[0:4])[0]
        offset = 4

        # Sequence type
        seq_type = data[offset]
        offset += 1

        # Only process compressed DNA (type 1)
        if seq_type != 1:
            continue

        # Parse header: 12 bytes in 3 groups of 4
        compressed_length = struct.unpack(">I", data[offset : offset + 4])[0]
        offset += 4
        uncompressed_length = struct.unpack(">I", data[offset : offset + 4])[0]
        offset += 4
        unknown_flag1 = struct.unpack(">I", data[offset : offset + 4])[0]
        offset += 4

        # Skip: 1B(01) + 4B(flag) + 1B(01) + 4B(length_repeat) = 10 bytes
        offset += 10

        # Calculate bits needed for sequence
        total_bases = uncompressed_length
        bits_needed = total_bases * 2
        bytes_needed = (bits_needed + 7) // 8  # Round up
        bits_in_last_byte = bits_needed % 8
        unused_bits_in_last_byte = (8 - bits_in_last_byte) % 8

        # Extract compressed DNA data
        seq_data = data[offset : offset + bytes_needed]

        # Get the last byte and show unused bits
        if len(seq_data) > 0:
            last_byte = seq_data[-1]
            last_byte_binary = format(last_byte, "08b")

            # Split into used and unused bits
            used_bits = (
                last_byte_binary[:bits_in_last_byte]
                if bits_in_last_byte > 0
                else last_byte_binary
            )
            unused_bits = (
                last_byte_binary[bits_in_last_byte:]
                if unused_bits_in_last_byte > 0
                else ""
            )

            print(
                f"node:{node_index:3d} | bases:{total_bases:5d} | bits:{bits_needed:5d} | "
                f"bytes:{bytes_needed:3d} | unused_bits:{unused_bits_in_last_byte} | "
                f"last_byte:[{used_bits}|{unused_bits}] = {last_byte_binary}"
            )


def main(filepath):
    """Main function - uncomment the analysis you want to run"""

    # list_all_blocks(filepath)
    parse_node_structure(filepath)
    # parse_history_nodes_to_json(filepath)
    # show_compressed_dna_padding(filepath)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <snapgene_file>")
        sys.exit(1)

    main(sys.argv[1])
