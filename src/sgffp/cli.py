#!/usr/bin/env python3
"""
Command-line interface for SGFF tools
"""

import sys
import json
import argparse
import struct

from .reader import SgffReader
from .writer import SgffWriter
from .internal import SgffObject


# Blocks that are unknown or not yet decoded
NEW_BLOCKS = [2, 3, 4, 9, 12, 15, 19, 20, 22, 23, 24, 25, 26, 27, 31]


def cmd_parse(args):
    """Parse SGFF file to JSON"""
    sgff = SgffReader.from_file(args.input)

    output = {
        "cookie": {
            "type_of_sequence": sgff.cookie.type_of_sequence,
            "export_version": sgff.cookie.export_version,
            "import_version": sgff.cookie.import_version,
        },
        "blocks": sgff.blocks,
    }

    if args.output:
        with open(args.output, "w") as f:
            json.dump(output, f, indent=2)
    else:
        print(json.dumps(output, indent=2))


def cmd_info(args):
    """Show file information"""
    sgff = SgffReader.from_file(args.input)

    print(f"SnapGene File: {args.input}")
    print(f"Export version: {sgff.cookie.export_version}")
    print(f"Import version: {sgff.cookie.import_version}")
    print(f"\nBlocks:")

    block_types = {}
    for key in sgff.blocks.keys():
        block_type = key.split(".")[0]
        block_types[block_type] = block_types.get(block_type, 0) + 1

    for block_type in sorted(block_types.keys(), key=int):
        count = block_types[block_type]
        print(f"  Type {block_type:>2}: {count} block(s)")


def cmd_filter(args):
    """Filter blocks and write new file"""
    sgff = SgffReader.from_file(args.input)

    # Parse keep list
    keep_types = [int(t.strip()) for t in args.keep.split(",")]

    # Filter blocks
    filtered = SgffObject(cookie=sgff.cookie)
    for key, value in sgff.blocks.items():
        block_type = int(key.split(".")[0])
        if block_type in keep_types:
            filtered.blocks[key] = value

    # Write output
    SgffWriter.to_file(filtered, args.output)
    print(f"Filtered file written to {args.output}")


def cmd_check(args):
    """Check for unknown/new block types"""

    # Read file and scan for blocks
    found_blocks = {}
    new_found = []

    with open(args.input, "rb") as f:
        # Skip header
        f.read(1 + 4 + 8)  # magic + length + title
        f.read(2 + 2 + 2)  # cookie

        # Read all blocks
        while True:
            type_byte = f.read(1)
            if not type_byte:
                break

            block_type = type_byte[0]
            block_length = struct.unpack(">I", f.read(4))[0]
            block_data = f.read(block_length)

            # Track all found blocks
            if block_type not in found_blocks:
                found_blocks[block_type] = []
            found_blocks[block_type].append(block_data)

            # Check if this is a new/unknown block
            if block_type in NEW_BLOCKS:
                if block_type not in new_found:
                    new_found.append(block_type)

    # Report findings
    for block_type in sorted(found_blocks.keys()):
        count = len(found_blocks[block_type])
        marker = "[NEW]" if block_type in NEW_BLOCKS else ""
        print(f"{block_type:>2}: {count:>2} {marker}")

    # Alert if new blocks found
    if new_found:
        print()
        if args.examine:
            for block_type in sorted(new_found):
                for block_data in found_blocks[block_type]:
                    print("NEW BLOCK!")
                    print(f"Type: {block_type}, Length: {len(block_data)}")
                    print(block_data.hex())
                    print()
        else:
            print("NEW BLOCK!")
            print(f"Types: {sorted(new_found)}")


def main():
    parser = argparse.ArgumentParser(description="SnapGene File Format tools")
    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Parse command
    parse_parser = subparsers.add_parser("parse", help="Parse SGFF to JSON")
    parse_parser.add_argument("input", help="Input SGFF file")
    parse_parser.add_argument(
        "-o", "--output", help="Output JSON file (default: stdout)"
    )

    # Info command
    info_parser = subparsers.add_parser("info", help="Show file information")
    info_parser.add_argument("input", help="Input SGFF file")

    # Check command
    check_parser = subparsers.add_parser("check", help="Check for unknown block types")
    check_parser.add_argument("input", help="Input SGFF file")
    check_parser.add_argument(
        "-e",
        "--examine",
        action="store_true",
        help="Dump raw content of new/unknown blocks",
    )

    # Filter command
    filter_parser = subparsers.add_parser("filter", help="Filter blocks")
    filter_parser.add_argument("input", help="Input SGFF file")
    filter_parser.add_argument(
        "-k", "--keep", required=True, help="Block types to keep (comma-separated)"
    )
    filter_parser.add_argument("-o", "--output", required=True, help="Output SGFF file")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    if args.command == "parse":
        cmd_parse(args)
    elif args.command == "info":
        cmd_info(args)
    elif args.command == "check":
        cmd_check(args)
    elif args.command == "filter":
        cmd_filter(args)


if __name__ == "__main__":
    main()
