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
from .parsers import SCHEME

# This list will be shortened when new blocks will be decoded
NEW_BLOCKS = [2, 3, 4, 9, 12, 15, 19, 20, 22, 23, 24, 25, 26, 27, 31]


def cmd_parse(args):
    """Parse SGFF file to JSON"""
    sgff = SgffReader.from_file(args.input)

    # Convert int keys to strings for JSON
    blocks_json = {str(k): v for k, v in sgff.blocks.items()}

    # Handle mystery bytes (not JSON serializable)
    for key, items in blocks_json.items():
        for item in items:
            if isinstance(item, dict) and "mystery" in item:
                item["mystery"] = item["mystery"].hex()

    output = {
        "cookie": {
            "type_of_sequence": sgff.cookie.type_of_sequence,
            "export_version": sgff.cookie.export_version,
            "import_version": sgff.cookie.import_version,
        },
        "blocks": blocks_json,
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
    print(f"Sequence type: {sgff.cookie.type_of_sequence}")
    print(f"Export version: {sgff.cookie.export_version}")
    print(f"Import version: {sgff.cookie.import_version}")
    print(f"\nBlocks:")

    for block_type in sorted(sgff.types):
        count = len(sgff.blocks[block_type])
        print(f"  Type {block_type:>2}: {count} block(s)")


def cmd_filter(args):
    """Filter blocks and write new file"""
    sgff = SgffReader.from_file(args.input)

    keep_types = {int(t.strip()) for t in args.keep.split(",")}

    filtered = SgffObject(cookie=sgff.cookie)
    for block_type in sgff.types:
        if block_type in keep_types:
            filtered.blocks[block_type] = sgff.blocks[block_type]

    SgffWriter.to_file(filtered, args.output)
    print(f"Filtered file written to {args.output}")


def cmd_check(args):
    """Check for unknown/new block types"""
    found_blocks = {}
    unknown = []

    with open(args.input, "rb") as f:
        f.read(1 + 4 + 8 + 6)  # Skip header + cookie

        while True:
            type_byte = f.read(1)
            if not type_byte:
                break

            block_type = type_byte[0]
            block_length = struct.unpack(">I", f.read(4))[0]
            block_data = f.read(block_length)

            if block_type not in found_blocks:
                found_blocks[block_type] = []
            found_blocks[block_type].append(block_data)

            if (
                block_type not in SCHEME
                and block_type in NEW_BLOCKS
                and block_type not in unknown
            ):
                unknown.append(block_type)

    if args.list:
        for block_type in sorted(found_blocks.keys()):
            count = len(found_blocks[block_type])
            marker = (
                "[NEW]" if block_type not in SCHEME and block_type in NEW_BLOCKS else ""
            )
            print(f"{block_type:>2}: {count:>2} {marker}")

    if unknown:
        print()
        if args.examine:
            for block_type in sorted(unknown):
                for block_data in found_blocks[block_type]:
                    print(f"NEW BLOCK: Type {block_type}, Length {len(block_data)}")
                    print(block_data.hex())
                    print()
        else:
            print(f"Unknown block types: {sorted(unknown)}")


def main():
    parser = argparse.ArgumentParser(description="SnapGene File Format tools")
    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Parse
    p = subparsers.add_parser("parse", help="Parse SGFF to JSON")
    p.add_argument("input", help="Input SGFF file")
    p.add_argument("-o", "--output", help="Output JSON file")

    # Info
    p = subparsers.add_parser("info", help="Show file information")
    p.add_argument("input", help="Input SGFF file")

    # Check
    p = subparsers.add_parser("check", help="Check for unknown block types")
    p.add_argument("input", help="Input SGFF file")
    p.add_argument("-e", "--examine", action="store_true", help="Dump unknown blocks")
    p.add_argument(
        "-l", "--list", action="store_true", help="List all blocks types and count"
    )

    # Filter
    p = subparsers.add_parser("filter", help="Filter blocks")
    p.add_argument("input", help="Input SGFF file")
    p.add_argument(
        "-k", "--keep", required=True, help="Block types to keep (comma-separated)"
    )
    p.add_argument("-o", "--output", required=True, help="Output SGFF file")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    commands = {
        "parse": cmd_parse,
        "info": cmd_info,
        "check": cmd_check,
        "filter": cmd_filter,
    }
    commands[args.command](args)


if __name__ == "__main__":
    main()
