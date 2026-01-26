#!/usr/bin/env python3
"""
Command-line interface for SGFF tools
"""

import sys
import json
import argparse
import struct
from importlib.metadata import version

from .reader import SgffReader
from .writer import SgffWriter
from .internal import SgffObject
from .parsers import SCHEME

KNOWN_BLOCKS = set(SCHEME.keys())


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

    seq_type_names = {1: "DNA", 2: "RNA", 3: "Protein"}
    seq_type_name = seq_type_names.get(sgff.cookie.type_of_sequence, "Unknown")

    print(f"File: {args.input}")
    print(f"Type: {seq_type_name} (export v{sgff.cookie.export_version}, import v{sgff.cookie.import_version})")

    # Sequence info
    seq = sgff.sequence
    if seq.length > 0:
        topo = seq.topology
        strand = "double-stranded" if seq.is_double_stranded else "single-stranded"
        print(f"Sequence: {seq.length} bp, {topo}, {strand}")

    # Features
    if sgff.has_features:
        print(f"Features: {len(sgff.features)}")

    # Primers
    if sgff.has_primers:
        print(f"Primers: {len(sgff.primers)}")

    # History
    if sgff.has_history:
        print(f"History: {len(sgff.history)} nodes")

    # Blocks summary
    print(f"Blocks: {', '.join(str(t) for t in sorted(sgff.types))}")


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

            if block_type not in KNOWN_BLOCKS and block_type not in unknown:
                unknown.append(block_type)

    if args.list:
        for block_type in sorted(found_blocks.keys()):
            count = len(found_blocks[block_type])
            marker = "[NEW]" if block_type not in KNOWN_BLOCKS else ""
            print(f"{block_type:>2}: {count:>2} {marker}")

    if unknown:
        if args.list:
            print()
        if args.dump:
            for block_type in sorted(unknown):
                for block_data in found_blocks[block_type]:
                    print(f"Block {block_type}: {len(block_data)} bytes")
                    print(block_data.hex())
                    print()
        else:
            print(f"Unknown blocks: {sorted(unknown)}")


def main():
    parser = argparse.ArgumentParser(description="SnapGene File Format tools")
    parser.add_argument(
        "-v", "--version", action="version", version=f"sgffp {version('sgffp')}"
    )
    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Parse
    p = subparsers.add_parser("parse", help="Parse SGFF to JSON")
    p.add_argument("input", help="Input SGFF file")
    p.add_argument("-o", "--output", help="Output JSON file")

    # Info
    p = subparsers.add_parser("info", help="Show file information")
    p.add_argument("input", help="Input SGFF file")

    # Check
    p = subparsers.add_parser(
        "check",
        help="Check for unknown block types",
        description="Silent by default. Use -l to list blocks or -d to dump unknown block data.",
    )
    p.add_argument("input", help="Input SGFF file")
    p.add_argument("-l", "--list", action="store_true", help="List all block types")
    p.add_argument("-d", "--dump", action="store_true", help="Dump unknown block data")

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
