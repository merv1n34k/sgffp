# Getting Started

## Installation

```bash
pip install sgffp
```

Or with [uv](https://docs.astral.sh/uv/):

```bash
uv add sgffp
```

For development:

```bash
git clone https://github.com/merv1n34k/sgffp.git
cd sgffp
uv sync --dev
```

Requires **Python 3.12+**.

## Read a File

```python
from sgffp import SgffReader

sgff = SgffReader.from_file("plasmid.dna")

# Sequence
print(sgff.sequence.value[:50])   # first 50 bases
print(sgff.sequence.length)       # e.g. 5420
print(sgff.sequence.topology)     # "circular" or "linear"

# Features
for f in sgff.features:
    print(f"{f.name} ({f.type}) {f.start}-{f.end}")

# Primers
for p in sgff.primers:
    print(f"{p.name}: {p.sequence}")
```

## Export to JSON

```python
import json
from sgffp import SgffReader

sgff = SgffReader.from_file("plasmid.dna")

blocks_json = {str(k): v for k, v in sgff.blocks.items()}
output = {
    "cookie": {
        "type_of_sequence": sgff.cookie.type_of_sequence,
        "export_version": sgff.cookie.export_version,
        "import_version": sgff.cookie.import_version,
    },
    "blocks": blocks_json,
}

print(json.dumps(output, indent=2))
```

Or use the CLI:

```bash
sff parse plasmid.dna -o plasmid.json
```

## Write Back

```python
from sgffp import SgffReader, SgffWriter

sgff = SgffReader.from_file("plasmid.dna")

# Modify something
sgff.notes.description = "Updated plasmid"
sgff.features.find_by_name("GFP").name = "eGFP"

# Write to a new file
SgffWriter.to_file(sgff, "updated.dna")
```

## Create from Scratch

Use the builder pattern to create a new SnapGene file without reading an existing one:

```python
from sgffp import SgffObject, SgffWriter

sgff = (
    SgffObject.new("ATGCATGCATGCATGC", topology="circular")
    .add_feature("GFP", "CDS", 0, 8)
    .add_feature("AmpR", "CDS", 8, 16, strand="-")
    .add_primer("fwd", "ATGCATGC", bind_position=0)
)

sgff.notes.description = "My plasmid"

SgffWriter.to_file(sgff, "new_plasmid.dna")
```

`SgffObject.new()` accepts:

| Parameter | Default | Options |
|-----------|---------|---------|
| `sequence` | `""` | DNA/RNA/protein string |
| `topology` | `"linear"` | `"linear"`, `"circular"` |
| `strandedness` | `"double"` | `"single"`, `"double"` |
| `sequence_type` | `"dna"` | `"dna"`, `"rna"`, `"protein"` |

All builder methods (`add_feature`, `add_primer`, `set_sequence`) return `self` for chaining.
