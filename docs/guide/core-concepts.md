# Core Concepts

## SgffObject

`SgffObject` is the central container for all SnapGene file data. It holds:

- **`cookie`** — file header metadata (`Cookie` dataclass)
- **`blocks`** — a `Dict[int, List[Any]]` mapping block type IDs to parsed data

```python
from sgffp import SgffReader

sgff = SgffReader.from_file("plasmid.dna")

# Header info
sgff.cookie.type_of_sequence  # 1 (DNA), 2 (Protein), 7 (RNA)
sgff.cookie.export_version    # e.g. 15
sgff.cookie.import_version    # e.g. 19

# Raw block access
sgff.types           # list of block IDs present, e.g. [0, 5, 6, 10]
sgff.block(0)        # first item of block type 0 (sequence dict)
sgff.type(10)        # BlockList for block type 10 (features)
```

## TLV Blocks

SnapGene files are a sequence of **Type-Length-Value** blocks after a 19-byte header. Each block has a 1-byte type ID and a 4-byte big-endian length.

The `blocks` dict groups parsed data by type ID. Most block types appear once, but some (like block 11 for history nodes or block 16 for traces) appear multiple times.

Key block types:

| ID | Content | Model |
|----|---------|-------|
| 0 | DNA sequence | `SgffSequence` |
| 1 | Compressed DNA | `SgffSequence` |
| 5 | Primers (XML) | `SgffPrimerList` |
| 6 | Notes (XML) | `SgffNotes` |
| 7 | History tree (LZMA XML) | `SgffHistory` |
| 8 | Properties (XML) | `SgffProperties` |
| 10 | Features (XML) | `SgffFeatureList` |
| 11 | History nodes (binary) | `SgffHistory` |
| 16 | Trace container | `SgffTraceList` |
| 17 | Alignable sequences (XML) | `SgffAlignmentList` |
| 21 | Protein sequence | `SgffSequence` |
| 32 | RNA sequence | `SgffSequence` |

See the [Format Spec](/format-spec) for the full block reference.

## Model Accessors

`SgffObject` provides cached property accessors for typed models:

```python
sgff.sequence    # SgffSequence
sgff.features    # SgffFeatureList
sgff.primers     # SgffPrimerList
sgff.notes       # SgffNotes
sgff.properties  # SgffProperties
sgff.history     # SgffHistory
sgff.alignments  # SgffAlignmentList
sgff.traces      # SgffTraceList
sgff.ops         # SgffOps
```

Models are **lazily loaded** — they parse block data on first access and cache the result. Each model writes changes back to the underlying `blocks` dict automatically when you modify it.

### Existence Checks

Use `has_*` properties to check if data exists without triggering a parse:

```python
if sgff.has_features:
    print(f"{len(sgff.features)} features")

if sgff.has_history:
    print(f"{len(sgff.history)} history nodes")
```

Available checks: `has_notes`, `has_properties`, `has_features`, `has_primers`, `has_history`, `has_alignments`, `has_traces`.

### Cache Invalidation

If you modify `blocks` directly (bypassing model accessors), call `invalidate()` to clear cached models:

```python
sgff.blocks[10] = [new_features_dict]
sgff.invalidate()  # clears all cached model instances
```

## Block-Level Access

For low-level manipulation, `SgffObject` provides direct block methods:

```python
# Read
sgff.block(6)         # first item of block 6 (notes), or None
sgff.type(11)         # BlockList for block 11 (history nodes)
sgff.type(11).first   # first history node
sgff.type(11).last    # last history node
len(sgff.type(11))    # count of history nodes

# Write
sgff.set(6, notes_dict)       # append to block 6
sgff.bset(6, notes_dict)      # replace entire block 6
sgff.remove(6)                # remove first item from block 6
sgff.bremove(6)               # remove entire block 6
```

## The `extras` Dict Pattern

Models use an `extras` dict to preserve unmodeled XML attributes for lossless round-trips. When a model serializes with `to_dict()`, it starts from `extras` and overlays the modeled fields:

```python
feature = sgff.features[0]
feature.name       # modeled field
feature.extras     # dict of unmodeled XML attributes

# Round-trip: extras flow through transparently
feature.to_dict()  # merges extras + modeled fields
```

This means attributes that sgffp doesn't explicitly model are preserved when reading and writing back.

## Builder Methods

`SgffObject` provides builder methods that return `self` for chaining:

```python
sgff = (
    SgffObject.new("ATGCATGC", topology="circular")
    .add_feature("ori", "rep_origin", 0, 4)
    .add_primer("seq_fwd", "ATGC", bind_position=0)
    .set_sequence("ATGCATGCATGC", operation="replace")
)
```

- **`add_feature(name, type, start, end, strand="+", **kwargs)`** — adds an annotation feature
- **`add_primer(name, sequence, bind_position=None, bind_strand="+")`** — adds a primer
- **`set_sequence(new_sequence, operation="replace", record_history=True)`** — sets the sequence and optionally records a history entry
