# SgffObject & I/O

## SgffObject

Container for all SnapGene file data.

```python
from sgffp import SgffObject
```

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `cookie` | `Cookie` | File header metadata |
| `blocks` | `Dict[int, List[Any]]` | Parsed block data keyed by type ID |

### Constructor

```python
SgffObject(cookie=Cookie(), blocks={})
```

### Class Methods

#### `new(sequence, *, topology, strandedness, sequence_type) → SgffObject`

Create a new SnapGene file with sensible defaults.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sequence` | `str` | `""` | DNA/RNA/protein sequence |
| `topology` | `str` | `"linear"` | `"linear"` or `"circular"` |
| `strandedness` | `str` | `"double"` | `"single"` or `"double"` |
| `sequence_type` | `str` | `"dna"` | `"dna"`, `"rna"`, or `"protein"` |

### Model Properties

All lazily loaded and cached. Call `invalidate()` to reset.

| Property | Type | Description |
|----------|------|-------------|
| `sequence` | `SgffSequence` | DNA/RNA/Protein sequence |
| `features` | `SgffFeatureList` | Annotation features |
| `primers` | `SgffPrimerList` | Primers |
| `notes` | `SgffNotes` | File notes |
| `properties` | `SgffProperties` | Sequence properties |
| `history` | `SgffHistory` | Edit history |
| `alignments` | `SgffAlignmentList` | Alignable sequences |
| `traces` | `SgffTraceList` | Chromatogram traces |
| `ops` | `SgffOps` | Operations API |

### Existence Checks

| Property | Checks block |
|----------|-------------|
| `has_notes` | 6 |
| `has_properties` | 8 |
| `has_features` | 10 |
| `has_primers` | 5 |
| `has_history` | 7, 11, 29, 30 |
| `has_alignments` | 17 |
| `has_traces` | 16 |

### Block Access Methods

| Method | Description |
|--------|-------------|
| `types` | List of block type IDs present |
| `type(block_id) → BlockList` | Get all items for a block type |
| `block(block_id) → Any \| None` | Get first item of a block type |
| `set(block_id, value)` | Append value to block type |
| `bset(block_id, value)` | Replace entire block type |
| `remove(block_id, idx=0) → bool` | Remove single item by index |
| `bremove(block_id) → bool` | Remove entire block type |
| `invalidate()` | Clear all cached model instances |

### Builder Methods

All return `self` for chaining.

#### `add_feature(name, type, start, end, strand="+", **kwargs) → SgffObject`

Add an annotation feature.

#### `add_primer(name, sequence, bind_position=None, bind_strand="+") → SgffObject`

Add a primer.

#### `set_sequence(new_sequence, operation="replace", record_history=True) → SgffObject`

Set a new sequence. If `record_history=True` and history exists, records the change in the history tree.

---

## Cookie

File header metadata.

```python
from sgffp import Cookie
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `type_of_sequence` | `int` | `1` | 1=DNA, 2=Protein, 7=RNA |
| `export_version` | `int` | `15` | SnapGene export version |
| `import_version` | `int` | `19` | SnapGene import version |

---

## BlockList

Wrapper for accessing items of a single block type. Returned by `sgff.type(block_id)`.

| Member | Description |
|--------|-------------|
| `type` | Block type ID |
| `first` | First item (or `None`) |
| `last` | Last item (or `None`) |
| `get(idx=0)` | Item at index (or `None`) |
| `len(bl)` | Number of items |
| `bl[idx]` | Item at index |
| `for item in bl` | Iterate items |

---

## SgffReader

Read and parse SnapGene files into `SgffObject`.

```python
from sgffp import SgffReader
```

| Method | Description |
|--------|-------------|
| `SgffReader(source)` | Create reader from path, `Path`, or `BinaryIO` |
| `reader.read() → SgffObject` | Parse and return |
| `SgffReader.from_file(path) → SgffObject` | Read from file path |
| `SgffReader.from_bytes(data) → SgffObject` | Read from `bytes` |

---

## SgffWriter

Write `SgffObject` to SnapGene file format.

```python
from sgffp import SgffWriter
```

| Method | Description |
|--------|-------------|
| `SgffWriter(target)` | Create writer to path, `Path`, or `BinaryIO` |
| `writer.write(sgff)` | Write `SgffObject` to target |
| `SgffWriter.to_file(sgff, path)` | Write to file path |
| `SgffWriter.to_bytes(sgff) → bytes` | Write to bytes |

Blocks are written in sorted order by type ID.
