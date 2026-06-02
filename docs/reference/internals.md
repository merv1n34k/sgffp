# Internals

## SCHEME Dict

The `SCHEME` dictionary in `sgffp.parsers` maps block type IDs to parser functions. This is the central dispatch table used by `parse_blocks()`.

```python
from sgffp.parsers import SCHEME
```

| Block ID | Parser | Content |
|----------|--------|---------|
| 0 | `parse_sequence` | DNA sequence |
| 1 | `parse_compressed_dna` | 2-bit compressed DNA |
| 5 | `parse_xml` | Primers (XML) |
| 6 | `parse_xml` | Notes (XML) |
| 7 | `parse_lzma_xml` | History tree (LZMA XML) |
| 8 | `parse_xml` | Sequence properties (XML) |
| 10 | `parse_features` | Features (XML + qualifier extraction) |
| 11 | `parse_history_node` | History node (binary) |
| 14 | `parse_xml` | Custom enzyme sets (XML) |
| 16 | `parse_trace_container` | Trace container (flags + nested TLV) |
| 17 | `parse_xml` | Alignable sequences (XML) |
| 18 | `parse_ztr` | ZTR trace data (inside block 16) |
| 20 | `parse_xml` | Strand colors (XML) |
| 21 | `parse_sequence` | Protein sequence |
| 28 | `parse_xml` | Enzyme visibilities (XML) |
| 29 | `parse_lzma_xml` | History modifier (LZMA XML) |
| 30 | `parse_lzma_nested` | History node content (LZMA nested TLV) |
| 32 | `parse_sequence` | RNA sequence |
| 34 | `parse_lzma_json` | RNA structure predictions (LZMA JSON) |

Blocks not in SCHEME (2, 3, 13) are skipped — SnapGene regenerates them on import.

---

## Parser Functions

All parsers accept `data: bytes` and return a parsed dict (or `None` on failure).

### `parse_blocks(stream) → Dict[int, List[Any]]`

Read TLV blocks from a binary stream and dispatch each to its SCHEME parser. Returns the top-level blocks dict used by `SgffObject`.

### `parse_sequence(data) → Dict`

Parse uncompressed sequence (blocks 0, 21, 32). Returns:

```python
{
    "sequence": str,
    "length": int,
    "topology": "linear" | "circular",
    "strandedness": "single" | "double",
    "dam_methylated": bool,
    "dcm_methylated": bool,
    "ecoki_methylated": bool,
}
```

### `parse_compressed_dna(data) → Dict`

Parse a compressed-DNA block (top-level type 1 or embedded in block 11). The payload is a chain of self-describing sections — plain DNA (`0x01`), IUPAC ambiguity (`0x02`), all-N runs (`0x03`) — followed by optional lowercase range pairs. Returns:

```python
{
    "sequence": str,          # decoded sequence with case applied
    "length": int,            # uncompressed character count
    "writer_stamp": int,      # opaque byte preserved for round-trip
}
```

See [format spec § Block 1](/format-spec#block-1-compressed-dna-sequence) for the full layout.

### `parse_xml(data) → Dict | None`

Parse XML block via `xmltodict` with clean JSON keys (strips `@` prefixes, converts `#text` to `_text`).

### `parse_lzma_xml(data) → Dict | None`

Decompress LZMA, then parse XML.

### `parse_lzma_json(data) → Any | None`

Decompress LZMA, then parse JSON.

### `parse_lzma_nested(data) → Dict[int, List] | None`

Decompress LZMA, then parse as nested TLV blocks.

### `parse_features(data) → Dict | None`

Parse feature XML with qualifier extraction and strand mapping. Returns:

```python
{
    "features": [{"name": str, "type": str, "strand": str, "start": int, "end": int, ...}],
    "wrapper_extras": {"nextValidID": str, ...},
}
```

### `parse_ztr(data) → Dict | None`

Parse ZTR-format Sanger sequencing trace. Handles chunks: BASE, BPOS, CNF4, SMP4, SAMP, TEXT, CLIP, COMM. Supports raw (format 0) and zlib (format 2) compression.

### `parse_trace_container(data) → Dict`

Parse block 16: 4-byte flags header + nested TLV blocks.

### `parse_history_node(data) → Dict`

Parse block 11 binary: node_index, sequence_type, sequence data, nested TLV blocks.

### `octet_to_dna(raw_data, base_count) → bytes`

Convert 2-bit encoded bytes to ASCII DNA. Encoding: G=00, A=01, T=10, C=11.

---

## SgffModel

Base class for all block-backed models.

```python
from sgffp.models.base import SgffModel
```

| Member | Description |
|--------|-------------|
| `BLOCK_IDS` | Tuple of relevant block type IDs |
| `exists → bool` | `True` if any relevant blocks exist |

### Protected Helpers

| Method | Description |
|--------|-------------|
| `_get_block(block_id) → Any \| None` | Get first item from block |
| `_set_block(block_id, value)` | Set block value (replaces) |
| `_get_blocks(block_id) → List` | Get all items from block |
| `_set_blocks(block_id, values)` | Set all block values |
| `_remove_block(block_id) → bool` | Remove block entirely |

---

## SgffListModel[T]

Generic base for list-backed models. Extends `SgffModel`.

```python
from sgffp.models.base import SgffListModel
```

| Member | Description |
|--------|-------------|
| `items → List[T]` | Lazily loaded item list |
| `add(item)` | Append item and sync |
| `remove(idx) → bool` | Remove by index and sync |
| `clear()` | Remove all items and sync |
| `len(model)` | Item count |
| `model[idx]` | Item at index |
| `for item in model` | Iterate items |

### Abstract Methods (subclass must implement)

| Method | Description |
|--------|-------------|
| `_load() → List[T]` | Parse items from block storage |
| `_sync()` | Write items back to block storage |

### Subclasses

- `SgffFeatureList` (block 10)
- `SgffPrimerList` (block 5)
- `SgffAlignmentList` (block 17)
- `SgffTraceList` (block 16)
