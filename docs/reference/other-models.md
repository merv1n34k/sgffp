# Other Models

## SgffNotes

File-level metadata (block 6). Inherits from `SgffModel`.

```python
from sgffp import SgffNotes
```

**BLOCK_IDS:** `(6,)`

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `data` | `Dict` | Raw notes dict |
| `description` | `str` | File description (settable) |
| `created` | `str \| None` | Creation timestamp |
| `last_modified` | `str \| None` | Last modification timestamp |

### Methods

| Method | Description |
|--------|-------------|
| `get(key, default=None)` | Get value by key |
| `set(key, value)` | Set value and sync |
| `remove(key) → bool` | Remove key and sync |

### Example

```python
sgff.notes.description = "My plasmid"
sgff.notes.set("Type", "Synthetic")
print(sgff.notes.created)         # "2024-01-15.12:00:00"
sgff.notes.remove("CustomMapLabel")
```

---

## SgffProperties

Sequence properties (block 8). Inherits from `SgffModel`.

```python
from sgffp import SgffProperties
```

**BLOCK_IDS:** `(8,)`

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `data` | `Dict` | Raw properties dict |

### Methods

| Method | Description |
|--------|-------------|
| `get(key, default=None)` | Get value by key |
| `set(key, value)` | Set value and sync |

### Example

```python
props = sgff.properties
print(props.get("UpstreamModification"))  # "Unmodified"
props.set("UpstreamStickiness", "4")
```

---

## SgffAlignmentList

Alignable sequences (block 17). Inherits from `SgffListModel[SgffAlignment]`.

```python
from sgffp import SgffAlignmentList, SgffAlignment
```

**BLOCK_IDS:** `(17,)`

Provides standard list operations: `add()`, `remove()`, `clear()`, indexing, iteration.

### SgffAlignment

Single alignable sequence.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `name` | `str` | `""` | Sequence name |
| `sequence` | `str` | `""` | Sequence string |
| `extras` | `Dict` | `{}` | Unmodeled attributes |

Methods: `from_dict(data)`, `to_dict()`.

---

## SgffTraceList

Chromatogram traces (block 16 containers). Inherits from `SgffListModel[SgffTrace]`.

```python
from sgffp import SgffTraceList
```

**BLOCK_IDS:** `(16,)`

Loads traces from block 16 containers — each container wraps a block 18 ZTR trace with optional block 8 properties. Block 18 never appears at the top level.

Provides standard list operations: `add()`, `remove()`, `clear()`, indexing, iteration.

### SgffTrace

Single sequence trace (Sanger chromatogram).

```python
from sgffp import SgffTrace
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `bases` | `str` | `""` | Base calls |
| `positions` | `List[int]` | `[]` | Base-to-sample positions |
| `confidence` | `List[int]` | `[]` | Confidence scores per base |
| `samples` | `SgffTraceSamples \| None` | `None` | ACGT channel intensities |
| `clip` | `SgffTraceClip \| None` | `None` | Quality clip boundaries |
| `text` | `Dict[str, str]` | `{}` | Metadata key-value pairs |
| `comments` | `List[str]` | `[]` | Comment strings |

#### Properties

| Property | Type | Description |
|----------|------|-------------|
| `sequence` | `str` | Alias for `bases` |
| `length` | `int` | Number of bases |
| `sample_count` | `int` | Number of sample points |

#### Methods

| Method | Description |
|--------|-------------|
| `get_metadata(key, default="")` | Get metadata value |
| `get_confidence_at(index) → int \| None` | Confidence at base index |
| `get_position_at(index) → int \| None` | Sample position at base index |
| `from_dict(data)` | Create from parsed dict |
| `to_dict()` | Serialize to dict |

### SgffTraceClip

Quality clip boundaries.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `left` | `int` | `0` | Left clip position |
| `right` | `int` | `0` | Right clip position |

### SgffTraceSamples

Trace sample intensities for each channel.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `a` | `List[int]` | `[]` | Adenine channel |
| `c` | `List[int]` | `[]` | Cytosine channel |
| `g` | `List[int]` | `[]` | Guanine channel |
| `t` | `List[int]` | `[]` | Thymine channel |

| Property | Description |
|----------|-------------|
| `length` | Number of sample points |

### Example

```python
for trace in sgff.traces:
    print(f"{trace.length} bases, {trace.sample_count} samples")
    if trace.clip:
        print(f"  clip: {trace.clip.left}..{trace.clip.right}")
    if trace.samples:
        print(f"  A channel: {trace.samples.a[:5]}...")
```

---

## SgffAttachmentList

File attachments (block 23). Inherits from `SgffListModel[SgffAttachment]`.

```python
from sgffp import SgffAttachmentList, SgffAttachment
```

**BLOCK_IDS:** `(23,)`

Provides standard list operations: `add()`, `remove()`, `clear()`, indexing, iteration. Additional lookup methods: `get_by_name(name)`, `get_by_id(file_id)`.

### SgffAttachment

Single file attachment embedded in a SnapGene file.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `id` | `int` | `0` | File ID (auto-assigned on add) |
| `name` | `str` | `""` | Filename |
| `data` | `bytes` | `b""` | Raw file content |
| `size` | `int` | `0` | File size in bytes |
| `mtime` | `str` | `""` | Modification timestamp |
| `compressible` | `str` | `"0"` | Whether SnapGene may compress |
| `extras` | `Dict` | `{}` | Unmodeled manifest attributes |

Methods: `from_blocks(file_block, manifest_entry)`, `to_manifest_dict()`, `to_file_block()`.

### Example

```python
for att in sgff.attachments:
    print(f"[{att.id}] {att.name}: {att.size} bytes")

sgff.add_attachment("gel.jpg", open("gel.jpg", "rb").read())
att = sgff.attachments.get_by_name("gel.jpg")
```

---

## SgffTraceAlignment

Trace alignment data (block 27). Inherits from `SgffModel`. Contains BGZF-compressed BAM data — alignments of sequencing traces against the reference sequence.

```python
from sgffp import SgffTraceAlignment, SgffBamRecord, SgffBamReference
```

**BLOCK_IDS:** `(27,)`

| Property | Type | Description |
|----------|------|-------------|
| `header` | `str` | SAM header text |
| `references` | `List[SgffBamReference]` | Reference sequences |
| `records` | `List[SgffBamRecord]` | Alignment records |
| `record_count` | `int` | Number of alignment records |
| `reference_count` | `int` | Number of reference sequences |

### SgffBamReference

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `name` | `str` | `""` | Reference name |
| `length` | `int` | `0` | Reference length in bp |

### SgffBamRecord

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `read_name` | `str` | `""` | Read/trace name (matches block 17 Sequence/@ID) |
| `flag` | `int` | `0` | SAM flag bits |
| `ref_id` | `int` | `0` | Reference sequence index |
| `pos` | `int` | `0` | 0-based mapping position |
| `mapq` | `int` | `255` | Mapping quality |
| `cigar` | `str` | `""` | CIGAR string (e.g. `3S154M6S`) |
| `sequence` | `str` | `""` | Read sequence |
| `quality` | `List[int]` | `[]` | Per-base quality scores |

Properties: `is_unmapped`, `is_reverse`, `length`.

### Example

```python
if sgff.has_trace_alignment:
    ta = sgff.trace_alignment
    print(f"SAM header: {ta.header}")
    for ref in ta.references:
        print(f"  Reference: {ref.name} ({ref.length} bp)")
    for rec in ta.records:
        strand = "reverse" if rec.is_reverse else "forward"
        print(f"  {rec.read_name}: {strand}, CIGAR={rec.cigar}, {rec.length}bp")
```
