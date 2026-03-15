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
