# Sequence & Features

## SgffSequence

Wrapper for DNA, RNA, and Protein sequence blocks.

```python
from sgffp import SgffSequence
```

**BLOCK_IDS:** `(0, 1, 21, 32)` — DNA, compressed DNA, protein, RNA

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `value` | `str` | Sequence string (settable) |
| `length` | `int` | Sequence length |
| `topology` | `str` | `"linear"` or `"circular"` (settable) |
| `strandedness` | `str` | `"single"` or `"double"` (settable) |
| `is_circular` | `bool` | `True` if circular topology |
| `is_double_stranded` | `bool` | `True` if double-stranded |
| `block_id` | `int \| None` | Active block type (0, 1, 21, 32) (settable) |
| `data` | `SgffSequenceData` | Underlying dataclass |

Setting `value`, `topology`, or `strandedness` automatically syncs back to the blocks dict.

### SgffSequenceData

Underlying dataclass for sequence data.

| Field | Type | Default |
|-------|------|---------|
| `value` | `str` | `""` |
| `topology` | `str` | `"linear"` |
| `strandedness` | `str` | `"single"` |
| `dam_methylated` | `bool` | `False` |
| `dcm_methylated` | `bool` | `False` |
| `ecoki_methylated` | `bool` | `False` |

Methods: `from_dict(data)`, `to_dict()`, property `length`.

---

## SgffFeatureList

List model for annotation features. Inherits from `SgffListModel[SgffFeature]`.

```python
from sgffp import SgffFeatureList
```

**BLOCK_IDS:** `(10,)`

### Methods

| Method | Description |
|--------|-------------|
| `find_by_name(name) → SgffFeature \| None` | Find first feature by name |
| `find_by_type(type_) → List[SgffFeature]` | Find all features of a type |
| `add(feature)` | Add a feature and sync |
| `remove(idx) → bool` | Remove by index and sync |
| `clear()` | Remove all features |
| `len(fl)` | Count of features |
| `fl[idx]` | Feature at index |
| `for f in fl` | Iterate features |

### SgffFeature

Single annotation feature.

```python
from sgffp import SgffFeature, SgffSegment
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `name` | `str` | | Feature name |
| `type` | `str` | | Feature type (e.g., `"CDS"`, `"promoter"`) |
| `strand` | `str` | `"+"` | `"+"`, `"-"`, `"."` (none), `"="` (both) |
| `segments` | `List[SgffSegment]` | `[]` | Location ranges |
| `qualifiers` | `Dict[str, Any]` | `{}` | Parsed qualifiers |
| `color` | `str \| None` | `None` | Color from first segment |
| `reading_frame` | `int \| None` | `None` | Reading frame offset (-3…-1, 1…3) |
| `extras` | `Dict` | `{}` | Unmodeled XML attributes |
| `raw_qualifiers` | `List \| None` | `None` | Raw qualifier dicts for lossless roundtrip |

#### Computed Properties

| Property | Type | Description |
|----------|------|-------------|
| `start` | `int` | Minimum start across all segments (0-based) |
| `end` | `int` | Maximum end across all segments (1-based) |
| `length` | `int` | `end - start` |

Methods: `from_dict(data)`, `to_dict()`.

### SgffSegment

Single feature location range.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `start` | `int` | | Start position (0-based) |
| `end` | `int` | | End position (1-based) |
| `color` | `str \| None` | `None` | Hex color (e.g., `"#00ff00"`) |
| `type` | `str` | `"standard"` | Segment type: `"standard"` or `"gap"` |
| `translated` | `bool` | `False` | Whether this segment is translated |
| `extras` | `Dict` | `{}` | Unmodeled XML attributes |

Methods: `from_dict(data)`, `to_dict()`.

---

## SgffPrimerList

List model for primers. Inherits from `SgffListModel[SgffPrimer]`.

```python
from sgffp import SgffPrimerList
```

**BLOCK_IDS:** `(5,)`

Provides the same list operations as `SgffFeatureList`: `add()`, `remove()`, `clear()`, indexing, iteration.

### SgffPrimer

Single primer definition.

```python
from sgffp import SgffPrimer
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `name` | `str` | `""` | Primer name |
| `sequence` | `str` | `""` | Primer sequence |
| `binding_sites` | `List[SgffBindingSite]` | `[]` | Parsed binding site data |
| `extras` | `Dict` | `{}` | Unmodeled XML attributes |

#### Computed Properties (backward compatible)

| Property | Type | Description |
|----------|------|-------------|
| `bind_position` | `int \| None` | Start position of first non-simplified binding site (0-based) |
| `bind_strand` | `str` | Strand of first non-simplified binding site (`"+"` or `"-"`) |
| `melting_temperature` | `float \| None` | Melting temperature of first non-simplified binding site |

#### Methods

| Method | Description |
|--------|-------------|
| `clear_binding_sites()` | Remove all binding site data so SnapGene recalculates it |
| `from_dict(data)` | Create from parsed dict |
| `to_dict()` | Serialize to dict |

### SgffBindingSite

A single primer binding site on the target sequence. SnapGene stores both detailed and simplified versions of each site.

```python
from sgffp import SgffBindingSite
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `start` | `int` | `0` | Start position (0-based) |
| `end` | `int` | `0` | End position (1-based) |
| `bound_strand` | `str` | `"+"` | `"+"` (forward/top) or `"-"` (reverse/bottom) |
| `annealed_bases` | `str` | `""` | Annealed base sequence |
| `melting_temperature` | `float \| None` | `None` | Melting temperature in °C |
| `simplified` | `bool` | `False` | Whether this is a simplified binding site |
| `extras` | `Dict` | `{}` | Unmodeled data (e.g., `Component` elements) |

Methods: `from_dict(data)`, `to_dict()`.
