# Extras Reference

Every model uses an `extras` dict to preserve unmodeled XML attributes through roundtrip. This page documents the known extras keys, their types, and meanings.

All extras values are strings unless noted otherwise.

---

## SgffFeature

Accessed via `feature.extras`.

| Key | Type | Description |
|-----|------|-------------|
| `recentID` | `str` | Internal SnapGene ID for the feature |
| `translationMW` | `str` | Molecular weight of translated protein (Da) |
| `hitsStopCodon` | `str` | `"1"` if translation hits a stop codon |
| `consecutiveTranslationNumbering` | `str` | `"1"` to use consecutive numbering |

> `readingFrame` was formerly in extras; it is now the typed field `reading_frame`.

## SgffSegment

Accessed via `segment.extras`.

| Key | Type | Description |
|-----|------|-------------|
| `name` | `str` | Optional segment name |

> `type` and `translated` were formerly in extras; they are now typed fields.

## SgffPrimer

Accessed via `primer.extras`.

| Key | Type | Description |
|-----|------|-------------|
| `recentID` | `str` | Internal SnapGene ID |
| `color` | `str` | Hex color (e.g., `"#FF0000"`) |
| `description` | `str` | User-provided description |
| `dateAdded` | `str` | Creation date |

## SgffBindingSite

Accessed via `binding_site.extras`.

| Key | Type | Description |
|-----|------|-------------|
| `Component` | `dict` | Nested structure with hybridization details (`hybridizedRange`, `bases`) |

## SgffHistoryTreeNode

Accessed via `node.extras`.

| Key | Type | Description |
|-----|------|-------------|
| `hideHistory` | `str` | `"1"` to hide node from history tree view |
| `customMapLabel` | `str` | Custom label for map view |
| `useCustomMapLabel` | `str` | `"1"` to enable custom map label |
| `upstreamStickiness` | `str` | Upstream sticky end sequence (e.g., `"AATT"`) |
| `downstreamStickiness` | `str` | Downstream sticky end sequence (e.g., `"TTAA"`) |
| `RegeneratedSite` | `dict` | Enzyme info (`{"name": "BamHI"}`) |

## SgffInputSummary

Accessed via `summary.extras`.

| Key | Type | Description |
|-----|------|-------------|
| `strainName` | `str` | Bacterial strain (e.g., `"DH5alpha"`) |
| `methylationChanges` | `str` | Methylation pattern (e.g., `"dam+"`) |
| `dam2` | `str` | `"1"` for Dam methylation flag |

## SgffAlignment

Accessed via `alignment.extras`.

| Key | Type | Description |
|-----|------|-------------|
| `ID` | `str` | Unique identifier |
| `sortOrder` | `str` | Display ordering index |
| `trimmedRange` | `str` | Quality trimming range (e.g., `"5-100"`) |
| `isTrace` | `str` | `"0"` or `"1"` |

---

## Wrapper-level Extras

List models store wrapper-level extras in `_wrapper_extras` (not part of the public API but preserved for lossless roundtrip).

### SgffPrimerList

| Key | Type | Description |
|-----|------|-------------|
| `HybridizationParams` | `dict` | Hybridization parameters (`minContinuous`, etc.) |
| `nextValidID` | `str` | Next available primer ID counter |

### SgffFeatureList

| Key | Type | Description |
|-----|------|-------------|
| `nextValidID` | `str` | Next available feature ID counter |
| `recycledIDs` | `str` | Comma-separated list of deleted feature IDs |

### SgffAlignmentList

| Key | Type | Description |
|-----|------|-------------|
| `trimStringency` | `str` | Alignment quality threshold (e.g., `"0.05"`) |
