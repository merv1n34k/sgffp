# SnapGene .dna File Format Specification

This document describes the binary format of SnapGene `.dna` files based on reverse engineering and analysis of files produced by SnapGene versions 5.x–7.x.

---

## File Header

Every `.dna` file starts with a fixed 19-byte header:

| Offset | Size | Value | Description |
|--------|------|-------|-------------|
| 0 | 1 | `0x09` (`\t`) | Magic byte |
| 1 | 4 | `0x00000E` (14) | Header length (big-endian uint32) |
| 5 | 8 | `SnapGene` | ASCII title |
| 13 | 2 | varies | `type_of_sequence` (big-endian uint16) |
| 15 | 2 | varies | `export_version` (big-endian uint16) |
| 17 | 2 | varies | `import_version` (big-endian uint16) |

**type_of_sequence values:**

| Value | Meaning |
|-------|---------|
| 1 | DNA |
| 2 | RNA |
| 3 | Protein |

---

## TLV Block Format

After the header, the file contains a sequence of TLV (Type-Length-Value) blocks:

| Field | Size | Description |
|-------|------|-------------|
| type | 1 byte | Block type ID (unsigned) |
| length | 4 bytes | Data length (big-endian uint32) |
| data | `length` bytes | Block payload |

Blocks appear in any order. Some types can appear multiple times (e.g., block 11 for each history node, block 16 for each trace). Unknown block types should be skipped.

---

## Block Reference

### Block 0 — DNA Sequence

Uncompressed DNA sequence with property flags.

| Offset | Size | Description |
|--------|------|-------------|
| 0 | 1 | Property flags byte |
| 1 | N | ASCII sequence data |

**Property flags (bitmask):**

| Bit | Mask | Meaning |
|-----|------|---------|
| 0 | `0x01` | Circular topology (0=linear) |
| 1 | `0x02` | Double-stranded (0=single) |
| 2 | `0x04` | Dam methylated |
| 3 | `0x08` | Dcm methylated |
| 4 | `0x10` | EcoKI methylated |

### Block 1 — Compressed DNA Sequence

2-bit encoded DNA used in history nodes for compact storage.

| Offset | Size | Description |
|--------|------|-------------|
| 0 | 4 | `compressed_length` (big-endian uint32) — total bytes of remaining data |
| 4 | 4 | `uncompressed_length` (big-endian uint32) — number of bases |
| 8 | 14 | Metadata header (see below) |
| 22 | N | 2-bit GATC-encoded sequence data |

**14-byte metadata header:**

| Offset | Size | Description |
|--------|------|-------------|
| 0 | 1 | `format_version` — typically 30 (0x1e) |
| 1 | 3 | Reserved (always 0x000000) |
| 4 | 1 | `strandedness_flag` — 1 = double-stranded |
| 5 | 3 | Reserved (always 0x000000) |
| 8 | 2 | `property_flags` (big-endian uint16) — 1 = default, 257 = extended |
| 10 | 2 | Reserved (always 0x0000) |
| 12 | 2 | `header_seq_length` (big-endian uint16) — matches `uncompressed_length` |

**2-bit encoding** (2 bits per base, 4 bases per byte, MSB first):

| Bits | Base |
|------|------|
| 00 | G |
| 01 | A |
| 10 | T |
| 11 | C |

Total bytes = ceil(uncompressed_length * 2 / 8).

### Block 5 — Primers (XML)

XML block parsed by `xmltodict`. Top-level element: `<Primers>`.

```xml
<Primers nextValidID="3">
  <HybridizationParams
    minContinuousMatchLen="10"
    allowMismatch="1"
    minMeltingTemperature="40"
    showAdditionalFivePrimeMatches="1"
    minimumFivePrimeAnnealing="15"
  />
  <Primer
    name="Forward"
    sequence="ATGCATGCATGC"
    bindingSite="100"
    strand="+"
  />
</Primers>
```

### Block 6 — Notes (XML)

File-level metadata. Top-level element: `<Notes>`.

```xml
<Notes>
  <UUID>550e8400-e29b-41d4-a716-446655440000</UUID>
  <Type>Synthetic</Type>
  <Description>My plasmid</Description>
  <ConfirmedExperimentally>0</ConfirmedExperimentally>
  <CustomMapLabel>plasmid.dna</CustomMapLabel>
  <UseCustomMapLabel>1</UseCustomMapLabel>
  <SequenceClass>UNA</SequenceClass>
  <Created UTC="2024-01-15.12:00:00">2024-01-15.12:00:00</Created>
  <LastModified UTC="2024-06-01.15:30:00">2024-06-01.15:30:00</LastModified>
</Notes>
```

### Block 7 — History Tree (LZMA-compressed XML)

The entire block is LZMA-compressed. After decompression, it contains XML with a recursive `<Node>` tree describing the cloning operation history. Root node represents the current file state; children are previous states (tree grows backward in time).

```xml
<HistoryTree>
  <Node ID="2" name="Final.dna" type="DNA" seqLen="5000"
        circular="1" strandedness="double" operation="insertFragment"
        upstreamModification="Unmodified" downstreamModification="Unmodified">
    <InputSummary manipulation="insert" val1="100" val2="4900"
                  name1="EcoRI" siteCount1="1" name2="BamHI" siteCount2="1"/>
    <Node ID="0" name="Vector.dna" type="DNA" seqLen="4000"
          circular="1" strandedness="double" operation="invalid"
          resurrectable="1" .../>
    <Node ID="1" name="Insert.dna" type="DNA" seqLen="1000"
          circular="0" strandedness="double" operation="invalid"
          resurrectable="1" .../>
  </Node>
</HistoryTree>
```

**Node fields:**

| Field | Type | Description |
|-------|------|-------------|
| `ID` | string | Unique node identifier (integer as string) |
| `name` | string | Sequence filename |
| `type` | string | `"DNA"`, `"RNA"`, or `"Protein"` |
| `seqLen` | string | Sequence length (integer as string) |
| `circular` | string | `"0"` (linear) or `"1"` (circular) |
| `strandedness` | string | `"single"` or `"double"` |
| `operation` | string | Operation that created this state (see below) |
| `upstreamModification` | string | e.g., `"Unmodified"`, `"FivePrimePhosphorylated"` |
| `downstreamModification` | string | e.g., `"Unmodified"` |
| `resurrectable` | string | `"1"` if the user can restore this state |

**Operation types:**

| Operation | Description |
|-----------|-------------|
| `invalid` | Original/imported file (leaf node) |
| `makeDna` | Created new DNA sequence |
| `makeRna` | Created new RNA sequence |
| `makeProtein` | Created new protein sequence |
| `amplifyFragment` | PCR amplification |
| `insertFragment` | Restriction/ligation cloning |
| `replace` | Sequence edit |
| `digest` | Restriction digest |
| `ligate` | Ligation |
| `gatewayLR` | Gateway LR reaction |
| `gatewayBP` | Gateway BP reaction |
| `gibsonAssembly` | Gibson assembly |
| `goldenGateAssembly` | Golden Gate assembly |
| `restrictionClone` | Restriction cloning |
| `taClone` | TA cloning |
| `topoClone` | TOPO cloning |
| `inFusion` | In-Fusion cloning |

**Nested elements:**

- `<Oligo>` — Primer used in PCR: `name`, `sequence`, `phosphorylated` (optional, `"0"` or `"1"`)
- `<InputSummary>` — Range/selection: `manipulation`, `val1`, `val2`, `name1`/`siteCount1` (enzymes)
- `<Parameter>` — Key-value pairs: `name`, `val`
- `<Node>` — Child nodes (single dict or list)

### Block 8 — Sequence Properties (XML)

Additional sequence properties. Top-level element: `<AdditionalSequenceProperties>`.

```xml
<AdditionalSequenceProperties>
  <UpstreamStickiness>0</UpstreamStickiness>
  <DownstreamStickiness>0</DownstreamStickiness>
  <UpstreamModification>Unmodified</UpstreamModification>
  <DownstreamModification>Unmodified</DownstreamModification>
</AdditionalSequenceProperties>
```

### Block 10 — Features (XML)

Annotation features. The XML `<Features>` block is parsed into a structured format with segment ranges, qualifiers, and strand mapping.

**Raw XML structure:**
```xml
<Features>
  <Feature name="GFP" type="CDS" directionality="1">
    <Segment range="101-820" color="#00ff00"/>
    <Q name="label"><V text="GFP"/></Q>
    <Q name="note"><V text="Green fluorescent protein"/></Q>
    <Q name="translation"><V text="MVSK..."/></Q>
  </Feature>
</Features>
```

**Strand mapping** (from `directionality` attribute):

| Value | Strand |
|-------|--------|
| `0` | `.` (none) |
| `1` | `+` (forward) |
| `2` | `-` (reverse) |
| `3` | `=` (both) |

**Parsed format** (after sgffp processing):
```python
{
    "features": [
        {
            "name": "GFP",
            "type": "CDS",
            "strand": "+",
            "start": 100,   # 0-based (range "101" - 1)
            "end": 820,     # 1-based end
            "color": "#00ff00",
            "segments": [{"range": "101-820", "color": "#00ff00"}],
            "qualifiers": {"label": "GFP", "note": "Green fluorescent protein"}
        }
    ]
}
```

### Block 11 — History Nodes (Binary)

Each block 11 entry stores a sequence snapshot for one history state. Multiple block 11 entries exist (one per non-root tree node). The root node's sequence lives in block 0.

**Binary layout:**

| Offset | Size | Description |
|--------|------|-------------|
| 0 | 4 | `node_index` (big-endian uint32) — links to tree node `ID` |
| 4 | 1 | `sequence_type` (see below) |
| 5+ | varies | Sequence data (type-dependent) |
| ... | varies | Nested TLV blocks (node content) |

**sequence_type values:**

| Value | Format | Description |
|-------|--------|-------------|
| 0 | Block 0 format | Uncompressed DNA (4B length + ASCII) |
| 1 | Block 1 format | Compressed DNA (recommended, SnapGene default) |
| 21 | Block 0 format | Protein sequence |
| 29 | — | Modifier-only (no sequence, empty content) |
| 32 | Block 0 format | RNA sequence |

For type 1 (compressed): the data following `sequence_type` is identical to block 1 format (4B compressed_length + 4B uncompressed_length + 14-byte metadata header + 2-bit data).

For types 0, 21, 32 (uncompressed): 4B length + ASCII sequence bytes.

For type 29: no sequence data; remaining bytes are nested TLV blocks (typically empty `{30: [{}]}`).

**Nested content** follows the sequence data as TLV blocks, stored under key `node_info`. The content is structured as `{30: [{...blocks...}]}` where the inner dict contains blocks 5, 6, 8, 10, 16, 17 (same format as top-level).

### Block 16 — Trace Container

Container wrapping a ZTR trace (block 18) with optional properties (block 8). Block 18 **never** appears as a standalone top-level block — it is always nested inside a block 16 container.

Multiple traces = multiple block 16 entries. Each container holds exactly one trace.

| Offset | Size | Description |
|--------|------|-------------|
| 0 | 4 | Flags (big-endian uint32): `0` = forward, `1` = reverse |
| 4 | N | Nested TLV blocks (block 18, optionally block 8) |

The nested data uses the same TLV format (1B type + 4B length + data).

### Block 17 — Alignable Sequences (XML)

Sequence alignment data. Top-level element: `<AlignableSequences>`.

```xml
<AlignableSequences trimStringency="Medium">
  <Sequence name="ref_seq" sequence="ATGC..."/>
</AlignableSequences>
```

### Block 18 — ZTR Trace Data

**Only appears inside block 16 containers**, never at the top level.

ZTR format (Staden Package) for Sanger sequencing chromatograms.

**Header:**

| Offset | Size | Description |
|--------|------|-------------|
| 0 | 8 | Magic: `\xaeZTR\r\n\x1a\n` |
| 8 | 2 | Version (typically `\x01\x02`) |

**Chunks** follow the header, each with:

| Field | Size | Description |
|-------|------|-------------|
| type | 4 | ASCII chunk type (e.g., `BASE`, `BPOS`) |
| metadata_length | 4 | Length of metadata (big-endian uint32) |
| metadata | N | Metadata bytes (e.g., channel name for SAMP) |
| data_length | 4 | Length of chunk data (big-endian uint32) |
| data | N | Chunk data (may be compressed) |

**Compression:** First byte of chunk data indicates format:
- `0x00` — raw (uncompressed)
- `0x02` — zlib compressed (skip 5 header bytes, decompress remainder, prepend `0x00`)

**Chunk types:**

| Type | Description | Data format (after decompression) |
|------|-------------|-----------------------------------|
| `BASE` | Base calls | Format byte + 1 padding + ASCII bases |
| `BPOS` | Base-to-sample positions | Format byte + 3 padding + big-endian uint32 per base |
| `CNF4` | Confidence scores | Format byte + 1 byte per base (confidence of called base) |
| `SMP4` | Combined ACGT samples | Format byte + 1 padding + interleaved big-endian uint16 (A, C, G, T channels sequential) |
| `SAMP` | Single channel samples | Metadata: channel letter (`A`/`C`/`G`/`T` + 3 null bytes). Data: format byte + 1 padding + big-endian uint16 per sample |
| `TEXT` | Metadata key-value pairs | Format byte + 1 padding + null-terminated key-value pairs |
| `CLIP` | Quality clip boundaries | Format byte + left (big-endian uint32) + right (big-endian uint32) |
| `COMM` | Comments | Format byte + ASCII text |

### Block 14 — Custom Enzyme Sets (XML)

User-defined enzyme groups for the restriction enzyme map. Top-level element: `<CustomEnzymeSets>`.

```xml
<CustomEnzymeSets>
  <CustomEnzymeSet type="1" name="MCS Enzymes"
    enzymeNames="Acc65I AscI BamHI EcoRI HindIII KpnI NotI SmaI XmaI"/>
</CustomEnzymeSets>
```

Empty sets: `<CustomEnzymeSets/>` (no user-defined groups).

### Block 20 — Strand Colors (XML)

Per-strand color highlighting of specific base ranges. Top-level element: `<StrandColors>`.

```xml
<StrandColors>
  <TopStrand><ColorRange range="314..354" colors="magenta"/></TopStrand>
  <BottomStrand><ColorRange range="318..358" colors="magenta"/></BottomStrand>
</StrandColors>
```

### Block 21 — Protein Sequence

Same format as Block 0. The property flags byte and ASCII sequence follow the same layout, but the sequence contains amino acid single-letter codes.

### Block 28 — Enzyme Visibilities (XML)

Override enzyme visibility in the restriction enzyme map. Top-level element: `<EnzymeVisibilities>`.

```xml
<EnzymeVisibilities vals=""/>
```

The `vals` attribute contains a comma-separated list of enzyme visibility overrides (empty when no overrides are set).

### Block 29 — History Modifier (LZMA-compressed XML)

Metadata-only history changes (no sequence modification). LZMA-compressed XML. Rarely used — most files don't have this block.

### Block 30 — History Node Content (LZMA-compressed TLV)

LZMA-compressed data containing nested TLV blocks. Used inside block 11 `node_info` to store the full file state (features, primers, notes, properties, etc.) at that history point.

After decompression, the data is parsed as a sequence of TLV blocks, producing a dict like:
```python
{
    5: [primers_data],      # Optional: Primers
    6: [notes_data],        # Optional: Notes
    8: [properties_data],   # Required: Properties
    10: [features_data],    # Optional: Features
    16: [trace_container],  # Optional: Trace containers
    17: [alignable_data],   # Required: Alignable sequences
}
```

### Block 32 — RNA Sequence

Same format as Block 0. Contains an RNA sequence (ACGU bases) with the same property flags byte.

### Block 34 — RNA Structure Predictions (LZMA-compressed JSON)

LZMA-compressed JSON containing RNA secondary structure prediction data (e.g., from ViennaRNA). Only present in files with RNA structure analysis.

After decompression, the JSON contains:

```json
{
  "bondedPairs": [[22, 30], [23, 29]],
  "parameters": {
    "noClosingGU": false,
    "noLP": false,
    "temperature": 37
  },
  "probabilityMatrix": [[0, 5, 0.077]],
  "revision": 1,
  "stats": {
    "ensembleDiversity": 213.88,
    "freeEnergyEnsemble": -97.31,
    "freeEnergyOptimal": -75.7,
    "freqOptimalEnsemble": 0.0
  },
  "suboptimalStructures": [{"bondedPairs": []}]
}
```

| Key | Type | Description |
|-----|------|-------------|
| `bondedPairs` | `[[int, int]]` | Optimal structure base pairs (0-indexed positions) |
| `parameters` | `object` | Folding parameters (temperature in °C) |
| `probabilityMatrix` | `[[i, j, p]]` | Base-pair probability matrix entries |
| `revision` | `int` | Data format revision |
| `stats` | `object` | Thermodynamic statistics (kcal/mol) |
| `suboptimalStructures` | `[object]` | Alternative folding structures |

### Blocks Not Parsed (Auto-Generated by SnapGene)

The following blocks are **not parsed** because SnapGene regenerates them on file import.

#### Block 2 — Enzyme Cut Position Index

Rare block (3 files in the SnapGene DB). Pre-computed per-enzyme cut positions, superseded by block 3.

| Offset | Size | Description |
|--------|------|-------------|
| 0 | 1 | Version (always 1) |
| 1 | 1032 | 258 × uint32 BE cut counts (indices 0–255 = enzymes from block 3, 256–257 = zero padding) |
| 1033 | N | 256 position groups, each: uint16 LE count + count × uint16 LE sorted base positions |

- `sum(cut_counts) + 3 = sequence_length` (verified all 3 files)
- Positions are **little-endian** uint16 (unlike the rest of the format which is big-endian)

#### Block 3 — Restriction Enzyme Recognition Map

Present in all files. Contains the built-in enzyme database and precomputed cut site index.

| Offset | Size | Description |
|--------|------|-------------|
| 0 | 1 | Version (always 1) |
| 1 | 4 | uint32 BE = text length (L) |
| 5 | L | Comma-separated IUPAC recognition sites (always 473 enzymes) |
| 5+L | 1+ | Null terminator(s) |
| ... | 1892 | Section 1: 473 × uint32 BE per-enzyme cut counts |
| ... | 946 | Section 2: 473 × uint16 BE per-enzyme cut counts (duplicate of section 1) |
| ... | N | Section 3: flat uint16 BE cut positions, grouped by enzyme order (total_cuts × 2 bytes) |
| ... | 28 | Tail: display/view metadata (varies per file) |

- All files share an identical 473-enzyme recognition site list (SnapGene's built-in database)
- The uint32 and uint16 count arrays are always identical (redundant)
- Position encoding: some values have high bits set — possibly strand indicator or Type IIS offset. Not fully decoded.

#### Block 13 — Enzyme Display/Filter Settings

Present in all files. Always exactly 345 bytes. Stores the enzyme map UI state.

| Offset | Size | Description |
|--------|------|-------------|
| 0 | 1 | Flag (0 or 1) |
| 1 | 1 | Flag (usually 0, 0xFF in rare cases) |
| 2 | 1 | Flag (0 or 1) |
| 3 | 1 | Always 0 |
| 4 | 1 | Always 1 |
| 5 | 1 | Always 0 |
| 6 | 1 | Flag |
| 7 | 1 | `display_limit` (75, 150, 200, or 50) |
| 8–16 | 9 | Always 0x00 |
| 17 | 1 | `filter_type`: 3=Unique 6+, 6=custom, 7=unsaved, 10=Unique Cutters |
| 18+ | N | Null-terminated ASCII filter name (e.g., "Unique 6+ Cutters") |
| ... | ... | Enzyme-specific configuration blob (fills remaining bytes to 345) |

---

## Implementation Example

```python
from sgffp import SgffObject, SgffWriter

# Create a SnapGene file using the builder pattern
sgff = (
    SgffObject.new("ATGCATGCATGCATGC", topology="circular")
    .add_feature("GFP", "CDS", 0, 8)
    .add_feature("AmpR", "CDS", 8, 16, strand="-")
    .add_primer("fwd", "ATGCATGC", bind_position=0)
)
sgff.notes.description = "My plasmid"

SgffWriter.to_file(sgff, "output.dna")
```

---

## Verification

```python
from sgffp import SgffReader

sgff = SgffReader.from_file('output.dna')

# Check basic structure
assert sgff.cookie.type_of_sequence == 1
assert sgff.sequence.length > 0
print(f"Blocks: {sorted(sgff.types)}")

# Check features
if sgff.has_features:
    for f in sgff.features:
        print(f"  {f.name} ({f.type}, {f.start}-{f.end})")

# Check history
if sgff.has_history:
    tree = sgff.history.tree
    if tree:
        for node in tree.walk():
            print(f"  [{node.id}] {node.name} ({node.operation})")
```
