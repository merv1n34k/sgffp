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
| 0 | 4 | `compressed_length` (big-endian uint32) — total bytes of compressed block |
| 4 | 4 | `uncompressed_length` (big-endian uint32) — number of bases |
| 8 | 14 | Mystery bytes (preserved for round-trip fidelity) |
| 22 | N | 2-bit GATC-encoded sequence data |

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

For type 1 (compressed): the data following `sequence_type` is identical to block 1 format (4B compressed_length + 4B uncompressed_length + 14 mystery bytes + 2-bit data).

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

ZTR (Zettabyte Trace Record) format for Sanger sequencing chromatograms.

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

### Block 21 — Protein Sequence

Same format as Block 0. The property flags byte and ASCII sequence follow the same layout, but the sequence contains amino acid single-letter codes.

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

---

## Implementation Example

```python
import uuid
from sgffp import SgffObject, Cookie, SgffWriter

# Create a minimal SnapGene file
cookie = Cookie(type_of_sequence=1, export_version=15, import_version=17)
sgff = SgffObject(cookie=cookie, blocks={})

# Block 0: DNA sequence (circular, double-stranded)
sgff.blocks[0] = [{
    "sequence": "ATGCATGCATGC",
    "length": 12,
    "topology": "circular",
    "strandedness": "double",
}]

# Block 6: Notes
sgff.blocks[6] = [{"Notes": {
    "UUID": str(uuid.uuid4()),
    "Type": "Synthetic",
    "CustomMapLabel": "my_plasmid.dna",
    "UseCustomMapLabel": "1",
}}]

# Block 8: Properties
sgff.blocks[8] = [{"AdditionalSequenceProperties": {
    "UpstreamStickiness": "0",
    "DownstreamStickiness": "0",
    "UpstreamModification": "Unmodified",
    "DownstreamModification": "Unmodified",
}}]

# Block 17: Alignable sequences
sgff.blocks[17] = [{"AlignableSequences": {"trimStringency": "Medium"}}]

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
