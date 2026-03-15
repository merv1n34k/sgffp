---
outline: [2, 3]
---

# SnapGene .dna File Format Specification

This document describes the binary format of SnapGene `.dna` files based on reverse engineering and analysis of files produced by SnapGene versions 5.x–7.x.

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

## TLV Block Format

After the header, the file contains a sequence of TLV (Type-Length-Value) blocks:

| Field | Size | Description |
|-------|------|-------------|
| type | 1 byte | Block type ID (unsigned) |
| length | 4 bytes | Data length (big-endian uint32) |
| data | `length` bytes | Block payload |

Blocks appear in any order. Some types can appear multiple times (e.g., block 11 for each history node, block 16 for each trace). Unknown block types should be skipped.

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
| `upstreamStickiness` | string | Upstream sticky end length |
| `downstreamStickiness` | string | Downstream sticky end length |
| `upstreamModification` | string | e.g., `"Unmodified"`, `"FivePrimePhosphorylated"` |
| `downstreamModification` | string | e.g., `"Unmodified"` |
| `resurrectable` | string | `"1"` if the user can restore this state |

**Attribute order convention:** SnapGene expects attributes in this order: `name`, `type`, `seqLen`, `strandedness`, `ID`, `circular`, `[resurrectable]`, `operation`.

**Operation types:**

| Operation | Category | Description |
|-----------|----------|-------------|
| `invalid` | base | Original/imported file (leaf node) |
| `makeDna` | create | Created new DNA sequence |
| `makeRna` | create | Created new RNA sequence |
| `makeProtein` | create | Created new protein sequence |
| `amplifyFragment` | cloning | PCR amplification |
| `insertFragment` | cloning | Single fragment insertion |
| `insertFragments` | cloning | Multiple fragment insertion |
| `replace` | edit | Sequence edit/substitution |
| `digest` | cloning | Restriction digest |
| `ligateFragments` | cloning | Ligation of fragments |
| `gatewayLRCloning` | cloning | Gateway LR reaction |
| `gatewayBPCloning` | cloning | Gateway BP reaction |
| `gibsonAssembly` | cloning | Gibson assembly |
| `goldenGateAssembly` | cloning | Golden Gate assembly |
| `restrictionCloning` | cloning | Restriction cloning |
| `taCloning` | cloning | TA cloning |
| `topoCloning` | cloning | TOPO cloning |
| `inFusionCloning` | cloning | In-Fusion cloning |
| `flip` | edit | Reverse complement |
| `newFileFromSelection` | edit | Extract subsequence to new file |
| `primerDirectedMutagenesis` | edit | Site-directed mutagenesis |
| `changeMethylation` | metadata | Change methylation status |
| `changePhosphorylation` | metadata | Change phosphorylation status |
| `changeStrandedness` | metadata | Change single/double stranded |
| `changeTopology` | metadata | Change linear/circular topology |

**Nested elements:**

- `<InputSummary>` — **Required on every non-leaf node** (SnapGene segfaults without it). Contains: `manipulation`, `val1`, `val2`, `name1`/`siteCount1` (enzymes). An empty `<InputSummary/>` is valid.
- `<Oligo>` — Primer used in PCR: `name`, `sequence`, `phosphorylated` (optional)
- `<Parameter>` — Key-value pairs: `name`, `val`
- `<RegeneratedSite>` — Restriction site regenerated by cloning
- `<HistoryColors>` — Strand coloring for history map display
- `<Features>` — Snapshot of features at this history state
- `<Primers>` — Primer binding sites used in this operation
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

Annotation features with segment ranges, qualifiers, and strand mapping.

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

Each block 11 entry stores a sequence snapshot for one history state. Multiple block 11 entries exist (one per non-root tree node).

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
| 29 | — | Modifier-only (no sequence data) |
| 32 | Block 0 format | RNA sequence |

### Block 14 — Custom Enzyme Sets (XML)

User-defined enzyme groups. Top-level element: `<CustomEnzymeSets>`.

```xml
<CustomEnzymeSets>
  <CustomEnzymeSet type="1" name="MCS Enzymes"
    enzymeNames="Acc65I AscI BamHI EcoRI HindIII KpnI NotI SmaI XmaI"/>
</CustomEnzymeSets>
```

### Block 16 — Trace Container

Container wrapping a ZTR trace (block 18) with optional properties (block 8). Block 18 **never** appears as a standalone top-level block.

Multiple traces = multiple block 16 entries.

| Offset | Size | Description |
|--------|------|-------------|
| 0 | 4 | Flags (big-endian uint32): `0` = forward, `1` = reverse |
| 4 | N | Nested TLV blocks (block 18, optionally block 8) |

### Block 17 — Alignable Sequences (XML)

Sequence alignment data. Top-level element: `<AlignableSequences>`.

```xml
<AlignableSequences trimStringency="Medium">
  <Sequence name="ref_seq" sequence="ATGC..."/>
</AlignableSequences>
```

### Block 18 — ZTR Trace Data

**Only appears inside block 16 containers.**

ZTR format (Staden Package) for Sanger sequencing chromatograms.

**Header:**

| Offset | Size | Description |
|--------|------|-------------|
| 0 | 8 | Magic: `\xaeZTR\r\n\x1a\n` |
| 8 | 2 | Version (typically `\x01\x02`) |

**Chunks** follow the header:

| Field | Size | Description |
|-------|------|-------------|
| type | 4 | ASCII chunk type |
| metadata_length | 4 | Metadata length (big-endian uint32) |
| metadata | N | Metadata bytes |
| data_length | 4 | Chunk data length (big-endian uint32) |
| data | N | Chunk data (may be compressed) |

**Compression:** First byte of chunk data: `0x00` = raw, `0x02` = zlib.

**Chunk types:**

| Type | Description | Data format |
|------|-------------|-------------|
| `BASE` | Base calls | Padding + ASCII bases |
| `BPOS` | Base-to-sample positions | Padding + big-endian uint32 per base |
| `CNF4` | Confidence scores | 1 byte per base |
| `SMP4` | Combined ACGT samples | Padding + big-endian uint16 (A,C,G,T sequential) |
| `SAMP` | Single channel samples | Metadata: channel letter. Data: big-endian uint16 |
| `TEXT` | Metadata key-value pairs | Null-terminated pairs |
| `CLIP` | Quality clip boundaries | Left uint32 + right uint32 |
| `COMM` | Comments | ASCII text |

### Block 20 — Strand Colors (XML)

Per-strand color highlighting. Top-level element: `<StrandColors>`.

```xml
<StrandColors>
  <TopStrand><ColorRange range="314..354" colors="magenta"/></TopStrand>
  <BottomStrand><ColorRange range="318..358" colors="magenta"/></BottomStrand>
</StrandColors>
```

### Block 21 — Protein Sequence

Same format as Block 0. Sequence contains amino acid single-letter codes.

### Block 28 — Enzyme Visibilities (XML)

Override enzyme visibility. Top-level element: `<EnzymeVisibilities>`.

```xml
<EnzymeVisibilities vals=""/>
```

### Block 29 — History Modifier (LZMA XML)

Metadata-only history changes. LZMA-compressed XML. Rarely present.

### Block 30 — History Node Content (LZMA TLV)

LZMA-compressed nested TLV blocks inside block 11 `node_info`. Contains the full file state (features, primers, notes, etc.) at that history point.

### Block 32 — RNA Sequence

Same format as Block 0. Contains RNA sequence (ACGU bases).

### Block 34 — RNA Structure Predictions (LZMA JSON)

LZMA-compressed JSON with RNA secondary structure prediction data.

```json
{
  "bondedPairs": [[22, 30], [23, 29]],
  "parameters": {"noClosingGU": false, "noLP": false, "temperature": 37},
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

### Blocks Not Parsed

Blocks 2, 3, and 13 are auto-generated by SnapGene and not parsed:

- **Block 2** — Enzyme cut position index (rare, superseded by block 3)
- **Block 3** — Restriction enzyme recognition map (473 built-in enzymes)
- **Block 13** — Enzyme display/filter settings (always 345 bytes)
