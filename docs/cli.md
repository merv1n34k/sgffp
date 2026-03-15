# CLI Reference

The `sff` command-line tool provides utilities for inspecting and manipulating SnapGene `.dna` files.

```bash
sff --version          # show version
sff <command> --help   # command-specific help
```

All read commands accept stdin by default — omit the input argument or pass `-` explicitly.

---

## `sff parse`

Export a SnapGene file to JSON.

```bash
sff parse [input] [-o OUTPUT]
```

| Argument | Description |
|----------|-------------|
| `input` | Input `.dna` file (default: stdin) |
| `-o`, `--output` | Write JSON to file instead of stdout |

### Examples

```bash
sff parse plasmid.dna                     # print JSON to stdout
sff parse plasmid.dna -o plasmid.json     # write to file
cat plasmid.dna | sff parse               # read from stdin
```

---

## `sff info`

Show file information summary.

```bash
sff info [input] [-v]
```

| Argument | Description |
|----------|-------------|
| `input` | Input `.dna` file (default: stdin) |
| `-v`, `--verbose` | Show detailed contents |

### Output

Basic mode shows: file type, sequence length/topology/strandedness, feature/primer/trace/history counts, block IDs present.

Verbose mode (`-v`) additionally shows: notes (description, created, modified), feature list with strand/type/positions, primer list, trace details, history tree.

### Examples

```bash
sff info plasmid.dna
# File: plasmid.dna
# Type: DNA (export v15, import v19)
# Sequence: 5420 bp, circular, double-stranded
# Features: 12
# Primers: 4
# History: 3 nodes
# Blocks: 0, 5, 6, 7, 8, 10, 11

sff info -v plasmid.dna
# ... (includes feature list, primer list, history tree)
```

---

## `sff tree`

Display the history tree in timeline order (chronological, leaves first).

```bash
sff tree [input] [-v]
```

| Argument | Description |
|----------|-------------|
| `input` | Input `.dna` file (default: stdin) |
| `-v`, `--verbose` | Show operation details (enzymes, oligos) |

### Output

Leaf nodes show as plain entries. Non-leaf nodes show the operation arrow. Resurrectable nodes are marked.

```bash
sff tree cloned.dna
#   [0] Vector.dna (4000bp, circular) [resurrectable]
#   [1] Insert.dna (1000bp, linear) [resurrectable]
#   └── insertFragment → [2] Final.dna (5000bp, circular)
```

Verbose mode shows `InputSummary` details (enzyme names, site counts) and oligo sequences.

---

## `sff check`

Detect unknown block types. Silent by default.

```bash
sff check [input] [-l] [-d]
```

| Argument | Description |
|----------|-------------|
| `input` | Input `.dna` file (default: stdin) |
| `-l`, `--list` | List all block types found |
| `-d`, `--dump` | Dump hex data of unknown blocks |

### Examples

```bash
sff check plasmid.dna -l
#  0:  1
#  3:  1
#  6:  1
# 10:  1
# 13:  1

sff check plasmid.dna -l -d
# (also dumps hex of any blocks not in SCHEME)
```

Unknown blocks are marked with `[NEW]`.

---

## `sff filter`

Keep specific blocks and write a new file.

```bash
sff filter <input> -k KEEP -o OUTPUT
```

| Argument | Description |
|----------|-------------|
| `input` | Input `.dna` file (required) |
| `-k`, `--keep` | Comma-separated block type IDs to keep (required) |
| `-o`, `--output` | Output `.dna` file (required) |

### Examples

```bash
# Keep only sequence and features
sff filter plasmid.dna -k 0,10 -o minimal.dna

# Keep sequence, features, and notes
sff filter plasmid.dna -k 0,6,10 -o stripped.dna
```
