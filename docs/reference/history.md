# History

## SgffHistory

Main history model managing the tree, nodes, and modifiers.

```python
from sgffp import SgffHistory
```

**BLOCK_IDS:** `(7, 11, 29, 30)`

### Tree Access (Block 7)

| Member | Description |
|--------|-------------|
| `tree → SgffHistoryTree \| None` | History tree structure (settable) |
| `get_tree_node(node_id) → SgffHistoryTreeNode \| None` | Get tree node by ID |
| `walk_tree() → Iterator[SgffHistoryTreeNode]` | Iterate all tree nodes depth-first |

### Node Access (Block 11)

| Member | Description |
|--------|-------------|
| `nodes → Dict[int, SgffHistoryNode]` | Sequence snapshots indexed by node_index |
| `get_node(index) → SgffHistoryNode \| None` | Get snapshot by index |
| `get_sequence_at(index) → str \| None` | Get sequence at a history node |
| `add_node(node) → int` | Add snapshot and sync |
| `remove_node(index) → bool` | Remove snapshot and sync |
| `update_node(index, **kwargs) → bool` | Update snapshot attributes |

### Modifiers (Block 29)

| Member | Description |
|--------|-------------|
| `modifiers → List[Dict]` | Modifier metadata |

### Tree Modification

| Method | Description |
|--------|-------------|
| `record_operation(blocks, new_seq, operation, name="", **kwargs)` | Snapshot current state and create new root |
| `snapshot_current_state(blocks) → SgffHistoryNode` | Create block 11 entry from main blocks |
| `update_for_new_sequence(new_seq)` | Update root node's `seqLen` without recording |
| `next_id() → int` | Next available tree node ID |
| `clear()` | Remove all history (blocks 7, 11, 29, 30) |

### Dunder Methods

| Method | Description |
|--------|-------------|
| `len(history)` | Number of block 11 nodes |
| `iter(history)` | Iterate over `SgffHistoryNode` values |

---

## SgffHistoryTree

Tree structure parsed from block 7.

```python
from sgffp import SgffHistoryTree
```

| Member | Description |
|--------|-------------|
| `root → SgffHistoryTreeNode \| None` | Root node (current state) |
| `nodes → Dict[int, SgffHistoryTreeNode]` | All nodes indexed by ID |
| `get(node_id) → SgffHistoryTreeNode \| None` | Lookup by ID |
| `walk(from_node=None) → Iterator` | Pre-order traversal (parent first) |
| `walk_reverse(from_node=None) → Iterator` | Post-order traversal (children first, chronological) |
| `ancestors(node_id) → List[SgffHistoryTreeNode]` | Chain from node to root |
| `to_dict() → Dict` | Serialize to dict |
| `from_dict(data) → SgffHistoryTree` | Create from parsed dict |
| `len(tree)` | Total node count |
| `iter(tree)` | Pre-order walk |

---

## SgffHistoryTreeNode

Single node in the history tree (block 7).

```python
from sgffp import SgffHistoryTreeNode
```

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `id` | `int` | Unique node ID |
| `name` | `str` | Sequence filename |
| `type` | `str` | `"DNA"`, `"RNA"`, or `"Protein"` |
| `seq_len` | `int` | Sequence length |
| `strandedness` | `str` | `"single"` or `"double"` |
| `circular` | `bool` | Circular topology |
| `operation` | `str` | Operation that created this state |
| `upstream_modification` | `str` | e.g., `"Unmodified"` |
| `downstream_modification` | `str` | e.g., `"Unmodified"` |
| `resurrectable` | `bool` | User can restore this state |

### Nested Data

| Field | Type | Description |
|-------|------|-------------|
| `oligos` | `List[SgffHistoryOligo]` | Primers used in operation |
| `parameters` | `Dict[str, str]` | Key-value parameters |
| `input_summaries` | `List[SgffInputSummary]` | Range/enzyme selections |
| `primers` | `Dict \| None` | Primer binding data |
| `history_colors` | `Dict \| None` | Strand coloring |
| `features` | `List[Dict]` | Feature snapshots |

### Tree Structure

| Field | Type | Description |
|-------|------|-------------|
| `children` | `List[SgffHistoryTreeNode]` | Child nodes |
| `parent` | `SgffHistoryTreeNode \| None` | Parent node |
| `extras` | `Dict` | Unmodeled attributes |

### Properties

| Property | Description |
|----------|-------------|
| `input_summary` | First `SgffInputSummary` (or `None`) |

Methods: `from_dict(data, parent=None)`, `to_dict()`.

---

## SgffHistoryNode

Sequence snapshot stored in block 11.

```python
from sgffp import SgffHistoryNode
```

### Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `index` | `int` | | Node index (links to tree node ID) |
| `sequence` | `str` | `""` | Sequence at this history point |
| `sequence_type` | `int` | `0` | 0=DNA, 1=compressed, 21=protein, 32=RNA |
| `length` | `int` | `0` | Sequence length |
| `content` | `SgffHistoryNodeContent \| None` | `None` | Nested content (features, etc.) |
| `format_version` | `int` | `30` | Compressed DNA format version |
| `strandedness_flag` | `int` | `1` | 1=double-stranded |
| `property_flags` | `int` | `1` | Compressed DNA property flags |
| `header_seq_length` | `int \| None` | `None` | Header sequence length |
| `tree_node` | `SgffHistoryTreeNode \| None` | `None` | Linked tree node |

### Model Accessors

Shortcuts to content models (delegate to `content`):

| Property | Type |
|----------|------|
| `features` | `SgffFeatureList` |
| `primers` | `SgffPrimerList` |
| `notes` | `SgffNotes` |
| `properties` | `SgffProperties` |
| `traces` | `SgffTraceList` |

Methods: `from_dict(data)`, `to_dict()`.

---

## SgffHistoryNodeContent

Content snapshot with model accessors. Wraps a blocks dict.

```python
from sgffp import SgffHistoryNodeContent
```

| Member | Description |
|--------|-------------|
| `blocks → Dict[int, List]` | Raw blocks dict |
| `block_types → List[int]` | Block types present |
| `exists → bool` | Any content present |
| `features` | `SgffFeatureList` |
| `primers` | `SgffPrimerList` |
| `notes` | `SgffNotes` |
| `properties` | `SgffProperties` |
| `alignments` | `SgffAlignmentList` |
| `traces` | `SgffTraceList` |
| `has_features`, `has_primers`, `has_notes`, `has_properties`, `has_alignments`, `has_traces` | Existence checks |

---

## SgffHistoryOligo

Primer/oligo used in a cloning operation.

```python
from sgffp import SgffHistoryOligo
```

| Field | Type | Default |
|-------|------|---------|
| `name` | `str` | |
| `sequence` | `str` | |
| `phosphorylated` | `bool` | `False` |

Methods: `from_dict(data)`, `to_dict()`.

---

## SgffInputSummary

Describes range/selection for an operation.

```python
from sgffp import SgffInputSummary
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `manipulation` | `str` | | e.g., `"select"`, `"insert"`, `"digest"` |
| `val1` | `int` | | Start position |
| `val2` | `int` | | End position |
| `enzymes` | `List[tuple]` | `[]` | `[(name, site_count), ...]` |
| `extras` | `Dict` | `{}` | Unmodeled attributes |

| Property | Description |
|----------|-------------|
| `enzyme_names → List[str]` | List of enzyme names |

Methods: `from_dict(data)`, `to_dict()`.

---

## HistoryOperation

String enum of known operation types.

```python
from sgffp import HistoryOperation
```

### Values by Category

**Base:**
| Enum | Value |
|------|-------|
| `INVALID` | `"invalid"` |

**Create:**
| Enum | Value |
|------|-------|
| `MAKE_DNA` | `"makeDna"` |
| `MAKE_RNA` | `"makeRna"` |
| `MAKE_PROTEIN` | `"makeProtein"` |

**Cloning:**
| Enum | Value |
|------|-------|
| `AMPLIFY` | `"amplifyFragment"` |
| `INSERT` | `"insertFragment"` |
| `INSERT_MULTI` | `"insertFragments"` |
| `DIGEST` | `"digest"` |
| `LIGATE` | `"ligateFragments"` |
| `GATEWAY_LR` | `"gatewayLRCloning"` |
| `GATEWAY_BP` | `"gatewayBPCloning"` |
| `GIBSON` | `"gibsonAssembly"` |
| `GOLDEN_GATE` | `"goldenGateAssembly"` |
| `RESTRICTION_CLONE` | `"restrictionCloning"` |
| `TA_CLONE` | `"taCloning"` |
| `TOPO_CLONE` | `"topoCloning"` |
| `IN_FUSION` | `"inFusionCloning"` |

**Edit:**
| Enum | Value |
|------|-------|
| `REPLACE` | `"replace"` |
| `FLIP` | `"flip"` |
| `NEW_FROM_SELECTION` | `"newFileFromSelection"` |
| `MUTAGENESIS` | `"primerDirectedMutagenesis"` |

**Metadata:**
| Enum | Value |
|------|-------|
| `CHANGE_METHYLATION` | `"changeMethylation"` |
| `CHANGE_PHOSPHORYLATION` | `"changePhosphorylation"` |
| `CHANGE_STRANDEDNESS` | `"changeStrandedness"` |
| `CHANGE_TOPOLOGY` | `"changeTopology"` |

Unknown operation strings are handled gracefully via `_missing_()`.

---

## SgffOps

Operations API for recording history. See the [History & Operations guide](/guide/history-operations) for usage examples.

All methods return `SgffObject` for chaining. Access via `sgff.ops`.

### Methods

**Metadata:** `change_methylation(**kw)`, `change_phosphorylation(**kw)`, `change_strandedness(**kw)`, `change_topology(**kw)`

**Create:** `make_dna(seq, *, topology=, strandedness=, **kw)`, `make_rna(seq, *, topology=, **kw)`, `make_protein(seq, **kw)`

**Edit:** `replace(seq, **kw)`, `insert_fragment(seq, **kw)`, `insert_fragments(seq, **kw)`, `flip(seq, **kw)`, `new_from_selection(seq, **kw)`, `mutagenesis(seq, **kw)`

**Cloning:** `amplify(seq, **kw)`, `digest(seq, **kw)`, `ligate(seq, **kw)`, `gateway_lr(seq, **kw)`, `gateway_bp(seq, **kw)`, `gibson(seq, **kw)`, `golden_gate(seq, **kw)`, `restriction_clone(seq, **kw)`, `ta_clone(seq, **kw)`, `topo_clone(seq, **kw)`, `in_fusion(seq, **kw)`

**Generic:** `custom(operation, new_seq, **kw)`

**Bulk:** `build_from_spec(nodes, final_sequence) → SgffObject`

**Edit node:** `edit_node(node_id, **kwargs) → SgffObject`
