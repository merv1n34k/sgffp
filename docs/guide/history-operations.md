# History & Operations

## What is SnapGene History?

SnapGene tracks a tree of cloning operations in the `.dna` file itself. The tree grows backward in time — the **root** is the current file state, and **children** are the input sequences that were combined to create it.

```
  [1] Vector.dna (4000bp, circular)
  [2] Insert.dna (1000bp, linear)
  └── insertFragment → [3] Final.dna (5000bp, circular)
```

History data lives in:
- **Block 7** — LZMA-compressed XML tree of `<Node>` elements
- **Block 11** — binary sequence snapshots for each non-root node
- **Block 29** — modifier metadata (rare)
- **Block 30** — LZMA-compressed content snapshots (features, primers, etc.)

## Reading History

```python
from sgffp import SgffReader

sgff = SgffReader.from_file("cloned.dna")

if sgff.has_history:
    tree = sgff.history.tree

    # Root = current state
    root = tree.root
    print(f"{root.name} ({root.seq_len}bp, {root.operation})")

    # Walk all nodes (pre-order: parent first)
    for node in tree.walk():
        print(f"  [{node.id}] {node.name} → {node.operation}")

    # Walk chronologically (post-order: children first)
    for node in tree.walk_reverse():
        print(f"  [{node.id}] {node.name}")

    # Access sequence snapshot at a specific node
    seq = sgff.history.get_sequence_at(index=1)
```

## SgffOps API

The `SgffOps` API provides chainable methods that automatically record history entries. Access it via `sgff.ops`:

```python
from sgffp import SgffObject

sgff = SgffObject.new("ATGCATGCATGC", topology="linear")

# Record a digest operation
sgff = sgff.ops.digest("ATGC")

# Chain multiple operations
sgff = sgff.ops.insert_fragment("ATGCATGCATGCGGGG").ops.change_topology()
```

Each method:
1. Snapshots the current state into block 11
2. Creates a new root node in the history tree
3. Updates the sequence
4. Returns `SgffObject` for chaining

### Operation Categories

**Metadata** (no sequence change):

| Method | Operation |
|--------|-----------|
| `change_methylation()` | `changeMethylation` |
| `change_phosphorylation()` | `changePhosphorylation` |
| `change_strandedness()` | `changeStrandedness` |
| `change_topology()` | `changeTopology` |

**Create:**

| Method | Operation |
|--------|-----------|
| `make_dna(seq, topology=, strandedness=)` | `makeDna` |
| `make_rna(seq, topology=)` | `makeRna` |
| `make_protein(seq)` | `makeProtein` |

**Edit:**

| Method | Operation |
|--------|-----------|
| `replace(new_seq)` | `replace` |
| `insert_fragment(new_seq)` | `insertFragment` |
| `insert_fragments(new_seq)` | `insertFragments` |
| `flip(new_seq)` | `flip` |
| `new_from_selection(new_seq)` | `newFileFromSelection` |
| `mutagenesis(new_seq)` | `primerDirectedMutagenesis` |

**Cloning:**

| Method | Operation |
|--------|-----------|
| `amplify(new_seq)` | `amplifyFragment` |
| `digest(new_seq)` | `digest` |
| `ligate(new_seq)` | `ligateFragments` |
| `gateway_lr(new_seq)` | `gatewayLRCloning` |
| `gateway_bp(new_seq)` | `gatewayBPCloning` |
| `gibson(new_seq)` | `gibsonAssembly` |
| `golden_gate(new_seq)` | `goldenGateAssembly` |
| `restriction_clone(new_seq)` | `restrictionCloning` |
| `ta_clone(new_seq)` | `taCloning` |
| `topo_clone(new_seq)` | `topoCloning` |
| `in_fusion(new_seq)` | `inFusionCloning` |

**Generic:**

```python
sgff.ops.custom("myOperation", new_sequence, name="result", circular=True)
```

### Passing Tree Metadata

All ops methods accept `**kwargs` forwarded to `record_operation()`. Common kwargs:

```python
sgff.ops.digest("ATGC",
    name="Digested.dna",
    InputSummary={
        "manipulation": "digest",
        "val1": "100",
        "val2": "200",
        "name1": "EcoRI",
        "siteCount1": "1",
    },
    Oligo={"name": "fwd_primer", "sequence": "ATGCATGC"},
)
```

## Building Trees from Spec

For complex multi-source cloning histories, use `build_from_spec()`:

```python
from sgffp import SgffObject, SgffReader

vector = SgffReader.from_file("vector.dna")
insert = SgffReader.from_file("insert.dna")

sgff = SgffObject.new()
sgff.ops.build_from_spec(
    nodes=[
        {
            "id": 1,
            "operation": "insertFragment",
            "name": "Final.dna",
            "sequence": "ATGC" * 100,
            "circular": True,
            "children": [2, 3],
            "InputSummary": {
                "manipulation": "insert",
                "val1": "100",
                "val2": "400",
            },
        },
        {"id": 2, "source": vector},   # import full history
        {"id": 3, "source": insert},   # import full history
    ],
    final_sequence="ATGC" * 100,
)
```

**Manual nodes** define the tree explicitly with `operation`, `sequence`, `children`, etc.

**Import nodes** use `"source": sgff_object` to pull in an existing file's full history subtree, remapping IDs automatically.

## Editing Nodes

Edit existing tree nodes in place with `edit_node()`:

```python
sgff.ops.edit_node(
    node_id=2,
    name="Renamed Vector",
    sequence="ATGC" * 50,         # updates block 11 snapshot + tree seqLen
    InputSummary={"manipulation": "select", "val1": "0", "val2": "100"},
)
```

Supported kwargs: `name`, `operation`, `circular`, `strandedness`, `sequence`, `InputSummary`, `Oligo`, `Parameter`.

## SnapGene Crash Rule

::: warning
Every non-leaf `<Node>` in the history tree **must** have an `<InputSummary>` child element. SnapGene segfaults (null pointer dereference) without it. An empty `<InputSummary/>` is valid.
:::

sgffp handles this automatically — both `record_operation()` and `build_from_spec()` inject empty `InputSummary` on non-leaf nodes when one isn't provided.
