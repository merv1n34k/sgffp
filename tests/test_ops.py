"""
Tests for SgffOps operations API
"""

import pytest

from sgffp import SgffObject, SgffReader, SgffWriter
from sgffp.models.history import HistoryOperation


def _make_sgff(seq: str = "ATCGATCG", *, with_history: bool = True) -> SgffObject:
    """Create a minimal SgffObject, optionally bootstrapping history."""
    sgff = SgffObject.new(sequence=seq)
    if with_history:
        sgff.ops._ensure_history()
    return sgff


# -------------------------------------------------------------------------
# Metadata operations (no sequence change)
# -------------------------------------------------------------------------


class TestSgffOpsMetadata:
    def test_change_methylation(self):
        sgff = _make_sgff()
        old_seq = sgff.sequence.value
        old_tree_len = len(sgff.history.tree)

        result = sgff.ops.change_methylation()

        assert result is sgff
        assert sgff.sequence.value == old_seq
        assert len(sgff.history.tree) == old_tree_len + 1
        assert sgff.history.tree.root.operation == HistoryOperation.CHANGE_METHYLATION

    def test_change_topology(self):
        sgff = _make_sgff()
        old_seq = sgff.sequence.value

        sgff.ops.change_topology()

        assert sgff.sequence.value == old_seq
        assert sgff.history.tree.root.operation == HistoryOperation.CHANGE_TOPOLOGY

    def test_change_strandedness(self):
        sgff = _make_sgff()
        sgff.ops.change_strandedness()
        assert (
            sgff.history.tree.root.operation == HistoryOperation.CHANGE_STRANDEDNESS
        )

    def test_change_phosphorylation(self):
        sgff = _make_sgff()
        sgff.ops.change_phosphorylation()
        assert (
            sgff.history.tree.root.operation
            == HistoryOperation.CHANGE_PHOSPHORYLATION
        )

    def test_metadata_chaining(self):
        sgff = _make_sgff()
        result = sgff.ops.change_methylation().ops.change_topology()
        assert result is sgff
        assert sgff.history.tree.root.operation == HistoryOperation.CHANGE_TOPOLOGY
        # Two operations added on top of the initial node
        assert len(sgff.history.tree) == 3


# -------------------------------------------------------------------------
# Edit operations
# -------------------------------------------------------------------------


class TestSgffOpsEdit:
    def test_replace(self):
        sgff = _make_sgff("ATCGATCG")
        new_seq = "GGGGCCCC"

        sgff.ops.replace(new_seq)

        assert sgff.sequence.value == new_seq
        assert sgff.history.tree.root.operation == HistoryOperation.REPLACE
        assert sgff.history.tree.root.seq_len == len(new_seq)

    def test_insert_fragment(self):
        sgff = _make_sgff("ATCG")
        new_seq = "ATCGGATCCATCG"

        sgff.ops.insert_fragment(new_seq)

        assert sgff.sequence.value == new_seq
        assert sgff.history.tree.root.operation == HistoryOperation.INSERT
        assert sgff.sequence.length == 13

    def test_flip(self):
        sgff = _make_sgff("ATCG")
        reversed_seq = "GCTA"

        sgff.ops.flip(reversed_seq)

        assert sgff.sequence.value == reversed_seq
        assert sgff.history.tree.root.operation == HistoryOperation.FLIP

    def test_mutagenesis(self):
        sgff = _make_sgff("ATCGATCG")
        mutated = "ATCAATCG"  # G→A point mutation

        sgff.ops.mutagenesis(mutated)

        assert sgff.sequence.value == mutated
        assert sgff.history.tree.root.operation == HistoryOperation.MUTAGENESIS

    def test_new_from_selection(self):
        sgff = _make_sgff("ATCGATCGATCG")
        selection = "GATCG"

        sgff.ops.new_from_selection(selection)

        assert sgff.sequence.value == selection
        assert (
            sgff.history.tree.root.operation == HistoryOperation.NEW_FROM_SELECTION
        )

    def test_insert_fragments(self):
        sgff = _make_sgff("ATCG")
        new_seq = "ATCGAAAGGGCCCTTTGATCG"

        sgff.ops.insert_fragments(new_seq)

        assert sgff.sequence.value == new_seq
        assert sgff.history.tree.root.operation == HistoryOperation.INSERT_MULTI


# -------------------------------------------------------------------------
# Cloning operations
# -------------------------------------------------------------------------


class TestSgffOpsCloning:
    def test_gibson(self):
        sgff = _make_sgff("ATCG")
        assembled = "ATCGATCGATCG"

        sgff.ops.gibson(assembled)

        assert sgff.sequence.value == assembled
        assert sgff.history.tree.root.operation == HistoryOperation.GIBSON

    def test_restriction_clone(self):
        sgff = _make_sgff("ATCG")
        cloned = "ATCGGATCCATCG"

        sgff.ops.restriction_clone(
            cloned,
            InputSummary={"manipulation": "insert", "val1": "4", "val2": "4"},
        )

        assert sgff.sequence.value == cloned
        assert (
            sgff.history.tree.root.operation == HistoryOperation.RESTRICTION_CLONE
        )

    def test_ligate(self):
        sgff = _make_sgff("ATCG")
        ligated = "ATCGATCGATCG"
        sgff.ops.ligate(ligated)
        assert sgff.history.tree.root.operation == HistoryOperation.LIGATE

    def test_amplify(self):
        sgff = _make_sgff("ATCG")
        sgff.ops.amplify("AATCGA")
        assert sgff.history.tree.root.operation == HistoryOperation.AMPLIFY

    def test_digest(self):
        sgff = _make_sgff("ATCGATCG")
        sgff.ops.digest("ATCG")
        assert sgff.history.tree.root.operation == HistoryOperation.DIGEST

    def test_gateway_lr(self):
        sgff = _make_sgff("ATCG")
        sgff.ops.gateway_lr("ATCGATCG")
        assert sgff.history.tree.root.operation == HistoryOperation.GATEWAY_LR

    def test_gateway_bp(self):
        sgff = _make_sgff("ATCG")
        sgff.ops.gateway_bp("ATCGATCG")
        assert sgff.history.tree.root.operation == HistoryOperation.GATEWAY_BP

    def test_golden_gate(self):
        sgff = _make_sgff("ATCG")
        sgff.ops.golden_gate("ATCGATCG")
        assert sgff.history.tree.root.operation == HistoryOperation.GOLDEN_GATE

    def test_ta_clone(self):
        sgff = _make_sgff("ATCG")
        sgff.ops.ta_clone("ATCGATCG")
        assert sgff.history.tree.root.operation == HistoryOperation.TA_CLONE

    def test_topo_clone(self):
        sgff = _make_sgff("ATCG")
        sgff.ops.topo_clone("ATCGATCG")
        assert sgff.history.tree.root.operation == HistoryOperation.TOPO_CLONE

    def test_in_fusion(self):
        sgff = _make_sgff("ATCG")
        sgff.ops.in_fusion("ATCGATCG")
        assert sgff.history.tree.root.operation == HistoryOperation.IN_FUSION


# -------------------------------------------------------------------------
# Create operations
# -------------------------------------------------------------------------


class TestSgffOpsCreate:
    def test_make_rna(self):
        sgff = _make_sgff("AUCG")
        sgff.ops.make_rna("AUCGAUCG")

        assert sgff.sequence.value == "AUCGAUCG"
        root = sgff.history.tree.root
        assert root.operation == HistoryOperation.MAKE_RNA
        assert root.type == "RNA"
        assert root.strandedness == "single"

    def test_make_dna_no_history(self):
        """Bootstraps history on a fresh file without existing history."""
        sgff = SgffObject.new(sequence="ATCG")
        assert not sgff.has_history

        sgff.ops.make_dna("GGCCAATT")

        assert sgff.has_history
        assert sgff.sequence.value == "GGCCAATT"
        root = sgff.history.tree.root
        assert root.operation == HistoryOperation.MAKE_DNA
        assert root.type == "DNA"

    def test_make_protein(self):
        sgff = _make_sgff("MVLK")
        sgff.ops.make_protein("MVLKAA")

        assert sgff.sequence.value == "MVLKAA"
        root = sgff.history.tree.root
        assert root.operation == HistoryOperation.MAKE_PROTEIN
        assert root.type == "Protein"
        assert root.strandedness == "single"
        assert root.circular is False


# -------------------------------------------------------------------------
# General / integration
# -------------------------------------------------------------------------


class TestSgffOpsGeneral:
    def test_custom_operation(self):
        sgff = _make_sgff("ATCG")
        sgff.ops.custom("myCustomOp", "GGCC")

        assert sgff.sequence.value == "GGCC"
        assert sgff.history.tree.root.operation == "myCustomOp"

    def test_ensure_history_bootstraps(self):
        """File without history gets one via _ensure_history."""
        sgff = SgffObject.new(sequence="ATCGATCG")
        assert not sgff.has_history

        sgff.ops._ensure_history()

        assert sgff.has_history
        assert sgff.history.tree is not None
        root = sgff.history.tree.root
        assert root.operation == HistoryOperation.INVALID
        assert root.seq_len == 8
        assert root.strandedness == "double"

    def test_ensure_history_idempotent(self):
        """Calling _ensure_history twice doesn't create duplicate trees."""
        sgff = SgffObject.new(sequence="ATCG")
        sgff.ops._ensure_history()
        tree_len = len(sgff.history.tree)

        sgff.ops._ensure_history()

        assert len(sgff.history.tree) == tree_len

    def test_ops_roundtrip(self):
        """Write → read → verify tree structure survives roundtrip."""
        sgff = _make_sgff("ATCGATCG")
        sgff.ops.insert_fragment("ATCGATCGATCG")

        data = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(data)

        assert restored.sequence.value == "ATCGATCGATCG"
        assert len(restored.history.tree) == 2
        assert restored.history.tree.root.operation == HistoryOperation.INSERT

    def test_ops_real_file(self, pib2_dna):
        """Operations on the pIB2 fixture file."""
        sgff = SgffReader.from_file(str(pib2_dna))
        old_tree_len = len(sgff.history.tree)
        old_seq = sgff.sequence.value

        sgff.ops.change_methylation()

        assert len(sgff.history.tree) == old_tree_len + 1
        assert sgff.sequence.value == old_seq
        assert sgff.history.tree.root.operation == HistoryOperation.CHANGE_METHYLATION

    def test_chaining_multiple_ops(self):
        """Chain 3+ operations."""
        sgff = _make_sgff("ATCG")
        result = (
            sgff.ops.change_methylation()
            .ops.insert_fragment("ATCGATCG")
            .ops.replace("GGCC")
        )

        assert result is sgff
        assert sgff.sequence.value == "GGCC"
        # Initial node + 3 operations = 4 total nodes
        assert len(sgff.history.tree) == 4
        assert sgff.history.tree.root.operation == HistoryOperation.REPLACE

    def test_ops_property_returns_same_instance(self):
        """ops property returns cached instance."""
        sgff = _make_sgff("ATCG")
        assert sgff.ops is sgff.ops

    def test_ops_invalidate_clears_cache(self):
        """invalidate() resets ops cache."""
        sgff = _make_sgff("ATCG")
        ops1 = sgff.ops
        sgff.invalidate()
        ops2 = sgff.ops
        assert ops1 is not ops2


# -------------------------------------------------------------------------
# build_from_spec
# -------------------------------------------------------------------------


class TestSgffOpsBuildFromSpec:
    def test_linear_tree(self):
        """3 nodes in a chain: makeDna → insert → replace."""
        sgff = SgffObject.new(sequence="")
        sgff.ops.build_from_spec(
            [
                {
                    "id": 1,
                    "operation": "replace",
                    "sequence": "GGCC",
                    "children": [2],
                },
                {
                    "id": 2,
                    "operation": "insertFragment",
                    "sequence": "ATCGATCG",
                    "children": [3],
                },
                {
                    "id": 3,
                    "operation": "makeDna",
                    "sequence": "ATCG",
                },
            ],
            final_sequence="GGCC",
        )

        assert sgff.sequence.value == "GGCC"
        tree = sgff.history.tree
        assert len(tree) == 3
        assert tree.root.id == 1
        assert tree.root.operation == "replace"
        assert len(tree.root.children) == 1
        assert tree.root.children[0].id == 2

    def test_branching_tree(self):
        """Root with 2 children (ligation of 2 fragments)."""
        sgff = SgffObject.new(sequence="")
        sgff.ops.build_from_spec(
            [
                {
                    "id": 1,
                    "operation": "ligateFragments",
                    "sequence": "ATCGATCG",
                    "name": "Ligated",
                    "children": [2, 3],
                },
                {"id": 2, "operation": "makeDna", "sequence": "ATCG"},
                {"id": 3, "operation": "makeDna", "sequence": "ATCG"},
            ],
            final_sequence="ATCGATCG",
        )

        tree = sgff.history.tree
        assert len(tree) == 3
        assert len(tree.root.children) == 2
        assert tree.root.name == "Ligated"

    def test_deep_tree(self):
        """5+ levels deep."""
        nodes = []
        for i in range(1, 7):
            node = {"id": i, "operation": "replace", "sequence": "A" * i}
            if i < 6:
                node["children"] = [i + 1]
            nodes.append(node)

        sgff = SgffObject.new(sequence="")
        sgff.ops.build_from_spec(nodes, final_sequence="AAAAAA")

        tree = sgff.history.tree
        assert len(tree) == 6
        # Walk from root to deepest leaf
        node = tree.root
        depth = 0
        while node.children:
            node = node.children[0]
            depth += 1
        assert depth == 5

    def test_build_from_spec_roundtrip(self):
        """Write → read → verify structure."""
        sgff = SgffObject.new(sequence="")
        sgff.ops.build_from_spec(
            [
                {
                    "id": 1,
                    "operation": "insertFragment",
                    "sequence": "ATCGATCG",
                    "children": [2],
                },
                {"id": 2, "operation": "makeDna", "sequence": "ATCG"},
            ],
            final_sequence="ATCGATCG",
        )

        data = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(data)

        assert restored.sequence.value == "ATCGATCG"
        tree = restored.history.tree
        assert len(tree) == 2
        assert tree.root.operation == "insertFragment"
        assert tree.root.children[0].operation == "makeDna"

    def test_build_sets_final_sequence(self):
        """Main sequence matches final_sequence."""
        sgff = SgffObject.new(sequence="")
        sgff.ops.build_from_spec(
            [{"id": 1, "operation": "makeDna", "sequence": "ATCG"}],
            final_sequence="ATCG",
        )
        assert sgff.sequence.value == "ATCG"

    def test_build_creates_block11_snapshots(self):
        """Non-root nodes have block 11 entries."""
        sgff = SgffObject.new(sequence="")
        sgff.ops.build_from_spec(
            [
                {
                    "id": 1,
                    "operation": "insertFragment",
                    "sequence": "ATCGATCG",
                    "children": [2],
                },
                {"id": 2, "operation": "makeDna", "sequence": "ATCG"},
            ],
            final_sequence="ATCGATCG",
        )

        # Node 2 (non-root) should have a block 11 snapshot
        node = sgff.history.get_node(2)
        assert node is not None
        assert node.sequence == "ATCG"
        assert node.length == 4

    def test_build_input_summary_injected(self):
        """Non-leaf nodes get InputSummary automatically."""
        sgff = SgffObject.new(sequence="")
        sgff.ops.build_from_spec(
            [
                {
                    "id": 1,
                    "operation": "insertFragment",
                    "sequence": "ATCGATCG",
                    "children": [2],
                },
                {"id": 2, "operation": "makeDna", "sequence": "ATCG"},
            ],
            final_sequence="ATCGATCG",
        )

        # Root (non-leaf) should have InputSummary
        raw = sgff.history._get_block(7)
        root_node = raw["HistoryTree"]["Node"]
        assert "InputSummary" in root_node

    def test_import_source_file(self, pib2_dna):
        """Import existing .dna file as subtree."""
        source = SgffReader.from_file(str(pib2_dna))
        source_tree_len = len(source.history.tree)

        sgff = SgffObject.new(sequence="")
        sgff.ops.build_from_spec(
            [
                {
                    "id": 1,
                    "operation": "insertFragment",
                    "sequence": source.sequence.value + "GGATCC",
                    "children": [2],
                },
                {"id": 2, "source": source},
            ],
            final_sequence=source.sequence.value + "GGATCC",
        )

        tree = sgff.history.tree
        # Root + all imported nodes from source
        assert len(tree) == 1 + source_tree_len

    def test_combine_three_files(self, pib2_dna):
        """Combine 3 SgffObjects with full history."""
        s1 = SgffObject.new(sequence="ATCG")
        s1.ops._ensure_history()
        s1.ops.insert_fragment("ATCGATCG")

        s2 = SgffObject.new(sequence="GGCC")
        s2.ops._ensure_history()

        s3 = SgffObject.new(sequence="TTAA")
        s3.ops._ensure_history()
        s3.ops.replace("TTAATTAA")

        combined_seq = "ATCGATCGGGCCTTAATTAA"
        result = SgffObject.new(sequence="")
        result.ops.build_from_spec(
            [
                {
                    "id": 100,
                    "operation": "ligateFragments",
                    "sequence": combined_seq,
                    "name": "Final",
                    "children": [200, 300, 400],
                },
                {"id": 200, "source": s1},
                {"id": 300, "source": s2},
                {"id": 400, "source": s3},
            ],
            final_sequence=combined_seq,
        )

        tree = result.history.tree
        assert tree.root.id == 100
        assert len(tree.root.children) == 3
        # s1 has 2 nodes, s2 has 1, s3 has 2 → total = 1 root + 5 imported
        assert len(tree) == 6

    def test_import_id_reassignment(self):
        """IDs don't conflict across sources."""
        # Both sources have ID 1 as root
        s1 = SgffObject.new(sequence="AAAA")
        s1.ops._ensure_history()  # root gets ID 1

        s2 = SgffObject.new(sequence="CCCC")
        s2.ops._ensure_history()  # root also gets ID 1

        sgff = SgffObject.new(sequence="")
        sgff.ops.build_from_spec(
            [
                {
                    "id": 10,
                    "operation": "ligateFragments",
                    "sequence": "AAAACCCC",
                    "children": [20, 30],
                },
                {"id": 20, "source": s1},
                {"id": 30, "source": s2},
            ],
            final_sequence="AAAACCCC",
        )

        tree = sgff.history.tree
        # All IDs should be unique
        all_ids = list(tree.nodes.keys())
        assert len(all_ids) == len(set(all_ids))

    def test_import_block11_copied(self):
        """Block 11 snapshots from sources are preserved."""
        s1 = SgffObject.new(sequence="ATCG")
        s1.ops._ensure_history()
        s1.ops.insert_fragment("ATCGATCG")
        # s1 now has a block 11 snapshot for the old root (ID 1)

        sgff = SgffObject.new(sequence="")
        sgff.ops.build_from_spec(
            [
                {
                    "id": 10,
                    "operation": "replace",
                    "sequence": "ATCGATCG",
                    "children": [20],
                },
                {"id": 20, "source": s1},
            ],
            final_sequence="ATCGATCG",
        )

        # The imported node's children should have block 11 snapshots
        assert len(sgff.history.nodes) >= 2

    def test_empty_nodes_raises(self):
        """Empty nodes list raises ValueError."""
        sgff = SgffObject.new(sequence="ATCG")
        with pytest.raises(ValueError, match="must not be empty"):
            sgff.ops.build_from_spec([], final_sequence="ATCG")

    def test_duplicate_id_raises(self):
        """Duplicate node IDs raise ValueError."""
        sgff = SgffObject.new(sequence="ATCG")
        with pytest.raises(ValueError, match="Duplicate"):
            sgff.ops.build_from_spec(
                [
                    {"id": 1, "operation": "makeDna", "sequence": "ATCG"},
                    {"id": 1, "operation": "makeDna", "sequence": "GGCC"},
                ],
                final_sequence="ATCG",
            )

    def test_multiple_roots_raises(self):
        """Multiple roots raise ValueError."""
        sgff = SgffObject.new(sequence="ATCG")
        with pytest.raises(ValueError, match="root"):
            sgff.ops.build_from_spec(
                [
                    {"id": 1, "operation": "makeDna", "sequence": "ATCG"},
                    {"id": 2, "operation": "makeDna", "sequence": "GGCC"},
                ],
                final_sequence="ATCG",
            )


# -------------------------------------------------------------------------
# edit_node
# -------------------------------------------------------------------------


class TestSgffOpsEditNode:
    def _make_two_node_tree(self):
        """Helper: returns sgff with root(id=2) → child(id=1)."""
        sgff = _make_sgff("ATCG")
        sgff.ops.insert_fragment("ATCGATCG")
        return sgff

    def test_edit_node_name(self):
        """Rename a node."""
        sgff = self._make_two_node_tree()
        root_id = sgff.history.tree.root.id

        result = sgff.ops.edit_node(root_id, name="Renamed")

        assert result is sgff
        assert sgff.history.tree.root.name == "Renamed"

    def test_edit_node_sequence(self):
        """Update block 11 snapshot + seqLen."""
        sgff = self._make_two_node_tree()
        # Edit the child node (non-root, has block 11 snapshot)
        child = sgff.history.tree.root.children[0]
        old_id = child.id

        sgff.ops.edit_node(old_id, sequence="GGCCGGCC")

        assert sgff.history.tree.get(old_id).seq_len == 8
        snapshot = sgff.history.get_node(old_id)
        assert snapshot is not None
        assert snapshot.sequence == "GGCCGGCC"
        assert snapshot.length == 8

    def test_edit_leaf_sequence(self):
        """Editing a leaf's sequence is safe."""
        sgff = self._make_two_node_tree()
        leaf = sgff.history.tree.root.children[0]
        assert len(leaf.children) == 0  # confirm it's a leaf

        sgff.ops.edit_node(leaf.id, sequence="TTTT")
        assert sgff.history.tree.get(leaf.id).seq_len == 4

    def test_edit_node_operation(self):
        """Change operation type."""
        sgff = self._make_two_node_tree()
        root_id = sgff.history.tree.root.id

        sgff.ops.edit_node(root_id, operation="gibsonAssembly")

        assert sgff.history.tree.root.operation == "gibsonAssembly"

    def test_edit_node_roundtrip(self):
        """Write → read → verify edits persist."""
        sgff = self._make_two_node_tree()
        root_id = sgff.history.tree.root.id
        sgff.ops.edit_node(root_id, name="Edited", operation="digest")

        data = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(data)

        assert restored.history.tree.root.name == "Edited"
        assert restored.history.tree.root.operation == "digest"

    def test_edit_nonexistent_node(self):
        """Raises ValueError for nonexistent node."""
        sgff = _make_sgff("ATCG")
        with pytest.raises(ValueError, match="not found"):
            sgff.ops.edit_node(9999, name="Nope")

    def test_edit_no_history(self):
        """Raises ValueError when no history exists."""
        sgff = SgffObject.new(sequence="ATCG")
        with pytest.raises(ValueError, match="No history"):
            sgff.ops.edit_node(1, name="Nope")

    def test_edit_node_circular(self):
        """Toggle circular flag."""
        sgff = self._make_two_node_tree()
        root_id = sgff.history.tree.root.id

        sgff.ops.edit_node(root_id, circular=True)
        assert sgff.history.tree.root.circular is True

        sgff.ops.edit_node(root_id, circular=False)
        assert sgff.history.tree.root.circular is False

    def test_edit_node_creates_snapshot_if_missing(self):
        """Editing sequence on a node without a snapshot creates one."""
        sgff = _make_sgff("ATCG")
        # Root node (id=1) has no block 11 snapshot
        root_id = sgff.history.tree.root.id
        assert sgff.history.get_node(root_id) is None

        sgff.ops.edit_node(root_id, sequence="GGGG")

        snapshot = sgff.history.get_node(root_id)
        assert snapshot is not None
        assert snapshot.sequence == "GGGG"
