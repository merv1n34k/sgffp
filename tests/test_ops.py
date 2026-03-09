"""
Tests for SgffOps operations API
"""

import pytest

from sgffp import SgffObject, SgffReader, SgffWriter, SgffOps
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
