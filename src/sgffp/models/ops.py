"""
First-class operations API for recording SnapGene history.

Provides ergonomic, chainable methods for common sequence operations
that automatically record history entries.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from .history import HistoryOperation

if TYPE_CHECKING:
    from ..internal import SgffObject


class SgffOps:
    """Operations API that records history automatically.

    All public methods return the parent SgffObject for chaining::

        sgff.ops.change_methylation().ops.insert_fragment(new_seq)
    """

    def __init__(self, sgff: SgffObject):
        self._sgff = sgff

    # -------------------------------------------------------------------------
    # Internal helpers
    # -------------------------------------------------------------------------

    def _ensure_history(self) -> None:
        """Bootstrap history for files without an existing tree (no block 7).

        Creates a minimal initial tree with a single root node using
        operation="invalid" (SnapGene convention for initial state).
        """
        if self._sgff.has_history and self._sgff.history.tree:
            return

        seq = self._sgff.sequence
        seq_type = "DNA"
        cookie_type = self._sgff.cookie.type_of_sequence
        if cookie_type == 7:
            seq_type = "RNA"
        elif cookie_type == 2:
            seq_type = "Protein"

        root_dict = {
            "name": "",
            "type": seq_type,
            "seqLen": str(seq.length),
            "strandedness": seq.strandedness if seq.strandedness else "double",
            "ID": "1",
            "circular": "1" if seq.topology == "circular" else "0",
            "operation": HistoryOperation.INVALID.value,
            "InputSummary": {},
        }

        tree_dict = {"HistoryTree": {"Node": root_dict}}
        history = self._sgff.history
        history._set_block(7, tree_dict)
        # Reset cached tree so it re-parses from the new block data
        history._tree = None
        history._linked = False
        _ = history.tree  # trigger parse

    def _record(self, operation: str, new_sequence: str, **tree_kwargs) -> SgffObject:
        """Record an operation in history and update the sequence.

        Args:
            operation: Operation string (e.g. "replace", "insertFragment").
            new_sequence: The new sequence after the operation.
            **tree_kwargs: Extra kwargs forwarded to record_operation
                (InputSummary, Oligo, name, type, strandedness, circular, etc.)
        """
        self._ensure_history()
        self._sgff.history.record_operation(
            self._sgff.blocks,
            new_sequence,
            operation,
            **tree_kwargs,
        )
        self._sgff.sequence.value = new_sequence
        return self._sgff

    # -------------------------------------------------------------------------
    # Metadata (no sequence change)
    # -------------------------------------------------------------------------

    def change_methylation(self, **kwargs) -> SgffObject:
        """Record a methylation change."""
        return self._record(
            HistoryOperation.CHANGE_METHYLATION.value,
            self._sgff.sequence.value,
            **kwargs,
        )

    def change_phosphorylation(self, **kwargs) -> SgffObject:
        """Record a phosphorylation change."""
        return self._record(
            HistoryOperation.CHANGE_PHOSPHORYLATION.value,
            self._sgff.sequence.value,
            **kwargs,
        )

    def change_strandedness(self, **kwargs) -> SgffObject:
        """Record a strandedness change."""
        return self._record(
            HistoryOperation.CHANGE_STRANDEDNESS.value,
            self._sgff.sequence.value,
            **kwargs,
        )

    def change_topology(self, **kwargs) -> SgffObject:
        """Record a topology change."""
        return self._record(
            HistoryOperation.CHANGE_TOPOLOGY.value,
            self._sgff.sequence.value,
            **kwargs,
        )

    # -------------------------------------------------------------------------
    # Create
    # -------------------------------------------------------------------------

    def make_dna(
        self,
        sequence: str,
        *,
        topology: str = "linear",
        strandedness: str = "double",
        **kwargs,
    ) -> SgffObject:
        """Create a new DNA sequence."""
        return self._record(
            HistoryOperation.MAKE_DNA.value,
            sequence,
            type="DNA",
            strandedness=strandedness,
            circular=topology == "circular",
            **kwargs,
        )

    def make_rna(
        self,
        sequence: str,
        *,
        topology: str = "linear",
        **kwargs,
    ) -> SgffObject:
        """Create a new RNA sequence."""
        return self._record(
            HistoryOperation.MAKE_RNA.value,
            sequence,
            type="RNA",
            strandedness="single",
            circular=topology == "circular",
            **kwargs,
        )

    def make_protein(self, sequence: str, **kwargs) -> SgffObject:
        """Create a new protein sequence."""
        return self._record(
            HistoryOperation.MAKE_PROTEIN.value,
            sequence,
            type="Protein",
            strandedness="single",
            circular=False,
            **kwargs,
        )

    # -------------------------------------------------------------------------
    # Edit
    # -------------------------------------------------------------------------

    def replace(self, new_sequence: str, **kwargs) -> SgffObject:
        """Replace the sequence."""
        return self._record(HistoryOperation.REPLACE.value, new_sequence, **kwargs)

    def insert_fragment(self, new_sequence: str, **kwargs) -> SgffObject:
        """Insert a fragment into the sequence."""
        return self._record(HistoryOperation.INSERT.value, new_sequence, **kwargs)

    def insert_fragments(self, new_sequence: str, **kwargs) -> SgffObject:
        """Insert multiple fragments into the sequence."""
        return self._record(HistoryOperation.INSERT_MULTI.value, new_sequence, **kwargs)

    def flip(self, new_sequence: str, **kwargs) -> SgffObject:
        """Flip/reverse a sequence region."""
        return self._record(HistoryOperation.FLIP.value, new_sequence, **kwargs)

    def new_from_selection(self, new_sequence: str, **kwargs) -> SgffObject:
        """Create new file from a selection."""
        return self._record(
            HistoryOperation.NEW_FROM_SELECTION.value, new_sequence, **kwargs
        )

    def mutagenesis(self, new_sequence: str, **kwargs) -> SgffObject:
        """Primer-directed mutagenesis."""
        return self._record(HistoryOperation.MUTAGENESIS.value, new_sequence, **kwargs)

    # -------------------------------------------------------------------------
    # Cloning
    # -------------------------------------------------------------------------

    def amplify(self, new_sequence: str, **kwargs) -> SgffObject:
        """Amplify a fragment."""
        return self._record(HistoryOperation.AMPLIFY.value, new_sequence, **kwargs)

    def digest(self, new_sequence: str, **kwargs) -> SgffObject:
        """Digest with restriction enzymes."""
        return self._record(HistoryOperation.DIGEST.value, new_sequence, **kwargs)

    def ligate(self, new_sequence: str, **kwargs) -> SgffObject:
        """Ligate fragments."""
        return self._record(HistoryOperation.LIGATE.value, new_sequence, **kwargs)

    def gateway_lr(self, new_sequence: str, **kwargs) -> SgffObject:
        """Gateway LR cloning."""
        return self._record(HistoryOperation.GATEWAY_LR.value, new_sequence, **kwargs)

    def gateway_bp(self, new_sequence: str, **kwargs) -> SgffObject:
        """Gateway BP cloning."""
        return self._record(HistoryOperation.GATEWAY_BP.value, new_sequence, **kwargs)

    def gibson(self, new_sequence: str, **kwargs) -> SgffObject:
        """Gibson assembly."""
        return self._record(HistoryOperation.GIBSON.value, new_sequence, **kwargs)

    def golden_gate(self, new_sequence: str, **kwargs) -> SgffObject:
        """Golden Gate assembly."""
        return self._record(HistoryOperation.GOLDEN_GATE.value, new_sequence, **kwargs)

    def restriction_clone(self, new_sequence: str, **kwargs) -> SgffObject:
        """Restriction cloning."""
        return self._record(
            HistoryOperation.RESTRICTION_CLONE.value, new_sequence, **kwargs
        )

    def ta_clone(self, new_sequence: str, **kwargs) -> SgffObject:
        """TA cloning."""
        return self._record(HistoryOperation.TA_CLONE.value, new_sequence, **kwargs)

    def topo_clone(self, new_sequence: str, **kwargs) -> SgffObject:
        """TOPO cloning."""
        return self._record(HistoryOperation.TOPO_CLONE.value, new_sequence, **kwargs)

    def in_fusion(self, new_sequence: str, **kwargs) -> SgffObject:
        """In-Fusion cloning."""
        return self._record(HistoryOperation.IN_FUSION.value, new_sequence, **kwargs)

    # -------------------------------------------------------------------------
    # Generic
    # -------------------------------------------------------------------------

    def custom(self, operation: str, new_sequence: str, **tree_kwargs) -> SgffObject:
        """Record any operation string with full control over tree kwargs."""
        return self._record(operation, new_sequence, **tree_kwargs)
