"""
First-class operations API for recording SnapGene history.

Provides ergonomic, chainable methods for common sequence operations
that automatically record history entries.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict, List

from .history import HistoryOperation, SgffHistoryNode

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

    # -------------------------------------------------------------------------
    # Bulk tree construction
    # -------------------------------------------------------------------------

    def build_from_spec(
        self, nodes: List[Dict[str, Any]], final_sequence: str
    ) -> SgffObject:
        """Build an entire history tree from a specification.

        Each node dict supports two modes:

        **Manual node** (define from scratch)::

            {"id": 1, "operation": "ligateFragments", "sequence": "ATCG...",
             "name": "Final", "children": [2, 3]}

        **Import node** (from existing SgffObject)::

            {"id": 2, "source": sgff1}

        Args:
            nodes: List of node dicts describing the tree.
            final_sequence: The final sequence to set on the SgffObject.

        Returns:
            SgffObject for chaining.
        """
        if not nodes:
            raise ValueError("nodes list must not be empty")

        # Index nodes by id
        spec_by_id: Dict[int, Dict] = {}
        for spec in nodes:
            nid = spec.get("id")
            if nid is None:
                raise ValueError("Every node must have an 'id' field")
            if nid in spec_by_id:
                raise ValueError(f"Duplicate node id: {nid}")
            spec_by_id[nid] = spec

        # Find root: the node not referenced as a child by any other node
        all_children: set = set()
        for spec in nodes:
            for cid in spec.get("children", []):
                all_children.add(cid)
        root_ids = [nid for nid in spec_by_id if nid not in all_children]
        if len(root_ids) != 1:
            raise ValueError(
                f"Expected exactly one root node, found {len(root_ids)}: {root_ids}"
            )
        root_id = root_ids[0]

        # Expand source nodes: flatten imported trees and collect snapshots
        id_counter = max(spec_by_id.keys())
        block11_snapshots: List[SgffHistoryNode] = []

        def _tree_node_to_spec(
            tree_node, new_id: int, children: List[int], sequence: str,
        ) -> Dict[str, Any]:
            """Build a spec dict from a SgffHistoryTreeNode, preserving rich metadata."""
            child_spec: Dict[str, Any] = {
                "id": new_id,
                "operation": tree_node.operation,
                "name": tree_node.name,
                "type": tree_node.type,
                "strandedness": tree_node.strandedness,
                "circular": tree_node.circular,
                "sequence": sequence,
                "children": children,
                "seq_len": tree_node.seq_len,
            }
            _merge_tree_node_metadata(child_spec, tree_node)
            return child_spec

        def _merge_tree_node_metadata(spec: Dict, tree_node) -> None:
            """Copy rich metadata from a tree node into a spec dict."""
            if tree_node.input_summaries:
                if len(tree_node.input_summaries) == 1:
                    spec.setdefault("InputSummary", tree_node.input_summaries[0].to_dict())
                else:
                    spec.setdefault("InputSummary", [s.to_dict() for s in tree_node.input_summaries])
            if tree_node.oligos:
                if len(tree_node.oligos) == 1:
                    spec.setdefault("Oligo", tree_node.oligos[0].to_dict())
                else:
                    spec.setdefault("Oligo", [o.to_dict() for o in tree_node.oligos])
            if tree_node.parameters:
                params = [{"name": k, "val": v} for k, v in tree_node.parameters.items()]
                spec.setdefault("Parameter", params[0] if len(params) == 1 else params)
            if tree_node.primers:
                spec.setdefault("Primers", tree_node.primers)
            if tree_node.history_colors:
                spec.setdefault("HistoryColors", tree_node.history_colors)
            if tree_node.features:
                feat = tree_node.features[0] if len(tree_node.features) == 1 else tree_node.features
                spec.setdefault("Features", {"Feature": feat})
            if tree_node.upstream_modification != "Unmodified":
                spec.setdefault("upstreamModification", tree_node.upstream_modification)
            if tree_node.downstream_modification != "Unmodified":
                spec.setdefault("downstreamModification", tree_node.downstream_modification)

        def _expand_source(spec: Dict) -> None:
            """Expand a source node: import its full history subtree."""
            nonlocal id_counter
            source: SgffObject = spec["source"]
            if not source.has_history or not source.history.tree:
                # Source has no history — treat as leaf with its sequence
                spec.pop("source")
                spec.setdefault("operation", HistoryOperation.INVALID.value)
                spec.setdefault("sequence", source.sequence.value)
                return

            src_tree = source.history.tree
            src_root = src_tree.root

            # Map old IDs → new IDs to avoid conflicts
            old_ids = list(src_tree.nodes.keys())
            id_map: Dict[int, int] = {}
            # The source root gets the spec's id
            id_map[src_root.id] = spec["id"]
            for old_id in old_ids:
                if old_id != src_root.id:
                    id_counter += 1
                    id_map[old_id] = id_counter

            # Copy block 11 snapshots from source with remapped IDs
            for idx, node in source.history.nodes.items():
                if idx in id_map:
                    new_node_dict = node.to_dict()
                    new_node_dict["node_index"] = id_map[idx]
                    block11_snapshots.append(SgffHistoryNode.from_dict(new_node_dict))

            # Build children list for this spec from the source tree
            def _remap_children(tree_node) -> List[int]:
                return [id_map[c.id] for c in tree_node.children]

            # Register child specs from source tree (non-root nodes)
            for old_id, tree_node in src_tree.nodes.items():
                if old_id == src_root.id:
                    continue
                new_id = id_map[old_id]
                child_spec = _tree_node_to_spec(
                    tree_node, new_id, _remap_children(tree_node),
                    source.history.get_sequence_at(old_id) or "",
                )
                spec_by_id[new_id] = child_spec

            # Update this spec from the source root
            spec.pop("source")
            spec.setdefault("operation", spec.get("operation", src_root.operation))
            spec.setdefault("name", src_root.name)
            spec.setdefault("type", src_root.type)
            spec.setdefault("strandedness", src_root.strandedness)
            spec.setdefault("circular", src_root.circular)
            spec.setdefault(
                "sequence", source.history.get_sequence_at(src_root.id)
                or source.sequence.value
            )
            spec["children"] = spec.get("children", []) + _remap_children(src_root)
            spec.setdefault("seq_len", src_root.seq_len)
            # Preserve rich metadata from source root
            _merge_tree_node_metadata(spec, src_root)

        # Expand all source nodes first
        for spec in list(spec_by_id.values()):
            if "source" in spec:
                _expand_source(spec)

        # Build tree dict recursively
        def _build_node_dict(nid: int) -> Dict[str, Any]:
            spec = spec_by_id[nid]
            seq = spec.get("sequence", "")
            seq_len = spec.get("seq_len", len(seq) if seq else 0)

            node_dict: Dict[str, Any] = {
                "name": spec.get("name", ""),
                "type": spec.get("type", "DNA"),
                "seqLen": str(seq_len),
                "strandedness": spec.get("strandedness", "double"),
                "ID": str(nid),
                "circular": "1" if spec.get("circular") else "0",
            }

            # resurrectable must appear before operation (SnapGene convention)
            if nid != root_id:
                node_dict["resurrectable"] = "1"
            node_dict["operation"] = spec.get(
                "operation", HistoryOperation.INVALID.value
            )

            # Forward optional tree attributes
            for key in (
                "upstreamModification",
                "downstreamModification",
                "InputSummary",
                "Oligo",
                "Parameter",
                "Primers",
                "HistoryColors",
                "Features",
            ):
                if key in spec:
                    node_dict[key] = spec[key]

            children_ids = spec.get("children", [])

            # Non-leaf nodes need InputSummary (SnapGene segfaults without it)
            if children_ids and "InputSummary" not in node_dict:
                node_dict["InputSummary"] = {}

            # Recurse into children
            if children_ids:
                child_dicts = [_build_node_dict(cid) for cid in children_ids]
                if len(child_dicts) == 1:
                    node_dict["Node"] = child_dicts[0]
                else:
                    node_dict["Node"] = child_dicts

            return node_dict

        root_dict = _build_node_dict(root_id)
        tree_dict = {"HistoryTree": {"Node": root_dict}}

        # Write block 7
        history = self._sgff.history
        history._set_block(7, tree_dict)
        history._tree = None
        history._linked = False
        _ = history.tree  # trigger parse

        # Create block 11 snapshots for non-root manual nodes
        for nid, spec in spec_by_id.items():
            if nid == root_id:
                continue
            # Skip if already covered by source import
            if any(s.index == nid for s in block11_snapshots):
                continue
            seq = spec.get("sequence", "")
            if not seq:
                continue
            snapshot_dict: Dict[str, Any] = {
                "node_index": nid,
                "sequence": seq,
                "sequence_type": 1,  # compressed DNA
                "length": len(seq),
                "format_version": 30,
                "strandedness_flag": 1 if spec.get("strandedness", "double") == "double" else 0,
                "property_flags": 1,
                "header_seq_length": min(len(seq), 65535),
            }
            block11_snapshots.append(SgffHistoryNode.from_dict(snapshot_dict))

        # Write block 11 snapshots
        if block11_snapshots:
            history._nodes = None  # reset cache
            for snapshot in block11_snapshots:
                history.add_node(snapshot)

        # Set final sequence and sync block 0 properties from root spec
        self._sgff.sequence.value = final_sequence
        root_spec = spec_by_id[root_id]
        if root_spec.get("circular"):
            self._sgff.sequence.topology = "circular"
        strandedness = root_spec.get("strandedness", "double")
        if strandedness:
            self._sgff.sequence.strandedness = strandedness

        return self._sgff

    # -------------------------------------------------------------------------
    # Node editing
    # -------------------------------------------------------------------------

    def edit_node(self, node_id: int, **kwargs) -> SgffObject:
        """Edit a history tree node in place.

        Supported kwargs:
            name, operation, circular, strandedness — update tree node (block 7)
            sequence — update block 11 snapshot + tree node seqLen
            InputSummary, Oligo, Parameter — update tree node nested elements

        No validation of upstream/downstream consistency — caller is
        responsible. Leaf nodes are always safe to edit.

        Args:
            node_id: The ID of the tree node to edit.
            **kwargs: Fields to update.

        Returns:
            SgffObject for chaining.

        Raises:
            ValueError: If the node doesn't exist.
        """
        history = self._sgff.history
        if not history.tree:
            raise ValueError("No history tree exists")

        tree_node = history.tree.get(node_id)
        if tree_node is None:
            raise ValueError(f"Node {node_id} not found in history tree")

        # Update tree node attributes
        if "name" in kwargs:
            tree_node.name = kwargs["name"]
        if "operation" in kwargs:
            tree_node.operation = kwargs["operation"]
        if "circular" in kwargs:
            tree_node.circular = kwargs["circular"]
        if "strandedness" in kwargs:
            tree_node.strandedness = kwargs["strandedness"]

        # Update nested elements
        if "InputSummary" in kwargs:
            from .history import SgffInputSummary

            data = kwargs["InputSummary"]
            if isinstance(data, dict):
                tree_node.input_summaries = [SgffInputSummary.from_dict(data)]
            elif isinstance(data, list):
                tree_node.input_summaries = [
                    SgffInputSummary.from_dict(d) for d in data
                ]

        if "Oligo" in kwargs:
            from .history import SgffHistoryOligo

            data = kwargs["Oligo"]
            if isinstance(data, dict):
                tree_node.oligos = [SgffHistoryOligo.from_dict(data)]
            elif isinstance(data, list):
                tree_node.oligos = [SgffHistoryOligo.from_dict(d) for d in data]

        if "Parameter" in kwargs:
            data = kwargs["Parameter"]
            if isinstance(data, dict):
                tree_node.parameters = {data.get("name", ""): data.get("val", "")}
            elif isinstance(data, list):
                tree_node.parameters = {
                    p.get("name", ""): p.get("val", "") for p in data
                }

        # Update sequence → block 11 snapshot + seqLen
        if "sequence" in kwargs:
            new_seq = kwargs["sequence"]
            tree_node.seq_len = len(new_seq)

            # Update or create block 11 snapshot
            existing = history.get_node(node_id)
            if existing:
                existing.sequence = new_seq
                existing.length = len(new_seq)
                existing.header_seq_length = min(len(new_seq), 65535)
                history._sync_nodes()
            else:
                snapshot_dict = {
                    "node_index": node_id,
                    "sequence": new_seq,
                    "sequence_type": 1,
                    "length": len(new_seq),
                    "format_version": 30,
                    "strandedness_flag": 1
                    if tree_node.strandedness == "double"
                    else 0,
                    "property_flags": 1,
                    "header_seq_length": min(len(new_seq), 65535),
                }
                history.add_node(SgffHistoryNode.from_dict(snapshot_dict))

        # Re-sync tree to block 7
        history._sync_tree()

        return self._sgff
