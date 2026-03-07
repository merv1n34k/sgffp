"""
History models for SnapGene edit history
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Any, Optional, Iterator

from .base import SgffModel
from .feature import SgffFeatureList
from .primer import SgffPrimerList
from .notes import SgffNotes
from .properties import SgffProperties
from .alignment import SgffAlignmentList
from .trace import SgffTraceList


# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------


class HistoryOperation(str, Enum):
    """Known history operation types"""

    INVALID = "invalid"  # Root/original import
    MAKE_DNA = "makeDna"
    MAKE_RNA = "makeRna"
    MAKE_PROTEIN = "makeProtein"
    AMPLIFY = "amplifyFragment"
    INSERT = "insertFragment"
    REPLACE = "replace"
    DIGEST = "digest"
    LIGATE = "ligate"
    GATEWAY_LR = "gatewayLR"
    GATEWAY_BP = "gatewayBP"
    GIBSON = "gibsonAssembly"
    GOLDEN_GATE = "goldenGateAssembly"
    RESTRICTION_CLONE = "restrictionClone"
    TA_CLONE = "taClone"
    TOPO_CLONE = "topoClone"
    IN_FUSION = "inFusion"

    @classmethod
    def _missing_(cls, value: str) -> "HistoryOperation":
        """Handle unknown operations."""
        return cls.INVALID


# -----------------------------------------------------------------------------
# Tree Components (Block 7)
# -----------------------------------------------------------------------------


@dataclass
class SgffHistoryOligo:
    """Primer/oligo used in a cloning operation"""

    name: str
    sequence: str
    phosphorylated: bool = False

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffHistoryOligo":
        return cls(
            name=data.get("name", ""),
            sequence=data.get("sequence", ""),
            phosphorylated=data.get("phosphorylated") == "1",
        )

    def to_dict(self) -> Dict:
        result = {"name": self.name, "sequence": self.sequence}
        if self.phosphorylated:
            result["phosphorylated"] = "1"
        return result


_INPUT_SUMMARY_FIXED_KEYS = frozenset({"manipulation", "val1", "val2"})


@dataclass
class SgffInputSummary:
    """Describes range/selection for an operation"""

    manipulation: str  # "select", "amplify", "insert", "replace", etc.
    val1: int  # start position
    val2: int  # end position
    enzymes: List[tuple] = field(default_factory=list)  # [(name, site_count), ...]
    extras: Dict = field(default_factory=dict, repr=False)

    @property
    def enzyme_names(self) -> List[str]:
        """List of enzyme names used in this operation"""
        return [name for name, _ in self.enzymes]

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffInputSummary":
        if not data:
            return cls(manipulation="", val1=0, val2=0)
        # Parse enzyme pairs (name1/siteCount1, name2/siteCount2, ...)
        enzymes = []
        consumed_keys = set(_INPUT_SUMMARY_FIXED_KEYS)
        i = 1
        while f"name{i}" in data:
            name = data[f"name{i}"]
            site_count = int(data.get(f"siteCount{i}", 0))
            enzymes.append((name, site_count))
            consumed_keys.add(f"name{i}")
            consumed_keys.add(f"siteCount{i}")
            i += 1

        extras = {k: v for k, v in data.items() if k not in consumed_keys}

        return cls(
            manipulation=data.get("manipulation", ""),
            val1=int(data.get("val1", 0)),
            val2=int(data.get("val2", 0)),
            enzymes=enzymes,
            extras=extras,
        )

    def to_dict(self) -> Dict:
        result = dict(self.extras)
        result["manipulation"] = self.manipulation
        result["val1"] = str(self.val1)
        result["val2"] = str(self.val2)
        for i, (name, site_count) in enumerate(self.enzymes, start=1):
            result[f"name{i}"] = name
            result[f"siteCount{i}"] = str(site_count)
        return result


_TREE_NODE_KNOWN_KEYS = frozenset({
    "ID", "name", "type", "seqLen", "strandedness", "circular",
    "operation", "upstreamModification", "downstreamModification",
    "resurrectable", "Oligo", "Parameter", "InputSummary",
    "Primers", "HistoryColors", "Features", "Node",
})


@dataclass
class SgffHistoryTreeNode:
    """Single node in the history tree."""

    id: int
    name: str
    type: str  # "DNA", "RNA", "Protein"
    seq_len: int
    strandedness: str  # "single", "double"
    circular: bool
    operation: str
    upstream_modification: str
    downstream_modification: str
    resurrectable: bool = False

    # Nested data
    oligos: List[SgffHistoryOligo] = field(default_factory=list)
    parameters: Dict[str, str] = field(default_factory=dict)
    input_summaries: List[SgffInputSummary] = field(default_factory=list)
    primers: Optional[Dict] = None
    history_colors: Optional[Dict] = None
    features: List[Dict] = field(default_factory=list)

    @property
    def input_summary(self) -> Optional[SgffInputSummary]:
        """First input summary."""
        return self.input_summaries[0] if self.input_summaries else None

    # Tree structure
    children: List["SgffHistoryTreeNode"] = field(default_factory=list, repr=False)
    parent: Optional["SgffHistoryTreeNode"] = field(default=None, repr=False)

    # Unmodeled attributes
    extras: Dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_dict(
        cls, data: Dict, parent: Optional["SgffHistoryTreeNode"] = None
    ) -> "SgffHistoryTreeNode":
        """Create tree node from parsed dict."""
        # Parse oligos
        oligos = []
        oligo_data = data.get("Oligo", [])
        if isinstance(oligo_data, dict):
            oligo_data = [oligo_data]
        for o in oligo_data:
            oligos.append(SgffHistoryOligo.from_dict(o))

        # Parse parameters
        parameters = {}
        param_data = data.get("Parameter")
        if param_data:
            if isinstance(param_data, dict):
                parameters[param_data.get("name", "")] = param_data.get("val", "")
            elif isinstance(param_data, list):
                for p in param_data:
                    parameters[p.get("name", "")] = p.get("val", "")

        # Parse input summaries (can be single dict or list)
        input_summaries = []
        input_data = data.get("InputSummary")
        if input_data:
            if isinstance(input_data, list):
                input_summaries = [SgffInputSummary.from_dict(d) for d in input_data]
            else:
                input_summaries = [SgffInputSummary.from_dict(input_data)]

        # Parse features
        features = []
        features_data = data.get("Features", {})
        if features_data:
            feat_list = features_data.get("Feature", [])
            if isinstance(feat_list, dict):
                feat_list = [feat_list]
            features = feat_list

        # Extras: everything not in known keys
        extras = {k: v for k, v in data.items() if k not in _TREE_NODE_KNOWN_KEYS}

        node = cls(
            id=int(data.get("ID", 0)),
            name=data.get("name", ""),
            type=data.get("type", "DNA"),
            seq_len=int(data.get("seqLen", 0)),
            strandedness=data.get("strandedness", "double"),
            circular=data.get("circular") == "1",
            operation=data.get("operation", HistoryOperation.INVALID),
            upstream_modification=data.get("upstreamModification", "Unmodified"),
            downstream_modification=data.get("downstreamModification", "Unmodified"),
            resurrectable=data.get("resurrectable") == "1",
            oligos=oligos,
            parameters=parameters,
            input_summaries=input_summaries,
            primers=data.get("Primers"),
            history_colors=data.get("HistoryColors"),
            features=features,
            parent=parent,
            extras=extras,
        )

        # Parse child nodes recursively
        child_data = data.get("Node")
        if child_data:
            if isinstance(child_data, dict):
                child_data = [child_data]
            for child_dict in child_data:
                child = SgffHistoryTreeNode.from_dict(child_dict, parent=node)
                node.children.append(child)

        return node

    def to_dict(self) -> Dict:
        """Serialize to dict, building forward from extras + named fields."""
        result: Dict[str, Any] = dict(self.extras)

        # Named fields
        result["name"] = self.name
        result["type"] = self.type
        result["seqLen"] = str(self.seq_len)
        result["strandedness"] = self.strandedness
        result["ID"] = str(self.id)
        result["circular"] = "1" if self.circular else "0"
        result["operation"] = self.operation

        if self.upstream_modification != "Unmodified":
            result["upstreamModification"] = self.upstream_modification
        if self.downstream_modification != "Unmodified":
            result["downstreamModification"] = self.downstream_modification

        if self.resurrectable:
            result["resurrectable"] = "1"

        if self.oligos:
            if len(self.oligos) == 1:
                result["Oligo"] = self.oligos[0].to_dict()
            else:
                result["Oligo"] = [o.to_dict() for o in self.oligos]

        if self.parameters:
            params = [{"name": k, "val": v} for k, v in self.parameters.items()]
            result["Parameter"] = params[0] if len(params) == 1 else params

        if self.input_summaries:
            if len(self.input_summaries) == 1:
                result["InputSummary"] = self.input_summaries[0].to_dict()
            else:
                result["InputSummary"] = [s.to_dict() for s in self.input_summaries]

        if self.primers:
            result["Primers"] = self.primers

        if self.history_colors:
            result["HistoryColors"] = self.history_colors

        if self.features:
            result["Features"] = {
                "Feature": self.features[0]
                if len(self.features) == 1
                else self.features
            }

        # Serialize children
        if self.children:
            if len(self.children) == 1:
                result["Node"] = self.children[0].to_dict()
            else:
                result["Node"] = [c.to_dict() for c in self.children]

        return result


class SgffHistoryTree:
    """History tree structure with traversal methods."""

    def __init__(self, data: Optional[Dict] = None):
        self._data = data
        self._root: Optional[SgffHistoryTreeNode] = None
        self._nodes_by_id: Optional[Dict[int, SgffHistoryTreeNode]] = None

    def _parse(self) -> None:
        """Parse tree from raw dict data."""
        if self._root is not None:
            return

        self._nodes_by_id = {}

        if not self._data:
            return

        # Block 7 contains {"HistoryTree": {"Node": {...}}}
        tree_data = self._data.get("HistoryTree", {})
        node_data = tree_data.get("Node")

        if node_data:
            self._root = SgffHistoryTreeNode.from_dict(node_data)
            self._index_nodes(self._root)

    def _index_nodes(self, node: SgffHistoryTreeNode) -> None:
        """Build index of all nodes by ID"""
        if self._nodes_by_id is None:
            self._nodes_by_id = {}
        self._nodes_by_id[node.id] = node
        for child in node.children:
            self._index_nodes(child)

    @property
    def root(self) -> Optional[SgffHistoryTreeNode]:
        """Current state (top of tree)"""
        self._parse()
        return self._root

    @property
    def nodes(self) -> Dict[int, SgffHistoryTreeNode]:
        """All nodes indexed by ID"""
        self._parse()
        return self._nodes_by_id or {}

    def get(self, node_id: int) -> Optional[SgffHistoryTreeNode]:
        """Get node by ID"""
        return self.nodes.get(node_id)

    def walk(
        self, from_node: Optional[SgffHistoryTreeNode] = None
    ) -> Iterator[SgffHistoryTreeNode]:
        """Iterate nodes in pre-order (parent before children)."""
        start = from_node or self.root
        if not start:
            return

        yield start
        for child in start.children:
            yield from self.walk(child)

    def walk_reverse(
        self, from_node: Optional[SgffHistoryTreeNode] = None
    ) -> Iterator[SgffHistoryTreeNode]:
        """Iterate nodes in post-order (children before parent)."""
        start = from_node or self.root
        if not start:
            return

        for child in start.children:
            yield from self.walk_reverse(child)
        yield start

    def ancestors(self, node_id: int) -> List[SgffHistoryTreeNode]:
        """Get ancestor chain from node up to root."""
        result = []
        node = self.get(node_id)
        while node:
            result.append(node)
            node = node.parent
        return result

    def to_dict(self) -> Dict:
        """Serialize to dict."""
        if not self.root:
            return {}
        return {"HistoryTree": {"Node": self.root.to_dict()}}

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffHistoryTree":
        """Create tree from parsed dict."""
        return cls(data)

    def __len__(self) -> int:
        return len(self.nodes)

    def __iter__(self) -> Iterator[SgffHistoryTreeNode]:
        return self.walk()

    def __repr__(self) -> str:
        return f"SgffHistoryTree(nodes={len(self)}, root={self.root.name if self.root else None})"


# -----------------------------------------------------------------------------
# Node Content (Block 11 / Block 30)
# -----------------------------------------------------------------------------


class SgffHistoryNodeContent:
    """Content snapshot for a history node with model accessors."""

    def __init__(self, blocks: Dict[int, List[Any]] = None):
        self._blocks = blocks or {}

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffHistoryNodeContent":
        """Create from parsed node_info dict."""
        # node_info can be {30: [{blocks}]} or {8: [...], 10: [...], ...}
        content_list = data.get(30, data.get("30", []))
        if content_list:
            content = content_list[0] if content_list else {}
            return cls(content)
        # Direct blocks dict (no block 30 wrapper)
        if data:
            return cls(data)
        return cls({})

    def to_dict(self) -> Dict:
        """Serialize to dict."""
        return self._blocks

    @property
    def blocks(self) -> Dict[int, List[Any]]:
        """Raw blocks dict access"""
        return self._blocks

    @property
    def block_types(self) -> List[int]:
        """List of block types present"""
        return list(self._blocks.keys())

    # -------------------------------------------------------------------------
    # Model accessors (same pattern as SgffObject)
    # -------------------------------------------------------------------------

    @property
    def features(self) -> SgffFeatureList:
        return SgffFeatureList(self._blocks)

    @property
    def primers(self) -> SgffPrimerList:
        return SgffPrimerList(self._blocks)

    @property
    def notes(self) -> SgffNotes:
        return SgffNotes(self._blocks)

    @property
    def properties(self) -> SgffProperties:
        return SgffProperties(self._blocks)

    @property
    def alignments(self) -> SgffAlignmentList:
        return SgffAlignmentList(self._blocks)

    @property
    def traces(self) -> SgffTraceList:
        return SgffTraceList(self._blocks)

    # -------------------------------------------------------------------------
    # Existence checks (same pattern as SgffObject)
    # -------------------------------------------------------------------------

    @property
    def exists(self) -> bool:
        """Check if any content is present"""
        return bool(self._blocks)

    @property
    def has_features(self) -> bool:
        return 10 in self._blocks

    @property
    def has_primers(self) -> bool:
        return 5 in self._blocks

    @property
    def has_notes(self) -> bool:
        return 6 in self._blocks

    @property
    def has_properties(self) -> bool:
        return 8 in self._blocks

    @property
    def has_alignments(self) -> bool:
        return 17 in self._blocks

    @property
    def has_traces(self) -> bool:
        return 16 in self._blocks

    # -------------------------------------------------------------------------
    # Dunder methods
    # -------------------------------------------------------------------------

    def __bool__(self) -> bool:
        return self.exists

    def __contains__(self, block_id: int) -> bool:
        return block_id in self._blocks

    def __repr__(self) -> str:
        return f"SgffHistoryNodeContent(blocks={self.block_types})"


@dataclass
class SgffHistoryNode:
    """Single history node containing a sequence snapshot."""

    index: int
    sequence: str = ""
    sequence_type: int = 0  # 0=DNA, 1=compressed DNA, 21=protein, 32=RNA
    length: int = 0
    content: Optional[SgffHistoryNodeContent] = None

    # Compressed DNA metadata header fields
    format_version: int = 30
    strandedness_flag: int = 1
    property_flags: int = 1
    header_seq_length: Optional[int] = None

    # Set by SgffHistory after loading
    tree_node: Optional[SgffHistoryTreeNode] = field(default=None, repr=False)

    @property
    def features(self) -> SgffFeatureList:
        """Annotation features."""
        return self.content.features if self.content else SgffFeatureList({})

    @property
    def primers(self) -> SgffPrimerList:
        """Primers."""
        return self.content.primers if self.content else SgffPrimerList({})

    @property
    def notes(self) -> SgffNotes:
        """File notes."""
        return self.content.notes if self.content else SgffNotes({})

    @property
    def properties(self) -> SgffProperties:
        """Sequence properties."""
        return self.content.properties if self.content else SgffProperties({})

    @property
    def traces(self) -> SgffTraceList:
        """Sequence traces."""
        return self.content.traces if self.content else SgffTraceList({})

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffHistoryNode":
        """Create from parsed dict."""
        content = None
        node_info = data.get("node_info")
        if node_info:
            content = SgffHistoryNodeContent.from_dict(node_info)

        return cls(
            index=data.get("node_index", 0),
            sequence=data.get("sequence", ""),
            sequence_type=data.get("sequence_type", 0),
            length=data.get("length", 0),
            content=content,
            format_version=data.get("format_version", 30),
            strandedness_flag=data.get("strandedness_flag", 1),
            property_flags=data.get("property_flags", 1),
            header_seq_length=data.get("header_seq_length"),
        )

    def to_dict(self) -> Dict:
        """Serialize to dict."""
        result: Dict[str, Any] = {
            "node_index": self.index,
            "sequence": self.sequence,
            "sequence_type": self.sequence_type,
        }

        if self.length:
            result["length"] = self.length

        # Compressed DNA metadata header fields
        if self.sequence_type == 1:
            result["format_version"] = self.format_version
            result["strandedness_flag"] = self.strandedness_flag
            result["property_flags"] = self.property_flags
            result["header_seq_length"] = (
                self.header_seq_length
                if self.header_seq_length is not None
                else self.length
            )

        if self.content and self.content.exists:
            result["node_info"] = self.content.to_dict()

        return result


# -----------------------------------------------------------------------------
# Main History Model
# -----------------------------------------------------------------------------


class SgffHistory(SgffModel):
    """SnapGene edit history with tree, nodes, and modifiers."""

    BLOCK_IDS = (7, 11, 29, 30)

    def __init__(self, blocks: Dict[int, List[Any]]):
        super().__init__(blocks)
        self._tree: Optional[SgffHistoryTree] = None
        self._nodes: Optional[Dict[int, SgffHistoryNode]] = None
        self._modifiers: Optional[List[Dict]] = None
        self._linked: bool = False

    # -------------------------------------------------------------------------
    # Tree Access (Block 7)
    # -------------------------------------------------------------------------

    @property
    def tree(self) -> Optional[SgffHistoryTree]:
        """History tree structure."""
        if self._tree is None:
            data = self._get_block(7)
            self._tree = SgffHistoryTree(data) if data else None
            self._link_tree_and_nodes()
        return self._tree

    @tree.setter
    def tree(self, value: SgffHistoryTree) -> None:
        self._tree = value
        self._linked = False
        self._link_tree_and_nodes()
        self._sync_tree()

    def get_tree_node(self, node_id: int) -> Optional[SgffHistoryTreeNode]:
        """Get tree node by ID"""
        return self.tree.get(node_id) if self.tree else None

    def walk_tree(self) -> Iterator[SgffHistoryTreeNode]:
        """Iterate all tree nodes depth-first"""
        if self.tree:
            yield from self.tree.walk()

    # -------------------------------------------------------------------------
    # Node Access (Block 11)
    # -------------------------------------------------------------------------

    @property
    def nodes(self) -> Dict[int, SgffHistoryNode]:
        """Sequence snapshots indexed by node_index."""
        if self._nodes is None:
            self._nodes = {}
            for data in self._get_blocks(11):
                node = SgffHistoryNode.from_dict(data)
                self._nodes[node.index] = node
            self._link_tree_and_nodes()
        return self._nodes

    def get_node(self, index: int) -> Optional[SgffHistoryNode]:
        """Get sequence node by index"""
        return self.nodes.get(index)

    def get_sequence_at(self, index: int) -> Optional[str]:
        """Get sequence at specific history node."""
        node = self.get_node(index)
        if node and node.sequence:
            return node.sequence

        # For root node, fall back to main sequence
        if self.tree and self.tree.root and self.tree.root.id == index:
            from .sequence import SgffSequence

            main_seq = SgffSequence(self._blocks)
            if main_seq.value:
                return main_seq.value

        return None

    def add_node(self, node: SgffHistoryNode) -> int:
        """Add a new history node"""
        self.nodes[node.index] = node
        # Link to tree if available
        if self._tree:
            tree_node = self._tree.get(node.index)
            if tree_node:
                node.tree_node = tree_node
        self._sync_nodes()
        return node.index

    def remove_node(self, index: int) -> bool:
        """Remove node by index"""
        if index in self.nodes:
            del self.nodes[index]
            self._sync_nodes()
            return True
        return False

    def update_node(self, index: int, **kwargs) -> bool:
        """Update node attributes"""
        node = self.get_node(index)
        if not node:
            return False

        for key, value in kwargs.items():
            if hasattr(node, key) and not key.startswith("_"):
                setattr(node, key, value)

        self._sync_nodes()
        return True

    # -------------------------------------------------------------------------
    # Modifiers (Block 29)
    # -------------------------------------------------------------------------

    @property
    def modifiers(self) -> List[Dict]:
        """Modifier metadata."""
        if self._modifiers is None:
            self._modifiers = self._get_blocks(29)
        return self._modifiers

    # -------------------------------------------------------------------------
    # Linking
    # -------------------------------------------------------------------------

    def _link_tree_and_nodes(self) -> None:
        """Link history nodes to their tree nodes."""
        if self._linked:
            return
        if self._tree is None or self._nodes is None:
            return

        for index, node in self._nodes.items():
            tree_node = self._tree.get(index)
            if tree_node:
                node.tree_node = tree_node

        self._linked = True

    # -------------------------------------------------------------------------
    # Sync
    # -------------------------------------------------------------------------

    def _sync_tree(self) -> None:
        """Write tree to block storage."""
        if self._tree:
            self._set_block(7, self._tree.to_dict())
        else:
            self._remove_block(7)

    def _sync_nodes(self) -> None:
        """Write nodes to block storage."""
        if self._nodes:
            node_dicts = [node.to_dict() for node in self._nodes.values()]
            self._set_blocks(11, node_dicts)
        else:
            self._remove_block(11)

    def _sync_modifiers(self) -> None:
        """Write modifiers to block storage."""
        if self._modifiers:
            self._set_blocks(29, self._modifiers)
        else:
            self._remove_block(29)

    # -------------------------------------------------------------------------
    # Sequence Update
    # -------------------------------------------------------------------------

    def update_for_new_sequence(self, new_sequence: str) -> None:
        """Update history tree to reflect a modified sequence.

        Updates the root node's seqLen to match the new sequence length.
        The existing history tree structure is preserved unchanged.
        """
        if not self.tree or not self.tree.root:
            return

        self.tree.root.seq_len = len(new_sequence)
        self._sync_tree()

    # -------------------------------------------------------------------------
    # Clear
    # -------------------------------------------------------------------------

    def clear(self) -> None:
        """Remove all history."""
        self._tree = None
        self._nodes = {}
        self._modifiers = []
        self._linked = False
        for bid in self.BLOCK_IDS:
            self._remove_block(bid)

    # -------------------------------------------------------------------------
    # Dunder methods
    # -------------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self.nodes)

    def __iter__(self) -> Iterator[SgffHistoryNode]:
        return iter(self.nodes.values())

    def __repr__(self) -> str:
        tree_info = f"tree={len(self.tree)}" if self.tree else "no tree"
        return f"SgffHistory(nodes={len(self.nodes)}, {tree_info})"
