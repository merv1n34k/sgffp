"""
History models for SnapGene edit history (blocks 7, 11, 29, 30)
"""

from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional

from .base import SgffModel


@dataclass
class SgffHistoryNode:
    """
    Single history node containing a sequence snapshot.

    Attributes:
        index: Node identifier in the history tree
        sequence: DNA/RNA/Protein sequence at this point
        sequence_type: 0=DNA, 1=compressed DNA, 21=protein, 32=RNA
        info: Metadata from nested block 30 (decompressed)
    """

    index: int
    sequence: str = ""
    sequence_type: int = 0
    info: Optional[Dict] = None
    _mystery: bytes = field(default_factory=lambda: b"\x00" * 14, repr=False)

    @classmethod
    def from_dict(cls, data: Dict) -> "SgffHistoryNode":
        """Create from parsed block 11 data"""
        return cls(
            index=data.get("node_index", 0),
            sequence=data.get("sequence", ""),
            sequence_type=data.get("sequence_type", 0),
            info=data.get("node_info"),
            _mystery=data.get("mystery", b"\x00" * 14),
        )

    def to_dict(self) -> Dict:
        """Convert to dict for block storage"""
        result = {
            "node_index": self.index,
            "sequence": self.sequence,
            "sequence_type": self.sequence_type,
        }
        if self.info:
            result["node_info"] = self.info
        if self._mystery:
            result["mystery"] = self._mystery
        return result


class SgffHistory(SgffModel):
    """
    SnapGene edit history.

    Wraps blocks 7 (tree), 11 (nodes), 29 (modifiers), 30 (content).
    """

    BLOCK_IDS = (7, 11, 29, 30)

    def __init__(self, blocks: Dict[int, List[Any]]):
        super().__init__(blocks)
        self._tree: Optional[Dict] = None
        self._nodes: Optional[Dict[int, SgffHistoryNode]] = None
        self._modifiers: Optional[List[Dict]] = None

    @property
    def tree(self) -> Optional[Dict]:
        """History tree structure from block 7"""
        if self._tree is None:
            data = self._get_block(7)
            self._tree = data if data else {}
        return self._tree

    @tree.setter
    def tree(self, value: Dict) -> None:
        self._tree = value
        self._sync_tree()

    @property
    def nodes(self) -> Dict[int, SgffHistoryNode]:
        """History nodes indexed by node_index"""
        if self._nodes is None:
            self._nodes = {}
            for data in self._get_blocks(11):
                node = SgffHistoryNode.from_dict(data)
                self._nodes[node.index] = node
        return self._nodes

    @property
    def modifiers(self) -> List[Dict]:
        """Modifier metadata from block 29"""
        if self._modifiers is None:
            self._modifiers = self._get_blocks(29)
        return self._modifiers

    def get_node(self, index: int) -> Optional[SgffHistoryNode]:
        """Get node by index"""
        return self.nodes.get(index)

    def get_sequence_at(self, index: int) -> Optional[str]:
        """Get sequence at specific history node"""
        node = self.get_node(index)
        return node.sequence if node else None

    def add_node(self, node: SgffHistoryNode) -> int:
        """Add a new history node"""
        self.nodes[node.index] = node
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
            if hasattr(node, key):
                setattr(node, key, value)

        self._sync_nodes()
        return True

    def clear(self) -> None:
        """Remove all history"""
        self._tree = {}
        self._nodes = {}
        self._modifiers = []
        for bid in self.BLOCK_IDS:
            self._remove_block(bid)

    def _sync_tree(self) -> None:
        """Write tree back to block 7"""
        if self._tree:
            self._set_block(7, self._tree)
        else:
            self._remove_block(7)

    def _sync_nodes(self) -> None:
        """Write nodes back to block 11"""
        if self._nodes:
            node_dicts = [node.to_dict() for node in self._nodes.values()]
            self._set_blocks(11, node_dicts)
        else:
            self._remove_block(11)

    def _sync_modifiers(self) -> None:
        """Write modifiers back to block 29"""
        if self._modifiers:
            self._set_blocks(29, self._modifiers)
        else:
            self._remove_block(29)

    def __len__(self) -> int:
        return len(self.nodes)

    def __iter__(self):
        return iter(self.nodes.values())

    def __repr__(self) -> str:
        return f"SgffHistory(nodes={len(self.nodes)}, has_tree={bool(self.tree)})"
