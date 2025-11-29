"""
Internal data structures for SGFF representation
"""

from dataclasses import dataclass, field
from typing import Dict, Any, Optional


@dataclass
class Cookie:
    """File header metadata"""

    type_of_sequence: int
    export_version: int
    import_version: int


@dataclass
class SgffObject:
    """
    Container for SnapGene file data

    Attributes:
        cookie: File metadata
        blocks: Dict mapping block type to parsed data
                For repeated blocks, keys are "11", "11.1", "11.2" etc
    """

    cookie: Cookie
    blocks: Dict[str, Any] = field(default_factory=dict)

    def get_block(self, block_type: int, index: int = 0) -> Optional[Any]:
        """Get block by type and optional index for repeated blocks"""
        if index == 0:
            return self.blocks.get(str(block_type))
        return self.blocks.get(f"{block_type}.{index}")

    def get_all_blocks(self, block_type: int) -> Dict[str, Any]:
        """Get all blocks of a specific type"""
        result = {}
        prefix = str(block_type)
        for key, value in self.blocks.items():
            if key == prefix or key.startswith(f"{prefix}."):
                result[key] = value
        return result

    def add_block(self, block_type: int, data: Any) -> str:
        """Add a block, handles automatic numbering for duplicates"""
        key = str(block_type)
        if key not in self.blocks:
            self.blocks[key] = data
            return key

        # Find next available index
        index = 1
        while f"{key}.{index}" in self.blocks:
            index += 1

        new_key = f"{key}.{index}"
        self.blocks[new_key] = data
        return new_key
