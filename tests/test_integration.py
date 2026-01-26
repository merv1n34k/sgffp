"""
Integration tests for SgffObject property accessors
"""

import pytest
from pathlib import Path
from sgffp import SgffReader, SgffWriter
from sgffp.models import SgffFeature, SgffSegment, SgffHistoryNode


DATA_DIR = Path(__file__).parent / "data"


class TestSgffObjectProperties:
    @pytest.fixture
    def sgff(self):
        """Load test file"""
        filepath = DATA_DIR / "pAG32.dna"
        if not filepath.exists():
            pytest.skip("Test file not available")
        return SgffReader.from_file(filepath)

    def test_sequence_property(self, sgff):
        """Access sequence via property"""
        seq = sgff.sequence
        assert seq.length > 0
        assert seq.value.startswith("A") or seq.value.startswith("T") or seq.value.startswith("G") or seq.value.startswith("C")

    def test_features_property(self, sgff):
        """Access features via property"""
        if not sgff.has_features:
            pytest.skip("No features in test file")
        features = sgff.features
        assert len(features) > 0
        assert features[0].name

    def test_block_method(self, sgff):
        """Raw block access via block()"""
        if 0 in sgff.blocks:
            raw = sgff.block(0)
            assert "sequence" in raw

    def test_has_properties(self, sgff):
        """Check has_* properties"""
        assert isinstance(sgff.has_history, bool)
        assert isinstance(sgff.has_features, bool)
        assert isinstance(sgff.has_primers, bool)


class TestModifyViaProperty:
    def test_modify_sequence_syncs_blocks(self):
        """Modifying via property updates underlying blocks"""
        blocks = {0: [{"sequence": "ATCG", "topology": "linear"}]}
        from sgffp.internal import SgffObject, Cookie

        sgff = SgffObject(cookie=Cookie(1, 1, 1), blocks=blocks)

        sgff.sequence.value = "GGGG"
        assert sgff.blocks[0][0]["sequence"] == "GGGG"

    def test_add_feature_syncs_blocks(self):
        """Adding feature via property updates underlying blocks"""
        from sgffp.internal import SgffObject, Cookie

        blocks = {10: [{"features": []}]}
        sgff = SgffObject(cookie=Cookie(1, 1, 1), blocks=blocks)

        sgff.features.add(SgffFeature(name="Test", type="gene"))
        assert len(sgff.blocks[10][0]["features"]) == 1

    def test_history_node_modification(self):
        """Modify history via property"""
        from sgffp.internal import SgffObject, Cookie

        blocks = {
            11: [{"node_index": 0, "sequence": "AAA", "sequence_type": 0}]
        }
        sgff = SgffObject(cookie=Cookie(1, 1, 1), blocks=blocks)

        sgff.history.update_node(0, sequence="BBB")
        assert sgff.blocks[11][0]["sequence"] == "BBB"


class TestHistoryRoundtrip:
    def test_history_node_serialize_deserialize(self):
        """History node can be written and read back"""
        from io import BytesIO
        from sgffp.internal import SgffObject, Cookie
        from sgffp.writer import SgffWriter

        blocks = {
            11: [
                {
                    "node_index": 0,
                    "sequence": "ATCGATCG",
                    "sequence_type": 0,
                }
            ]
        }
        sgff = SgffObject(cookie=Cookie(1, 1, 1), blocks=blocks)

        # Write
        output = SgffWriter.to_bytes(sgff)
        assert len(output) > 0

        # Read back
        sgff2 = SgffReader.from_bytes(output)
        assert sgff2.history.get_node(0).sequence == "ATCGATCG"

    def test_compressed_history_node_roundtrip(self):
        """Compressed history node round-trips correctly"""
        from sgffp.internal import SgffObject, Cookie
        from sgffp.writer import SgffWriter

        blocks = {
            11: [
                {
                    "node_index": 5,
                    "sequence": "GATCGATCGATCGATC",
                    "sequence_type": 1,
                    "mystery": b"\x00" * 14,
                }
            ]
        }
        sgff = SgffObject(cookie=Cookie(1, 1, 1), blocks=blocks)

        output = SgffWriter.to_bytes(sgff)
        sgff2 = SgffReader.from_bytes(output)

        node = sgff2.history.get_node(5)
        assert node.sequence == "GATCGATCGATCGATC"
        assert node.sequence_type == 1

    def test_history_tree_lzma_roundtrip(self):
        """History tree (block 7) LZMA round-trips correctly"""
        from sgffp.internal import SgffObject, Cookie
        from sgffp.writer import SgffWriter

        blocks = {
            7: [{"HistoryTree": {"Node": {"id": "1", "name": "test"}}}]
        }
        sgff = SgffObject(cookie=Cookie(1, 1, 1), blocks=blocks)

        output = SgffWriter.to_bytes(sgff)
        sgff2 = SgffReader.from_bytes(output)

        tree = sgff2.history.tree
        assert "HistoryTree" in tree
