"""
End-to-end roundtrip tests: read → write → read

Note: Some block types (11 - history nodes, 7/29/30 - compressed history)
are not fully round-trippable due to complex nested structures.
Tests focus on block types that the writer supports.
"""

import pytest

from sgffp.reader import SgffReader
from sgffp.writer import SgffWriter
from sgffp.internal import SgffObject

# Block types that are fully round-trippable
ROUNDTRIP_BLOCKS = {0, 1, 5, 6, 8, 10, 17, 21, 32}

# Block types with complex structures that may not round-trip
# 7, 29, 30 = LZMA compressed history
# 11 = history node with nested blocks
# 18 = ZTR trace format (complex binary)
SKIP_BLOCKS = {7, 11, 18, 29, 30}


def filter_roundtrippable_blocks(sgff):
    """Create a copy with only round-trippable blocks"""
    filtered = SgffObject(cookie=sgff.cookie)
    for block_type, items in sgff.blocks.items():
        if block_type not in SKIP_BLOCKS:
            filtered.blocks[block_type] = items
    return filtered


def compare_sequences(original, restored):
    """Compare sequence blocks between original and restored"""
    for block_type in [0, 1, 21, 32]:
        if block_type in original.blocks:
            assert block_type in restored.blocks, f"Block {block_type} missing after roundtrip"

            orig_seq = original.blocks[block_type][0].get("sequence", "")
            rest_seq = restored.blocks[block_type][0].get("sequence", "")
            assert orig_seq == rest_seq, f"Sequence mismatch in block {block_type}"


def compare_properties(original, restored):
    """Compare sequence properties"""
    for block_type in [0, 21, 32]:
        if block_type in original.blocks:
            orig = original.blocks[block_type][0]
            rest = restored.blocks[block_type][0]

            assert orig.get("topology") == rest.get("topology")
            assert orig.get("strandedness") == rest.get("strandedness")
            assert orig.get("dam_methylated") == rest.get("dam_methylated")
            assert orig.get("dcm_methylated") == rest.get("dcm_methylated")
            assert orig.get("ecoki_methylated") == rest.get("ecoki_methylated")


# =============================================================================
# Individual File Roundtrip Tests
# =============================================================================


class TestRoundtrip:
    def test_roundtrip_test_dna(self, test_dna):
        """test.dna survives read→write→read (supported blocks)"""
        original = SgffReader.from_file(test_dna)
        filtered = filter_roundtrippable_blocks(original)
        written = SgffWriter.to_bytes(filtered)
        restored = SgffReader.from_bytes(written)

        assert restored.cookie == filtered.cookie
        assert set(restored.blocks.keys()) == set(filtered.blocks.keys())

    def test_roundtrip_test2_dna(self, test2_dna):
        """test2.dna survives read→write→read"""
        original = SgffReader.from_file(test2_dna)
        filtered = filter_roundtrippable_blocks(original)
        written = SgffWriter.to_bytes(filtered)
        restored = SgffReader.from_bytes(written)

        assert restored.cookie == filtered.cookie

    def test_roundtrip_test3_dna(self, test3_dna):
        """test3.dna survives read→write→read"""
        original = SgffReader.from_file(test3_dna)
        filtered = filter_roundtrippable_blocks(original)
        written = SgffWriter.to_bytes(filtered)
        restored = SgffReader.from_bytes(written)

        assert restored.cookie == filtered.cookie

    def test_roundtrip_test_rna(self, test_rna):
        """test.rna survives read→write→read"""
        original = SgffReader.from_file(test_rna)
        filtered = filter_roundtrippable_blocks(original)
        written = SgffWriter.to_bytes(filtered)
        restored = SgffReader.from_bytes(written)

        assert restored.cookie == filtered.cookie

    def test_roundtrip_test_prot(self, test_prot):
        """test.prot survives read→write→read"""
        original = SgffReader.from_file(test_prot)
        filtered = filter_roundtrippable_blocks(original)
        written = SgffWriter.to_bytes(filtered)
        restored = SgffReader.from_bytes(written)

        assert restored.cookie == filtered.cookie


# =============================================================================
# Content Preservation Tests
# =============================================================================


class TestRoundtripContent:
    def test_roundtrip_sequence_preserved(self, test_dna):
        """Sequence identical after roundtrip"""
        original = SgffReader.from_file(test_dna)
        filtered = filter_roundtrippable_blocks(original)
        written = SgffWriter.to_bytes(filtered)
        restored = SgffReader.from_bytes(written)

        compare_sequences(filtered, restored)

    def test_roundtrip_properties_preserved(self, test_dna):
        """Sequence properties identical after roundtrip"""
        original = SgffReader.from_file(test_dna)
        filtered = filter_roundtrippable_blocks(original)
        written = SgffWriter.to_bytes(filtered)
        restored = SgffReader.from_bytes(written)

        compare_properties(filtered, restored)

    def test_roundtrip_block_counts(self, test_dna):
        """Block counts preserved after roundtrip"""
        original = SgffReader.from_file(test_dna)
        filtered = filter_roundtrippable_blocks(original)
        written = SgffWriter.to_bytes(filtered)
        restored = SgffReader.from_bytes(written)

        for block_type in filtered.blocks:
            assert len(filtered.blocks[block_type]) == len(
                restored.blocks[block_type]
            ), f"Block {block_type} count mismatch"


# =============================================================================
# Parametrized Tests
# =============================================================================


class TestRoundtripParametrized:
    @pytest.fixture(
        params=["test_dna", "test2_dna", "test3_dna", "test_rna", "test_prot"]
    )
    def sample_file(self, request, test_dna, test2_dna, test3_dna, test_rna, test_prot):
        files = {
            "test_dna": test_dna,
            "test2_dna": test2_dna,
            "test3_dna": test3_dna,
            "test_rna": test_rna,
            "test_prot": test_prot,
        }
        return files[request.param]

    def test_cookie_preserved(self, sample_file):
        """Cookie preserved for all files"""
        original = SgffReader.from_file(sample_file)
        filtered = filter_roundtrippable_blocks(original)
        written = SgffWriter.to_bytes(filtered)
        restored = SgffReader.from_bytes(written)

        assert filtered.cookie == restored.cookie

    def test_block_types_preserved(self, sample_file):
        """Block types preserved for all files"""
        original = SgffReader.from_file(sample_file)
        filtered = filter_roundtrippable_blocks(original)
        written = SgffWriter.to_bytes(filtered)
        restored = SgffReader.from_bytes(written)

        assert set(filtered.blocks.keys()) == set(restored.blocks.keys())

    def test_readable_output(self, sample_file):
        """Written output is valid SnapGene format"""
        original = SgffReader.from_file(sample_file)
        filtered = filter_roundtrippable_blocks(original)
        written = SgffWriter.to_bytes(filtered)

        # Valid header
        assert written[:1] == b"\t"
        assert written[5:13] == b"SnapGene"

        # Readable without error
        restored = SgffReader.from_bytes(written)
        assert restored is not None


# =============================================================================
# Multiple Roundtrip Tests
# =============================================================================


class TestMultipleRoundtrips:
    def test_double_roundtrip(self, test_dna):
        """File survives two roundtrips"""
        original = SgffReader.from_file(test_dna)
        filtered = filter_roundtrippable_blocks(original)

        # First roundtrip
        written1 = SgffWriter.to_bytes(filtered)
        restored1 = SgffReader.from_bytes(written1)

        # Second roundtrip
        written2 = SgffWriter.to_bytes(restored1)
        restored2 = SgffReader.from_bytes(written2)

        assert filtered.cookie == restored2.cookie
        assert set(filtered.blocks.keys()) == set(restored2.blocks.keys())

    def test_triple_roundtrip_stable(self, test_dna):
        """Output stabilizes after multiple roundtrips"""
        original = SgffReader.from_file(test_dna)
        filtered = filter_roundtrippable_blocks(original)

        written1 = SgffWriter.to_bytes(filtered)
        restored1 = SgffReader.from_bytes(written1)

        written2 = SgffWriter.to_bytes(restored1)
        restored2 = SgffReader.from_bytes(written2)

        written3 = SgffWriter.to_bytes(restored2)

        # After stabilization, output should be identical
        assert written2 == written3
