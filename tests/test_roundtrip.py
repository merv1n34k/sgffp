"""
End-to-end roundtrip tests: read → write → read

Blocks 7 (history tree) and 11 (history nodes) are fully round-trippable.
Block 18 (ZTR trace) requires special handling and is tested separately.
"""

import pytest

from sgffp.reader import SgffReader
from sgffp.writer import SgffWriter
from sgffp.internal import SgffObject

# Block types that are fully round-trippable
ROUNDTRIP_BLOCKS = {0, 1, 5, 6, 7, 8, 10, 11, 17, 21, 32}

# Block types with complex structures that may not round-trip
# 18 = ZTR trace format (complex binary, tested separately)
SKIP_BLOCKS = {18}


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
            assert block_type in restored.blocks, (
                f"Block {block_type} missing after roundtrip"
            )

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


# =============================================================================
# De Novo Creation Roundtrip Tests
# =============================================================================


class TestDeNovoRoundtrip:
    def test_new_write_read_roundtrip(self):
        """Create via .new(), write to bytes, read back, verify"""
        sgff = SgffObject.new("ATCGATCGATCG")
        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)

        assert restored.cookie.type_of_sequence == 1
        assert restored.cookie.export_version == 15
        assert restored.cookie.import_version == 19
        assert restored.sequence.value == "ATCGATCGATCG"
        assert restored.sequence.strandedness == "double"

    def test_new_circular_roundtrip(self):
        """Circular topology survives roundtrip"""
        sgff = SgffObject.new("ATCGATCG", topology="circular")
        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)

        assert restored.sequence.is_circular is True
        assert restored.sequence.topology == "circular"

    def test_new_with_features_roundtrip(self):
        """Features survive roundtrip"""
        sgff = SgffObject.new("ATCGATCGATCGATCG")
        sgff.add_feature("GFP", "CDS", 0, 8)
        sgff.add_feature("AmpR", "CDS", 8, 16, strand="-")

        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)

        assert len(restored.features) == 2
        assert restored.features[0].name == "GFP"
        assert restored.features[0].type == "CDS"
        assert restored.features[0].start == 0
        assert restored.features[0].end == 8
        assert restored.features[1].name == "AmpR"
        assert restored.features[1].strand == "-"

    def test_new_with_notes_roundtrip(self):
        """Notes survive roundtrip"""
        sgff = SgffObject.new("ATCG")
        sgff.notes.description = "Test plasmid"

        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)

        assert restored.has_notes
        assert restored.notes.description == "Test plasmid"

    def test_new_with_primers_roundtrip(self):
        """Primers survive roundtrip"""
        sgff = SgffObject.new("ATCGATCG")
        sgff.add_primer("fwd", "ATCG", bind_position=0)
        sgff.add_primer("rev", "GATC", bind_position=4, bind_strand="-")

        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)

        assert len(restored.primers) == 2
        assert restored.primers[0].name == "fwd"
        assert restored.primers[0].sequence == "ATCG"
        assert restored.primers[1].name == "rev"
        assert restored.primers[1].bind_strand == "-"

    def test_chained_build_roundtrip(self):
        """Full chained build survives roundtrip"""
        sgff = (
            SgffObject.new("ATCGATCGATCGATCG", topology="circular")
            .add_feature("GFP", "CDS", 0, 8)
            .add_primer("fwd", "ATCG", bind_position=0)
        )
        sgff.notes.description = "Built with builder"

        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)

        assert restored.sequence.value == "ATCGATCGATCGATCG"
        assert restored.sequence.is_circular is True
        assert len(restored.features) == 1
        assert restored.features[0].name == "GFP"
        assert len(restored.primers) == 1
        assert restored.primers[0].name == "fwd"
        assert restored.notes.description == "Built with builder"

    def test_new_double_roundtrip_stable(self):
        """De novo file output stabilizes after first roundtrip"""
        sgff = (
            SgffObject.new("ATCGATCG", topology="circular")
            .add_feature("test", "CDS", 0, 8)
        )

        written1 = SgffWriter.to_bytes(sgff)
        restored1 = SgffReader.from_bytes(written1)
        written2 = SgffWriter.to_bytes(restored1)
        restored2 = SgffReader.from_bytes(written2)
        written3 = SgffWriter.to_bytes(restored2)

        # Content preserved across roundtrips
        assert restored1.sequence.value == "ATCGATCG"
        assert restored1.features[0].name == "test"
        # Output stabilizes after first parse
        assert written2 == written3


# =============================================================================
# _raw Removal Verification
# =============================================================================


class TestNoRawField:
    def test_no_model_has_raw_field(self):
        """No dataclass model should have _raw in its fields"""
        import dataclasses
        from sgffp import models

        model_classes = [
            models.SgffPrimer,
            models.SgffAlignment,
            models.SgffSegment,
            models.SgffFeature,
            models.SgffHistoryTreeNode,
            models.SgffHistoryOligo,
            models.SgffInputSummary,
            models.SgffHistoryNode,
        ]
        for cls in model_classes:
            if dataclasses.is_dataclass(cls):
                field_names = {f.name for f in dataclasses.fields(cls)}
                assert "_raw" not in field_names, (
                    f"{cls.__name__} still has _raw field"
                )


class TestRoundtripExtrasPreserved:
    def test_roundtrip_feature_extras_preserved(self):
        """Feature extras survive write→read roundtrip"""
        from sgffp.models import SgffFeature, SgffSegment

        sgff = SgffObject.new("ATCGATCGATCGATCG")
        feature = SgffFeature(
            name="GFP",
            type="CDS",
            strand="+",
            segments=[SgffSegment(start=0, end=8, color="#00FF00")],
            qualifiers={"note": "test"},
        )
        sgff.features.add(feature)

        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)

        assert len(restored.features) == 1
        rest = restored.features[0]
        assert rest.name == "GFP"
        assert rest.type == "CDS"
        assert rest.strand == "+"
        assert rest.qualifiers["note"] == "test"

    def test_roundtrip_primer_extras_preserved(self, test_dna):
        """Primer extras survive write→read roundtrip"""
        original = SgffReader.from_file(test_dna)
        filtered = filter_roundtrippable_blocks(original)

        if 5 not in filtered.blocks:
            pytest.skip("No primers in test file")

        written = SgffWriter.to_bytes(filtered)
        restored = SgffReader.from_bytes(written)

        orig_primers = filtered.primers
        rest_primers = restored.primers

        assert len(orig_primers) == len(rest_primers)
        for orig, rest in zip(orig_primers, rest_primers):
            assert orig.name == rest.name
            assert orig.sequence == rest.sequence

    def test_roundtrip_qualifier_types_preserved(self):
        """Int qualifiers stay as ints through roundtrip"""
        sgff = SgffObject.new("ATCGATCG")
        from sgffp.models import SgffFeature, SgffSegment

        feature = SgffFeature(
            name="test",
            type="CDS",
            strand="+",
            segments=[SgffSegment(start=0, end=8)],
            qualifiers={"codon_start": 1, "note": "test gene"},
        )
        sgff.features.add(feature)

        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)

        rest_feature = restored.features[0]
        assert rest_feature.qualifiers["codon_start"] == 1
        assert isinstance(rest_feature.qualifiers["codon_start"], int)
        assert rest_feature.qualifiers["note"] == "test gene"


# =============================================================================
# History Roundtrip Tests
# =============================================================================


class TestHistoryRoundtrip:
    def test_history_tree_roundtrip(self, test_dna):
        """Block 7 (history tree) survives full-file roundtrip"""
        original = SgffReader.from_file(test_dna)
        assert 7 in original.blocks

        written = SgffWriter.to_bytes(original)
        restored = SgffReader.from_bytes(written)

        assert 7 in restored.blocks
        assert original.blocks[7] == restored.blocks[7]

    def test_history_nodes_roundtrip(self, test_dna):
        """Block 11 (history nodes) survives full-file roundtrip"""
        original = SgffReader.from_file(test_dna)
        assert 11 in original.blocks

        written = SgffWriter.to_bytes(original)
        restored = SgffReader.from_bytes(written)

        assert 11 in restored.blocks
        assert original.blocks[11] == restored.blocks[11]

    def test_history_complex_roundtrip(self):
        """pIB2 file (9 tree nodes, 7 history nodes) roundtrips correctly"""
        from pathlib import Path

        pib2 = Path("data/samples/pIB2-SEC13-mEGFP.dna")
        if not pib2.exists():
            pytest.skip("pIB2 sample file not available")

        original = SgffReader.from_file(pib2)
        assert len(original.history.tree) == 9
        assert len(original.history.nodes) == 7

        written = SgffWriter.to_bytes(original)
        restored = SgffReader.from_bytes(written)

        assert len(restored.history.tree) == 9
        assert len(restored.history.nodes) == 7
        assert original.blocks[7] == restored.blocks[7]
        assert original.blocks[11] == restored.blocks[11]

    def test_history_roundtrip_stability(self, test_dna):
        """History output stable after 2nd roundtrip"""
        original = SgffReader.from_file(test_dna)
        written1 = SgffWriter.to_bytes(original)
        restored1 = SgffReader.from_bytes(written1)
        written2 = SgffWriter.to_bytes(restored1)

        assert written1 == written2

    def test_block_29_synthetic_roundtrip(self):
        """Block 29 (LZMA XML modifier) roundtrips via synthetic data"""
        modifier = {
            "HistoryModifier": {
                "Node": {
                    "ID": "1",
                    "name": "test.dna",
                    "type": "DNA",
                    "seqLen": "100",
                }
            }
        }
        sgff = SgffObject.new("ATCGATCG")
        sgff.blocks[29] = [modifier]

        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)

        assert 29 in restored.blocks
        assert restored.blocks[29][0] == modifier

    def test_block_30_synthetic_roundtrip(self):
        """Block 30 (LZMA nested TLV) roundtrips via synthetic data"""
        nested = {
            6: [{"Notes": {"Note": "snapshot note"}}],
            8: [{"SequenceProperties": {"topology": "circular"}}],
        }
        sgff = SgffObject.new("ATCGATCG")
        sgff.blocks[30] = [nested]

        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)

        assert 30 in restored.blocks
        content = restored.blocks[30][0]
        assert 6 in content
        assert 8 in content
