"""
Tests for internal data structures: Cookie, BlockList, SgffObject
"""

import pytest
from sgffp.internal import Cookie, BlockList, SgffObject


# =============================================================================
# Cookie Tests
# =============================================================================


class TestCookie:
    def test_cookie_creation(self):
        """Create Cookie with values"""
        cookie = Cookie(type_of_sequence=1, export_version=2, import_version=3)
        assert cookie is not None

    def test_cookie_fields(self):
        """Access all Cookie fields"""
        cookie = Cookie(type_of_sequence=1, export_version=16, import_version=8)
        assert cookie.type_of_sequence == 1
        assert cookie.export_version == 16
        assert cookie.import_version == 8

    def test_cookie_equality(self):
        """Cookies with same values are equal"""
        c1 = Cookie(type_of_sequence=1, export_version=2, import_version=3)
        c2 = Cookie(type_of_sequence=1, export_version=2, import_version=3)
        assert c1 == c2

    def test_cookie_inequality(self):
        """Cookies with different values are not equal"""
        c1 = Cookie(type_of_sequence=1, export_version=2, import_version=3)
        c2 = Cookie(type_of_sequence=2, export_version=2, import_version=3)
        assert c1 != c2


# =============================================================================
# BlockList Tests
# =============================================================================


class TestBlockList:
    def test_blocklist_first(self):
        """Get first item from BlockList"""
        bl = BlockList(0, ["first", "second", "third"])
        assert bl.first == "first"

    def test_blocklist_last(self):
        """Get last item from BlockList"""
        bl = BlockList(0, ["first", "second", "third"])
        assert bl.last == "third"

    def test_blocklist_get(self):
        """Get item by index"""
        bl = BlockList(0, ["a", "b", "c"])
        assert bl.get(0) == "a"
        assert bl.get(1) == "b"
        assert bl.get(2) == "c"

    def test_blocklist_get_negative(self):
        """Get item with negative index"""
        bl = BlockList(0, ["a", "b", "c"])
        assert bl.get(-1) == "c"

    def test_blocklist_empty(self):
        """Empty BlockList returns None for first/last/get"""
        bl = BlockList(0, [])
        assert bl.first is None
        assert bl.last is None
        assert bl.get() is None

    def test_blocklist_len(self):
        """Length operation on BlockList"""
        bl = BlockList(0, ["a", "b", "c"])
        assert len(bl) == 3

    def test_blocklist_len_empty(self):
        """Length of empty BlockList is 0"""
        bl = BlockList(0, [])
        assert len(bl) == 0

    def test_blocklist_iter(self):
        """Iterate over BlockList"""
        bl = BlockList(0, ["a", "b", "c"])
        items = list(bl)
        assert items == ["a", "b", "c"]

    def test_blocklist_getitem(self):
        """Index access with []"""
        bl = BlockList(0, ["x", "y", "z"])
        assert bl[0] == "x"
        assert bl[1] == "y"
        assert bl[2] == "z"

    def test_blocklist_type_property(self):
        """BlockList exposes its type ID"""
        bl = BlockList(10, ["item"])
        assert bl.type == 10


# =============================================================================
# SgffObject Tests
# =============================================================================


class TestSgffObject:
    @pytest.fixture
    def sample_cookie(self):
        return Cookie(type_of_sequence=1, export_version=2, import_version=3)

    @pytest.fixture
    def sample_object(self, sample_cookie):
        obj = SgffObject(cookie=sample_cookie)
        obj.blocks = {0: ["seq1"], 10: ["feat1", "feat2"]}
        return obj

    def test_sgffobject_creation(self, sample_cookie):
        """Create SgffObject with cookie"""
        obj = SgffObject(cookie=sample_cookie)
        assert obj.cookie == sample_cookie
        assert obj.blocks == {}

    def test_sgffobject_types(self, sample_object):
        """List block type IDs present"""
        types = sample_object.types
        assert 0 in types
        assert 10 in types
        assert len(types) == 2

    def test_sgffobject_type_accessor(self, sample_object):
        """Get BlockList by type ID"""
        bl = sample_object.type(10)
        assert isinstance(bl, BlockList)
        assert len(bl) == 2
        assert bl.first == "feat1"

    def test_sgffobject_type_accessor_missing(self, sample_object):
        """Accessing missing type returns empty BlockList"""
        bl = sample_object.type(99)
        assert len(bl) == 0
        assert bl.first is None

    def test_sgffobject_set(self, sample_object):
        """Append to existing block type"""
        sample_object.set(10, "feat3")
        assert len(sample_object.blocks[10]) == 3
        assert "feat3" in sample_object.blocks[10]

    def test_sgffobject_set_new_type(self, sample_object):
        """Create new block type via set"""
        sample_object.set(5, "primer1")
        assert 5 in sample_object.types
        assert sample_object.blocks[5] == ["primer1"]

    def test_sgffobject_remove(self, sample_object):
        """Remove single item from block type"""
        result = sample_object.remove(10, 0)
        assert result is True
        assert len(sample_object.blocks[10]) == 1
        assert sample_object.blocks[10][0] == "feat2"

    def test_sgffobject_remove_missing_type(self, sample_object):
        """Remove from non-existent type returns False"""
        result = sample_object.remove(99, 0)
        assert result is False

    def test_sgffobject_remove_invalid_index(self, sample_object):
        """Remove with invalid index returns False"""
        result = sample_object.remove(10, 999)
        assert result is False

    def test_sgffobject_remove_cleans_empty(self, sample_object):
        """Removing last item deletes the block type"""
        sample_object.remove(0, 0)  # Remove only item from type 0
        assert 0 not in sample_object.blocks

    def test_sgffobject_bset(self, sample_object):
        """Replace entire block content with bset"""
        sample_object.bset(10, ["new1", "new2"])
        assert sample_object.blocks[10] == ["new1", "new2"]

    def test_sgffobject_bset_single_value(self, sample_object):
        """bset wraps single value in list"""
        sample_object.bset(5, "single")
        assert sample_object.blocks[5] == ["single"]

    def test_sgffobject_bremove(self, sample_object):
        """Remove entire block type"""
        result = sample_object.bremove(10)
        assert result is True
        assert 10 not in sample_object.blocks

    def test_sgffobject_bremove_missing(self, sample_object):
        """bremove on missing type returns False"""
        result = sample_object.bremove(99)
        assert result is False


# =============================================================================
# SgffObject.new() Factory Tests
# =============================================================================


class TestSgffObjectNew:
    def test_new_default(self):
        """new() creates object with default cookie and empty blocks"""
        sgff = SgffObject.new()
        assert sgff.cookie.type_of_sequence == 1
        assert sgff.cookie.export_version == 15
        assert sgff.cookie.import_version == 19
        assert sgff.blocks == {}

    def test_new_with_sequence(self):
        """new() with sequence creates block 0"""
        sgff = SgffObject.new("ATCGATCG")
        assert 0 in sgff.blocks
        assert sgff.sequence.value == "ATCGATCG"
        assert sgff.sequence.length == 8

    def test_new_circular(self):
        """new() propagates topology and strandedness"""
        sgff = SgffObject.new("ATCG", topology="circular", strandedness="single")
        assert sgff.sequence.topology == "circular"
        assert sgff.sequence.is_circular is True
        assert sgff.sequence.strandedness == "single"

    def test_new_rna(self):
        """new() with rna uses block 32 and type_of_sequence=7"""
        sgff = SgffObject.new("AUCG", sequence_type="rna")
        assert sgff.cookie.type_of_sequence == 7
        assert 32 in sgff.blocks
        assert sgff.sequence.block_id == 32

    def test_new_protein(self):
        """new() with protein uses block 21 and type_of_sequence=2"""
        sgff = SgffObject.new("MVLSPADKTNVKAAWG", sequence_type="protein")
        assert sgff.cookie.type_of_sequence == 2
        assert 21 in sgff.blocks
        assert sgff.sequence.block_id == 21

    def test_new_invalid_sequence_type(self):
        """new() raises ValueError for invalid sequence_type"""
        with pytest.raises(ValueError, match="Invalid sequence_type"):
            SgffObject.new("ATCG", sequence_type="invalid")


# =============================================================================
# Convenience Builder Method Tests
# =============================================================================


class TestSgffObjectBuilder:
    def test_add_feature(self):
        """add_feature adds a feature to the features list"""
        sgff = SgffObject.new("ATCGATCG")
        sgff.add_feature("GFP", "CDS", 0, 8)
        assert len(sgff.features) == 1
        assert sgff.features[0].name == "GFP"
        assert sgff.features[0].type == "CDS"
        assert sgff.features[0].start == 0
        assert sgff.features[0].end == 8

    def test_add_primer(self):
        """add_primer adds a primer to the primers list"""
        sgff = SgffObject.new("ATCGATCG")
        sgff.add_primer("fwd", "ATCG", bind_position=0)
        assert len(sgff.primers) == 1
        assert sgff.primers[0].name == "fwd"
        assert sgff.primers[0].sequence == "ATCG"
        assert sgff.primers[0].bind_position == 0

    def test_chaining(self):
        """Builder methods return self for chaining"""
        sgff = (
            SgffObject.new("ATCGATCG", topology="circular")
            .add_feature("GFP", "CDS", 0, 4)
            .add_feature("AmpR", "CDS", 4, 8, strand="-")
            .add_primer("fwd", "ATCG", bind_position=0)
        )
        assert len(sgff.features) == 2
        assert len(sgff.primers) == 1
        assert sgff.features[0].name == "GFP"
        assert sgff.features[1].name == "AmpR"
        assert sgff.features[1].strand == "-"

    def test_add_feature_creates_block(self):
        """add_feature on empty object creates block 10"""
        sgff = SgffObject.new("ATCG")
        assert 10 not in sgff.blocks
        sgff.add_feature("test", "misc_feature", 0, 4)
        assert 10 in sgff.blocks

    def test_add_primer_creates_block(self):
        """add_primer on empty object creates block 5"""
        sgff = SgffObject.new("ATCG")
        assert 5 not in sgff.blocks
        sgff.add_primer("fwd", "ATCG")
        assert 5 in sgff.blocks
