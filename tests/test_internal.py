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
