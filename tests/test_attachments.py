"""
Tests for block 23 (file attachments) support
"""

import struct
import zlib

from sgffp.parsers import parse_attachment
from sgffp.models.attachment import SgffAttachment, SgffAttachmentList
from sgffp.internal import SgffObject
from sgffp.writer import SgffWriter
from sgffp.reader import SgffReader
from io import BytesIO


# =============================================================================
# Helpers
# =============================================================================


def _make_manifest_bytes(files):
    """Build raw manifest block bytes from a list of file dicts."""
    file_elements = []
    for f in files:
        attrs = " ".join(f'{k}="{v}"' for k, v in f.items())
        file_elements.append(f"<File {attrs}/>")
    xml_str = f'<Files>{"".join(file_elements)}</Files>'
    compressed = zlib.compress(xml_str.encode("utf-8"))
    decompressed_size = len(xml_str.encode("utf-8"))
    return struct.pack(">I", 0) + struct.pack(">I", decompressed_size) + compressed


def _make_file_block_bytes(file_id, data):
    """Build raw file data block bytes."""
    return struct.pack(">I", file_id) + data


# =============================================================================
# TestParseAttachment
# =============================================================================


class TestParseAttachment:
    def test_parse_file_data(self):
        data = _make_file_block_bytes(1, b"hello world")
        result = parse_attachment(data)
        assert result is not None
        assert result["_type"] == "file"
        assert result["id"] == 1
        assert result["data"] == b"hello world"

    def test_parse_manifest(self):
        data = _make_manifest_bytes([
            {"id": "1", "name": "test.jpg", "size": "100", "compressible": "0"},
        ])
        result = parse_attachment(data)
        assert result is not None
        assert result["_type"] == "manifest"
        manifest = result["manifest"]
        files = manifest["Files"]["File"]
        assert files["name"] == "test.jpg"
        assert files["id"] == "1"

    def test_parse_manifest_multiple_files(self):
        data = _make_manifest_bytes([
            {"id": "1", "name": "a.jpg", "size": "100", "compressible": "0"},
            {"id": "2", "name": "b.png", "size": "200", "compressible": "1"},
        ])
        result = parse_attachment(data)
        assert result is not None
        files = result["manifest"]["Files"]["File"]
        assert isinstance(files, list)
        assert len(files) == 2

    def test_reject_short_input(self):
        assert parse_attachment(b"\x00\x01") is None

    def test_reject_short_manifest(self):
        # 4 zero bytes but no decompressed_size
        assert parse_attachment(b"\x00\x00\x00\x00") is None

    def test_file_data_empty_payload(self):
        data = _make_file_block_bytes(5, b"")
        result = parse_attachment(data)
        assert result["_type"] == "file"
        assert result["id"] == 5
        assert result["data"] == b""


# =============================================================================
# TestSgffAttachment
# =============================================================================


class TestSgffAttachment:
    def test_from_blocks_with_manifest(self):
        fb = {"_type": "file", "id": 2, "data": b"binary data"}
        entry = {"id": "2", "name": "image.jpg", "size": "11", "mtime": "2024-01-01", "compressible": "0"}
        att = SgffAttachment.from_blocks(fb, entry)
        assert att.id == 2
        assert att.name == "image.jpg"
        assert att.data == b"binary data"
        assert att.size == 11
        assert att.mtime == "2024-01-01"
        assert att.compressible == "0"

    def test_from_blocks_without_manifest(self):
        fb = {"_type": "file", "id": 3, "data": b"abc"}
        att = SgffAttachment.from_blocks(fb)
        assert att.id == 3
        assert att.name == ""
        assert att.size == 3

    def test_to_manifest_dict(self):
        att = SgffAttachment(id=1, name="doc.pdf", size=500, mtime="2024-06-01", compressible="1")
        d = att.to_manifest_dict()
        assert d["id"] == "1"
        assert d["name"] == "doc.pdf"
        assert d["size"] == "500"
        assert d["mtime"] == "2024-06-01"
        assert d["compressible"] == "1"

    def test_to_file_block(self):
        att = SgffAttachment(id=7, data=b"payload")
        fb = att.to_file_block()
        assert fb["_type"] == "file"
        assert fb["id"] == 7
        assert fb["data"] == b"payload"

    def test_extras_roundtrip(self):
        fb = {"_type": "file", "id": 1, "data": b"x"}
        entry = {"id": "1", "name": "f.txt", "size": "1", "compressible": "0", "custom": "value"}
        att = SgffAttachment.from_blocks(fb, entry)
        assert att.extras == {"custom": "value"}
        d = att.to_manifest_dict()
        assert d["custom"] == "value"

    def test_to_manifest_dict_no_mtime(self):
        att = SgffAttachment(id=1, name="f.txt", size=10, compressible="0")
        d = att.to_manifest_dict()
        assert "mtime" not in d


# =============================================================================
# TestSgffAttachmentList
# =============================================================================


class TestSgffAttachmentList:
    def test_empty(self):
        blocks = {}
        al = SgffAttachmentList(blocks)
        assert len(al) == 0
        assert list(al) == []

    def test_load_single(self):
        manifest_parsed = parse_attachment(_make_manifest_bytes([
            {"id": "1", "name": "test.bin", "size": "5", "compressible": "0"},
        ]))
        file_parsed = parse_attachment(_make_file_block_bytes(1, b"hello"))
        blocks = {23: [file_parsed, manifest_parsed]}

        al = SgffAttachmentList(blocks)
        assert len(al) == 1
        assert al[0].name == "test.bin"
        assert al[0].data == b"hello"
        assert al[0].id == 1

    def test_load_multiple(self):
        manifest = parse_attachment(_make_manifest_bytes([
            {"id": "1", "name": "a.bin", "size": "3", "compressible": "0"},
            {"id": "2", "name": "b.bin", "size": "4", "compressible": "0"},
        ]))
        f1 = parse_attachment(_make_file_block_bytes(1, b"abc"))
        f2 = parse_attachment(_make_file_block_bytes(2, b"defg"))
        blocks = {23: [f1, f2, manifest]}

        al = SgffAttachmentList(blocks)
        assert len(al) == 2
        assert al[0].name == "a.bin"
        assert al[1].name == "b.bin"

    def test_add(self):
        blocks = {}
        al = SgffAttachmentList(blocks)
        att = SgffAttachment(name="new.txt", data=b"content")
        al.add(att)
        assert len(al) == 1
        assert att.id == 1
        assert att.size == 7
        assert 23 in blocks

    def test_add_auto_increment_id(self):
        blocks = {}
        al = SgffAttachmentList(blocks)
        al.add(SgffAttachment(name="a.txt", data=b"a"))
        al.add(SgffAttachment(name="b.txt", data=b"b"))
        assert al[0].id == 1
        assert al[1].id == 2

    def test_remove(self):
        blocks = {}
        al = SgffAttachmentList(blocks)
        al.add(SgffAttachment(name="a.txt", data=b"a"))
        al.add(SgffAttachment(name="b.txt", data=b"b"))
        al.remove(0)
        assert len(al) == 1
        assert al[0].name == "b.txt"

    def test_clear(self):
        blocks = {}
        al = SgffAttachmentList(blocks)
        al.add(SgffAttachment(name="a.txt", data=b"a"))
        al.clear()
        assert len(al) == 0
        assert 23 not in blocks

    def test_get_by_name(self):
        blocks = {}
        al = SgffAttachmentList(blocks)
        al.add(SgffAttachment(name="find_me.txt", data=b"here"))
        assert al.get_by_name("find_me.txt") is not None
        assert al.get_by_name("missing.txt") is None

    def test_get_by_id(self):
        blocks = {}
        al = SgffAttachmentList(blocks)
        al.add(SgffAttachment(name="f.txt", data=b"x"))
        assert al.get_by_id(1) is not None
        assert al.get_by_id(99) is None


# =============================================================================
# TestAttachmentRoundtrip
# =============================================================================


class TestAttachmentRoundtrip:
    def _roundtrip(self, sgff):
        """Write and re-read an SgffObject."""
        data = SgffWriter.to_bytes(sgff)
        return SgffReader(BytesIO(data)).read()

    def test_single_file_roundtrip(self):
        sgff = SgffObject.new("ATCG")
        sgff.add_attachment("test.jpg", b"\xff\xd8\xff\xe0" + b"\x00" * 100)

        sgff2 = self._roundtrip(sgff)
        assert sgff2.has_attachments
        assert len(sgff2.attachments) == 1
        att = sgff2.attachments[0]
        assert att.name == "test.jpg"
        assert att.data == b"\xff\xd8\xff\xe0" + b"\x00" * 100

    def test_multiple_files_roundtrip(self):
        sgff = SgffObject.new("ATCG")
        sgff.add_attachment("a.bin", b"alpha")
        sgff.add_attachment("b.bin", b"bravo")

        sgff2 = self._roundtrip(sgff)
        assert len(sgff2.attachments) == 2
        names = {a.name for a in sgff2.attachments}
        assert names == {"a.bin", "b.bin"}

    def test_binary_data_integrity(self):
        raw = bytes(range(256)) * 10
        sgff = SgffObject.new("ATCG")
        sgff.add_attachment("binary.dat", raw)

        sgff2 = self._roundtrip(sgff)
        assert sgff2.attachments[0].data == raw

    def test_double_roundtrip_stability(self):
        sgff = SgffObject.new("ATCG")
        sgff.add_attachment("stable.bin", b"deterministic content")

        bytes1 = SgffWriter.to_bytes(sgff)
        sgff2 = SgffReader(BytesIO(bytes1)).read()
        bytes2 = SgffWriter.to_bytes(sgff2)
        sgff3 = SgffReader(BytesIO(bytes2)).read()
        bytes3 = SgffWriter.to_bytes(sgff3)

        assert bytes2 == bytes3

    def test_empty_attachment_data(self):
        sgff = SgffObject.new("ATCG")
        sgff.add_attachment("empty.txt", b"")

        sgff2 = self._roundtrip(sgff)
        assert len(sgff2.attachments) == 1
        assert sgff2.attachments[0].data == b""


# =============================================================================
# TestAttachmentIntegration
# =============================================================================


class TestAttachmentIntegration:
    def test_has_attachments_false(self):
        sgff = SgffObject.new("ATCG")
        assert not sgff.has_attachments

    def test_has_attachments_true(self):
        sgff = SgffObject.new("ATCG")
        sgff.add_attachment("f.txt", b"x")
        assert sgff.has_attachments

    def test_invalidate_cache(self):
        sgff = SgffObject.new("ATCG")
        sgff.add_attachment("f.txt", b"x")
        _ = sgff.attachments  # populate cache
        sgff.invalidate()
        # Should reload from blocks
        assert len(sgff.attachments) == 1

    def test_chained_builder(self):
        sgff = (
            SgffObject.new("ATCG")
            .add_attachment("a.txt", b"aaa")
            .add_attachment("b.txt", b"bbb")
        )
        assert len(sgff.attachments) == 2

    def test_attachment_with_features(self):
        sgff = (
            SgffObject.new("ATCGATCG")
            .add_feature("gene1", "CDS", 0, 4)
            .add_attachment("img.png", b"\x89PNG")
        )
        data = SgffWriter.to_bytes(sgff)
        sgff2 = SgffReader(BytesIO(data)).read()
        assert sgff2.has_features
        assert sgff2.has_attachments
        assert len(sgff2.features) == 1
        assert len(sgff2.attachments) == 1
