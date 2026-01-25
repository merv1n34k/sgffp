"""
Tests for SgffReader class
"""

import struct
from io import BytesIO
from pathlib import Path

import pytest

from sgffp.reader import SgffReader
from sgffp.internal import SgffObject, Cookie


def make_minimal_sgff(sequence=b"ATCG", cookie_type=1, export_ver=16, import_ver=8):
    """Create minimal valid SnapGene file bytes"""
    buf = BytesIO()

    # Header
    buf.write(b"\t")  # Magic byte
    buf.write(struct.pack(">I", 14))  # Header length
    buf.write(b"SnapGene")  # Title

    # Cookie
    buf.write(struct.pack(">H", cookie_type))
    buf.write(struct.pack(">H", export_ver))
    buf.write(struct.pack(">H", import_ver))

    # Block 0: sequence
    block_data = bytes([0]) + sequence  # props + sequence
    buf.write(bytes([0]))  # Block type
    buf.write(struct.pack(">I", len(block_data)))
    buf.write(block_data)

    return buf.getvalue()


# =============================================================================
# Constructor Tests
# =============================================================================


class TestSgffReaderConstructor:
    def test_from_file_path_str(self, test_dna):
        """Read using string path"""
        sgff = SgffReader.from_file(str(test_dna))
        assert isinstance(sgff, SgffObject)

    def test_from_file_path_pathlib(self, test_dna):
        """Read using Path object"""
        sgff = SgffReader.from_file(test_dna)
        assert isinstance(sgff, SgffObject)

    def test_from_bytes(self):
        """Read from bytes directly"""
        data = make_minimal_sgff()
        sgff = SgffReader.from_bytes(data)
        assert isinstance(sgff, SgffObject)

    def test_from_bytesio(self):
        """Read from BytesIO stream"""
        data = make_minimal_sgff()
        stream = BytesIO(data)
        reader = SgffReader(stream)
        sgff = reader.read()
        assert isinstance(sgff, SgffObject)


# =============================================================================
# Error Handling Tests
# =============================================================================


class TestSgffReaderErrors:
    def test_invalid_magic_byte(self):
        """Wrong magic byte raises ValueError"""
        # Replace magic byte \t with \x00
        data = b"\x00" + make_minimal_sgff()[1:]
        with pytest.raises(ValueError, match="wrong magic byte"):
            SgffReader.from_bytes(data)

    def test_invalid_header_length(self):
        """Wrong header length raises ValueError"""
        buf = BytesIO()
        buf.write(b"\t")
        buf.write(struct.pack(">I", 10))  # Wrong length (should be 14)
        buf.write(b"SnapGene")
        buf.write(b"\x00" * 6)  # Cookie

        with pytest.raises(ValueError, match="wrong header"):
            SgffReader.from_bytes(buf.getvalue())

    def test_invalid_title(self):
        """Missing 'SnapGene' raises ValueError"""
        buf = BytesIO()
        buf.write(b"\t")
        buf.write(struct.pack(">I", 14))
        buf.write(b"NotSnap!")  # Wrong title
        buf.write(b"\x00" * 6)

        with pytest.raises(ValueError, match="wrong header"):
            SgffReader.from_bytes(buf.getvalue())

    def test_truncated_file(self):
        """Incomplete file raises error"""
        # Only magic byte, no header
        with pytest.raises(Exception):
            SgffReader.from_bytes(b"\t")

    def test_empty_file(self):
        """Empty file raises error"""
        with pytest.raises(Exception):
            SgffReader.from_bytes(b"")


# =============================================================================
# Cookie Parsing Tests
# =============================================================================


class TestSgffReaderCookie:
    def test_cookie_parsing(self):
        """Cookie fields correctly extracted"""
        data = make_minimal_sgff(cookie_type=1, export_ver=16, import_ver=8)
        sgff = SgffReader.from_bytes(data)

        assert sgff.cookie.type_of_sequence == 1
        assert sgff.cookie.export_version == 16
        assert sgff.cookie.import_version == 8

    def test_cookie_different_values(self):
        """Cookie with different values"""
        data = make_minimal_sgff(cookie_type=2, export_ver=32, import_ver=16)
        sgff = SgffReader.from_bytes(data)

        assert sgff.cookie.type_of_sequence == 2
        assert sgff.cookie.export_version == 32
        assert sgff.cookie.import_version == 16


# =============================================================================
# Block Parsing Tests
# =============================================================================


class TestSgffReaderBlocks:
    def test_blocks_populated(self):
        """Blocks dict contains parsed data"""
        data = make_minimal_sgff(sequence=b"ATCG")
        sgff = SgffReader.from_bytes(data)

        assert 0 in sgff.blocks
        assert sgff.blocks[0][0]["sequence"] == "ATCG"

    def test_blocks_sequence_properties(self):
        """Sequence block properties parsed"""
        data = make_minimal_sgff(sequence=b"GATC")
        sgff = SgffReader.from_bytes(data)

        seq_block = sgff.blocks[0][0]
        assert seq_block["length"] == 4
        assert seq_block["topology"] == "linear"
        assert seq_block["strandedness"] == "single"


# =============================================================================
# Real File Tests
# =============================================================================


class TestSgffReaderRealFiles:
    def test_read_test_dna(self, test_dna):
        """Read test.dna successfully"""
        sgff = SgffReader.from_file(test_dna)
        assert isinstance(sgff, SgffObject)
        assert sgff.cookie is not None

    def test_read_test2_dna(self, test2_dna):
        """Read test2.dna successfully"""
        sgff = SgffReader.from_file(test2_dna)
        assert isinstance(sgff, SgffObject)

    def test_read_test3_dna(self, test3_dna):
        """Read test3.dna successfully"""
        sgff = SgffReader.from_file(test3_dna)
        assert isinstance(sgff, SgffObject)

    def test_read_test_rna(self, test_rna):
        """Read test.rna successfully"""
        sgff = SgffReader.from_file(test_rna)
        assert isinstance(sgff, SgffObject)

    def test_read_test_prot(self, test_prot):
        """Read test.prot successfully"""
        sgff = SgffReader.from_file(test_prot)
        assert isinstance(sgff, SgffObject)

    def test_real_file_has_sequence(self, test_dna):
        """Real file contains sequence block"""
        sgff = SgffReader.from_file(test_dna)
        # DNA files should have block 0 or 1
        has_sequence = 0 in sgff.blocks or 1 in sgff.blocks
        assert has_sequence

    def test_real_file_sequence_length(self, test_dna):
        """Real file has non-empty sequence"""
        sgff = SgffReader.from_file(test_dna)
        if 0 in sgff.blocks:
            assert sgff.blocks[0][0]["length"] > 0
        elif 1 in sgff.blocks:
            assert sgff.blocks[1][0]["length"] > 0
