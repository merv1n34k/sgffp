"""
Tests for SgffWriter class
"""

import struct
import tempfile
from io import BytesIO
from pathlib import Path

import pytest

from sgffp.writer import SgffWriter
from sgffp.reader import SgffReader
from sgffp.internal import SgffObject, Cookie


@pytest.fixture
def sample_cookie():
    return Cookie(type_of_sequence=1, export_version=16, import_version=8)


@pytest.fixture
def sample_sgff(sample_cookie):
    obj = SgffObject(cookie=sample_cookie)
    obj.blocks = {
        0: [
            {
                "sequence": "ATCGATCG",
                "length": 8,
                "topology": "linear",
                "strandedness": "single",
                "dam_methylated": False,
                "dcm_methylated": False,
                "ecoki_methylated": False,
            }
        ]
    }
    return obj


# =============================================================================
# Constructor Tests
# =============================================================================


class TestSgffWriterConstructor:
    def test_to_file(self, sample_sgff):
        """Write to file path"""
        with tempfile.NamedTemporaryFile(suffix=".dna", delete=False) as f:
            filepath = f.name

        try:
            SgffWriter.to_file(sample_sgff, filepath)
            assert Path(filepath).exists()
            assert Path(filepath).stat().st_size > 0
        finally:
            Path(filepath).unlink(missing_ok=True)

    def test_to_bytes(self, sample_sgff):
        """Get bytes output"""
        data = SgffWriter.to_bytes(sample_sgff)
        assert isinstance(data, bytes)
        assert len(data) > 0

    def test_to_bytesio(self, sample_sgff):
        """Write to BytesIO stream"""
        stream = BytesIO()
        writer = SgffWriter(stream)
        writer.write(sample_sgff)
        assert stream.tell() > 0


# =============================================================================
# Header Tests
# =============================================================================


class TestSgffWriterHeader:
    def test_write_header_magic_byte(self, sample_sgff):
        """Correct magic byte (\t)"""
        data = SgffWriter.to_bytes(sample_sgff)
        assert data[0] == ord("\t")

    def test_write_header_length(self, sample_sgff):
        """Header length is 14"""
        data = SgffWriter.to_bytes(sample_sgff)
        length = struct.unpack(">I", data[1:5])[0]
        assert length == 14

    def test_write_header_title(self, sample_sgff):
        """Header contains 'SnapGene'"""
        data = SgffWriter.to_bytes(sample_sgff)
        title = data[5:13]
        assert title == b"SnapGene"


# =============================================================================
# Cookie Tests
# =============================================================================


class TestSgffWriterCookie:
    def test_write_cookie(self, sample_sgff):
        """Cookie serialized correctly"""
        data = SgffWriter.to_bytes(sample_sgff)

        # Cookie starts at byte 13 (after magic + length + title)
        cookie_start = 13
        type_of_seq = struct.unpack(">H", data[cookie_start : cookie_start + 2])[0]
        export_ver = struct.unpack(">H", data[cookie_start + 2 : cookie_start + 4])[0]
        import_ver = struct.unpack(">H", data[cookie_start + 4 : cookie_start + 6])[0]

        assert type_of_seq == 1
        assert export_ver == 16
        assert import_ver == 8


# =============================================================================
# Sequence Block Tests
# =============================================================================


class TestSgffWriterSequence:
    def test_write_sequence_block_0(self, sample_cookie):
        """DNA sequence with properties byte"""
        obj = SgffObject(cookie=sample_cookie)
        obj.blocks = {
            0: [
                {
                    "sequence": "ATCG",
                    "topology": "linear",
                    "strandedness": "single",
                    "dam_methylated": False,
                    "dcm_methylated": False,
                    "ecoki_methylated": False,
                }
            ]
        }

        data = SgffWriter.to_bytes(obj)

        # Read it back to verify
        sgff = SgffReader.from_bytes(data)
        assert sgff.blocks[0][0]["sequence"] == "ATCG"

    def test_write_sequence_block_21(self, sample_cookie):
        """Protein sequence"""
        obj = SgffObject(cookie=sample_cookie)
        obj.blocks = {
            21: [
                {
                    "sequence": "MKTL",
                    "topology": "linear",
                    "strandedness": "single",
                    "dam_methylated": False,
                    "dcm_methylated": False,
                    "ecoki_methylated": False,
                }
            ]
        }

        data = SgffWriter.to_bytes(obj)
        sgff = SgffReader.from_bytes(data)
        assert 21 in sgff.blocks

    def test_write_sequence_block_32(self, sample_cookie):
        """RNA sequence"""
        obj = SgffObject(cookie=sample_cookie)
        obj.blocks = {
            32: [
                {
                    "sequence": "AUCG",
                    "topology": "linear",
                    "strandedness": "single",
                    "dam_methylated": False,
                    "dcm_methylated": False,
                    "ecoki_methylated": False,
                }
            ]
        }

        data = SgffWriter.to_bytes(obj)
        sgff = SgffReader.from_bytes(data)
        assert 32 in sgff.blocks

    def test_topology_flags_circular(self, sample_cookie):
        """Circular topology preserved"""
        obj = SgffObject(cookie=sample_cookie)
        obj.blocks = {
            0: [
                {
                    "sequence": "ATCG",
                    "topology": "circular",
                    "strandedness": "single",
                    "dam_methylated": False,
                    "dcm_methylated": False,
                    "ecoki_methylated": False,
                }
            ]
        }

        data = SgffWriter.to_bytes(obj)
        sgff = SgffReader.from_bytes(data)
        assert sgff.blocks[0][0]["topology"] == "circular"

    def test_topology_flags_linear(self, sample_cookie):
        """Linear topology preserved"""
        obj = SgffObject(cookie=sample_cookie)
        obj.blocks = {
            0: [
                {
                    "sequence": "ATCG",
                    "topology": "linear",
                    "strandedness": "double",
                    "dam_methylated": False,
                    "dcm_methylated": False,
                    "ecoki_methylated": False,
                }
            ]
        }

        data = SgffWriter.to_bytes(obj)
        sgff = SgffReader.from_bytes(data)
        assert sgff.blocks[0][0]["topology"] == "linear"

    def test_methylation_flags(self, sample_cookie):
        """DAM/DCM/EcoKI methylation flags preserved"""
        obj = SgffObject(cookie=sample_cookie)
        obj.blocks = {
            0: [
                {
                    "sequence": "ATCG",
                    "topology": "linear",
                    "strandedness": "single",
                    "dam_methylated": True,
                    "dcm_methylated": True,
                    "ecoki_methylated": True,
                }
            ]
        }

        data = SgffWriter.to_bytes(obj)
        sgff = SgffReader.from_bytes(data)

        seq = sgff.blocks[0][0]
        assert seq["dam_methylated"] is True
        assert seq["dcm_methylated"] is True
        assert seq["ecoki_methylated"] is True


# =============================================================================
# Block Ordering Tests
# =============================================================================


class TestSgffWriterOrdering:
    def test_block_ordering(self, sample_cookie):
        """Blocks sorted by type ID"""
        obj = SgffObject(cookie=sample_cookie)
        obj.blocks = {
            10: [{"features": []}],
            0: [
                {
                    "sequence": "ATCG",
                    "topology": "linear",
                    "strandedness": "single",
                    "dam_methylated": False,
                    "dcm_methylated": False,
                    "ecoki_methylated": False,
                }
            ],
            6: [{"Notes": {}}],
        }

        data = SgffWriter.to_bytes(obj)

        # Parse block types in order
        stream = BytesIO(data)
        stream.read(19)  # Skip header + cookie

        block_types = []
        while True:
            type_byte = stream.read(1)
            if not type_byte:
                break
            block_types.append(type_byte[0])
            length = struct.unpack(">I", stream.read(4))[0]
            stream.read(length)

        # Should be sorted: 0, 6, 10
        assert block_types == sorted(block_types)


# =============================================================================
# DNA Encoding Tests
# =============================================================================


class TestDnaToOctet:
    def test_dna_to_octet_gatc(self):
        """G=0, A=1, T=2, C=3 encoding"""
        writer = SgffWriter(BytesIO())
        result = writer._dna_to_octet("GATC")
        # G=00, A=01, T=10, C=11 -> 0b00011011 = 27
        assert result == bytes([0b00011011])

    def test_dna_to_octet_all_g(self):
        """All G encodes to zeros"""
        writer = SgffWriter(BytesIO())
        result = writer._dna_to_octet("GGGG")
        assert result == bytes([0b00000000])

    def test_dna_to_octet_all_c(self):
        """All C encodes to all 1s"""
        writer = SgffWriter(BytesIO())
        result = writer._dna_to_octet("CCCC")
        assert result == bytes([0b11111111])

    def test_dna_to_octet_padding(self):
        """Odd-length sequences handled"""
        writer = SgffWriter(BytesIO())
        # "AT" = A(01), T(10) = 0b01100000 (padded with zeros)
        result = writer._dna_to_octet("AT")
        assert result == bytes([0b01100000])

    def test_dna_to_octet_unknown_base(self):
        """Unknown bases default to G (0)"""
        writer = SgffWriter(BytesIO())
        result = writer._dna_to_octet("NXYZ")
        # All unknown = all zeros = 0b00000000
        assert result == bytes([0b00000000])

    def test_dna_to_octet_lowercase(self):
        """Lowercase bases handled"""
        writer = SgffWriter(BytesIO())
        result = writer._dna_to_octet("gatc")
        assert result == bytes([0b00011011])


# =============================================================================
# Compressed DNA Tests
# =============================================================================


class TestSgffWriterCompressedDna:
    def test_write_compressed_dna(self, sample_cookie):
        """2-bit encoding (block 1)"""
        obj = SgffObject(cookie=sample_cookie)
        obj.blocks = {
            1: [
                {
                    "sequence": "GATCGATC",
                    "length": 8,
                    "mystery": b"\x00" * 14,  # Mystery bytes should be 14 bytes
                }
            ]
        }

        data = SgffWriter.to_bytes(obj)
        sgff = SgffReader.from_bytes(data)

        assert 1 in sgff.blocks
        assert sgff.blocks[1][0]["sequence"] == "GATCGATC"


# =============================================================================
# XML Block Tests
# =============================================================================


class TestSgffWriterXml:
    def test_write_xml_blocks(self, sample_cookie):
        """XML blocks (5, 6, 8, 17)"""
        obj = SgffObject(cookie=sample_cookie)
        obj.blocks = {
            6: [{"Notes": {"Note": "test note"}}],
        }

        data = SgffWriter.to_bytes(obj)
        sgff = SgffReader.from_bytes(data)

        assert 6 in sgff.blocks


# =============================================================================
# Features Block Tests
# =============================================================================


class TestSgffWriterFeatures:
    def test_write_features_block_10(self, sample_cookie):
        """Features with qualifiers"""
        obj = SgffObject(cookie=sample_cookie)
        obj.blocks = {
            0: [
                {
                    "sequence": "ATCGATCG",
                    "topology": "linear",
                    "strandedness": "single",
                    "dam_methylated": False,
                    "dcm_methylated": False,
                    "ecoki_methylated": False,
                }
            ],
            10: [
                {
                    "features": [
                        {
                            "name": "test_gene",
                            "type": "gene",
                            "strand": "+",
                            "segments": [{"@range": "1-8"}],
                            "qualifiers": {"note": "A test gene"},
                        }
                    ]
                }
            ],
        }

        data = SgffWriter.to_bytes(obj)
        sgff = SgffReader.from_bytes(data)

        assert 10 in sgff.blocks


# =============================================================================
# Real File Roundtrip Tests
# =============================================================================


class TestSgffWriterRealFiles:
    def test_write_read_test_dna(self, test_dna):
        """Read test.dna, write, read again (supported blocks only)"""
        original = SgffReader.from_file(test_dna)

        # Filter out blocks that can't be round-tripped (history nodes, etc.)
        skip_blocks = {7, 11, 29, 30}
        filtered = SgffObject(cookie=original.cookie)
        for block_type, items in original.blocks.items():
            if block_type not in skip_blocks:
                filtered.blocks[block_type] = items

        data = SgffWriter.to_bytes(filtered)

        # Should be valid SnapGene format
        assert data[:1] == b"\t"
        assert data[5:13] == b"SnapGene"

        # Should be readable
        restored = SgffReader.from_bytes(data)
        assert restored is not None
