"""
Tests for parser functions in sgffp.parsers
"""

import struct
import lzma
from io import BytesIO

import pytest

from sgffp.parsers import (
    parse_blocks,
    read_header,
    octet_to_dna,
    parse_sequence,
    parse_compressed_dna,
    parse_xml,
    parse_lzma_xml,
    parse_lzma_nested,
    parse_features,
    parse_ztr,
    parse_history_node,
    SCHEME,
    ZTR_MAGIC,
)


# =============================================================================
# Helper Functions Tests
# =============================================================================


class TestReadHeader:
    def test_read_header_valid(self):
        """Read valid TLV header"""
        data = bytes([10]) + struct.pack(">I", 100)
        stream = BytesIO(data)
        block_type, block_length = read_header(stream)
        assert block_type == 10
        assert block_length == 100

    def test_read_header_empty(self):
        """Empty stream returns None, None"""
        stream = BytesIO(b"")
        block_type, block_length = read_header(stream)
        assert block_type is None
        assert block_length is None


class TestOctetToDna:
    def test_octet_to_dna_gatc(self):
        """2-bit encoding: G=0, A=1, T=2, C=3"""
        # 0b00011011 = G(00), A(01), T(10), C(11) = GATC
        raw = bytes([0b00011011])
        result = octet_to_dna(raw, 4)
        assert result == b"GATC"

    def test_octet_to_dna_partial(self):
        """Decode partial byte (fewer bases than byte holds)"""
        raw = bytes([0b00011011])
        result = octet_to_dna(raw, 2)
        assert result == b"GA"

    def test_octet_to_dna_multiple_bytes(self):
        """Decode multiple bytes"""
        # GATC GATC
        raw = bytes([0b00011011, 0b00011011])
        result = octet_to_dna(raw, 8)
        assert result == b"GATCGATC"

    def test_octet_to_dna_all_g(self):
        """All zeros = all G"""
        raw = bytes([0b00000000])
        result = octet_to_dna(raw, 4)
        assert result == b"GGGG"

    def test_octet_to_dna_all_c(self):
        """All 3s = all C"""
        raw = bytes([0b11111111])
        result = octet_to_dna(raw, 4)
        assert result == b"CCCC"


# =============================================================================
# Sequence Parser Tests
# =============================================================================


class TestParseSequence:
    def test_parse_sequence_basic(self):
        """Parse basic DNA sequence"""
        # props=0, sequence="ATCG"
        data = bytes([0]) + b"ATCG"
        result = parse_sequence(data)
        assert result["sequence"] == "ATCG"
        assert result["length"] == 4

    def test_parse_sequence_circular(self):
        """Circular topology flag (0x01)"""
        data = bytes([0x01]) + b"ATCG"
        result = parse_sequence(data)
        assert result["topology"] == "circular"

    def test_parse_sequence_linear(self):
        """Linear topology (no flag)"""
        data = bytes([0x00]) + b"ATCG"
        result = parse_sequence(data)
        assert result["topology"] == "linear"

    def test_parse_sequence_double_stranded(self):
        """Double stranded flag (0x02)"""
        data = bytes([0x02]) + b"ATCG"
        result = parse_sequence(data)
        assert result["strandedness"] == "double"

    def test_parse_sequence_single_stranded(self):
        """Single stranded (no flag)"""
        data = bytes([0x00]) + b"ATCG"
        result = parse_sequence(data)
        assert result["strandedness"] == "single"

    def test_parse_sequence_dam_methylated(self):
        """DAM methylation flag (0x04)"""
        data = bytes([0x04]) + b"ATCG"
        result = parse_sequence(data)
        assert result["dam_methylated"] is True

    def test_parse_sequence_dcm_methylated(self):
        """DCM methylation flag (0x08)"""
        data = bytes([0x08]) + b"ATCG"
        result = parse_sequence(data)
        assert result["dcm_methylated"] is True

    def test_parse_sequence_ecoki_methylated(self):
        """EcoKI methylation flag (0x10)"""
        data = bytes([0x10]) + b"ATCG"
        result = parse_sequence(data)
        assert result["ecoki_methylated"] is True

    def test_parse_sequence_all_flags(self):
        """All flags set"""
        # 0x1F = 0b00011111 = all 5 flags
        data = bytes([0x1F]) + b"ATCG"
        result = parse_sequence(data)
        assert result["topology"] == "circular"
        assert result["strandedness"] == "double"
        assert result["dam_methylated"] is True
        assert result["dcm_methylated"] is True
        assert result["ecoki_methylated"] is True

    def test_parse_sequence_empty(self):
        """Empty sequence"""
        data = bytes([0])
        result = parse_sequence(data)
        assert result["sequence"] == ""
        assert result["length"] == 0


class TestParseCompressedDna:
    def test_parse_compressed_dna(self):
        """Parse 2-bit compressed DNA"""
        # Build compressed block:
        # 4 bytes: compressed_length
        # 4 bytes: uncompressed_length (base count)
        # 14 bytes: mystery bytes
        # N bytes: 2-bit encoded sequence

        sequence = "GATCGATC"  # 8 bases
        base_count = 8

        # Encode sequence: GATC = 0b00011011
        encoded = bytes([0b00011011, 0b00011011])

        compressed_length = 4 + 14 + len(encoded)  # uncompressed_len + mystery + encoded

        data = (
            struct.pack(">I", compressed_length)
            + struct.pack(">I", base_count)
            + (b"\x00" * 14)  # mystery bytes (14 bytes, but only 10 preserved)
            + encoded
        )

        result = parse_compressed_dna(data)
        assert result["sequence"] == "GATCGATC"
        assert result["length"] == 8

    def test_parse_compressed_dna_mystery_bytes(self):
        """Mystery bytes are preserved"""
        mystery = b"\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a"
        sequence = "GATC"
        base_count = 4
        encoded = bytes([0b00011011])

        compressed_length = 4 + 14 + len(encoded)

        data = (
            struct.pack(">I", compressed_length)
            + struct.pack(">I", base_count)
            + mystery
            + (b"\x00" * 4)  # remaining mystery bytes
            + encoded
        )

        result = parse_compressed_dna(data)
        assert result["mystery"] == mystery  # First 10 bytes preserved


# =============================================================================
# XML Parser Tests
# =============================================================================


class TestParseXml:
    def test_parse_xml_valid(self):
        """Valid XML parsed to dict"""
        xml = b"<Root><Child>value</Child></Root>"
        result = parse_xml(xml)
        assert result is not None
        assert "Root" in result
        assert result["Root"]["Child"] == "value"

    def test_parse_xml_with_attributes(self):
        """XML with attributes"""
        xml = b'<Root attr="test"><Child>value</Child></Root>'
        result = parse_xml(xml)
        assert result["Root"]["attr"] == "test"

    def test_parse_xml_malformed(self):
        """Malformed XML returns None"""
        xml = b"<Root><Unclosed>"
        result = parse_xml(xml)
        assert result is None

    def test_parse_xml_empty(self):
        """Empty data returns None"""
        result = parse_xml(b"")
        assert result is None

    def test_parse_xml_not_xml(self):
        """Non-XML data returns None"""
        result = parse_xml(b"this is not xml")
        assert result is None


class TestParseLzmaXml:
    def test_parse_lzma_xml_valid(self):
        """LZMA-compressed XML"""
        xml = b"<Root><Item>test</Item></Root>"
        compressed = lzma.compress(xml)
        result = parse_lzma_xml(compressed)
        assert result is not None
        assert result["Root"]["Item"] == "test"

    def test_parse_lzma_xml_invalid(self):
        """Invalid LZMA data returns None"""
        result = parse_lzma_xml(b"not lzma data")
        assert result is None


class TestParseLzmaNested:
    def test_parse_lzma_nested(self):
        """LZMA with nested TLV blocks"""
        # Create a simple TLV block: type 0 (sequence), length 5, data
        tlv_data = bytes([0]) + struct.pack(">I", 5) + bytes([0]) + b"ATCG"

        compressed = lzma.compress(tlv_data)
        result = parse_lzma_nested(compressed)

        assert result is not None
        assert 0 in result  # Sequence block type
        assert result[0][0]["sequence"] == "ATCG"

    def test_parse_lzma_nested_invalid(self):
        """Invalid data returns None"""
        result = parse_lzma_nested(b"invalid")
        assert result is None


# =============================================================================
# Feature Parser Tests
# =============================================================================


class TestParseFeatures:
    def test_parse_features_single(self):
        """Parse single feature"""
        xml = b"""<Features>
            <Feature name="test" type="gene" directionality="1">
                <Segment range="1-100" color="#FF0000"/>
            </Feature>
        </Features>"""
        result = parse_features(xml)
        assert result is not None
        assert "features" in result
        assert len(result["features"]) == 1
        assert result["features"][0]["name"] == "test"
        assert result["features"][0]["type"] == "gene"

    def test_parse_features_multiple(self):
        """Parse multiple features"""
        xml = b"""<Features>
            <Feature name="feat1" type="gene" directionality="1">
                <Segment range="1-50"/>
            </Feature>
            <Feature name="feat2" type="CDS" directionality="2">
                <Segment range="51-100"/>
            </Feature>
        </Features>"""
        result = parse_features(xml)
        assert len(result["features"]) == 2

    def test_parse_features_strands(self):
        """All strand types: 0=., 1=+, 2=-, 3==="""
        strands = [("0", "."), ("1", "+"), ("2", "-"), ("3", "=")]
        for dir_val, expected in strands:
            xml = f"""<Features>
                <Feature name="test" type="gene" directionality="{dir_val}">
                    <Segment range="1-10"/>
                </Feature>
            </Features>""".encode()
            result = parse_features(xml)
            assert result["features"][0]["strand"] == expected

    def test_parse_features_qualifiers(self):
        """Qualifier extraction"""
        xml = b"""<Features>
            <Feature name="test" type="gene" directionality="1">
                <Segment range="1-100"/>
                <Q name="note"><V text="A note"/></Q>
                <Q name="product"><V text="Test protein"/></Q>
            </Feature>
        </Features>"""
        result = parse_features(xml)
        quals = result["features"][0]["qualifiers"]
        assert "note" in quals
        assert quals["note"] == "A note"
        assert quals["product"] == "Test protein"

    def test_parse_features_segments(self):
        """Multi-segment features"""
        xml = b"""<Features>
            <Feature name="test" type="gene" directionality="1">
                <Segment range="1-50"/>
                <Segment range="100-150"/>
            </Feature>
        </Features>"""
        result = parse_features(xml)
        feat = result["features"][0]
        assert len(feat["segments"]) == 2
        assert feat["start"] == 0  # 1-1 = 0 (0-indexed)
        assert feat["end"] == 150

    def test_parse_features_colors(self):
        """Color attribute extracted"""
        xml = b"""<Features>
            <Feature name="test" type="gene" directionality="1">
                <Segment range="1-100" color="#FF0000"/>
            </Feature>
        </Features>"""
        result = parse_features(xml)
        assert result["features"][0]["color"] == "#FF0000"

    def test_parse_features_no_features(self):
        """Empty features element returns empty list"""
        xml = b"<Features></Features>"
        result = parse_features(xml)
        assert result == {"features": []}

    def test_parse_features_invalid_xml(self):
        """Invalid XML returns None"""
        result = parse_features(b"not xml")
        assert result is None


# =============================================================================
# ZTR Parser Tests
# =============================================================================


class TestParseZtr:
    def test_parse_ztr_valid_magic(self):
        """Valid ZTR magic bytes accepted"""
        # Minimal ZTR: magic + version (2 bytes)
        data = ZTR_MAGIC + b"\x01\x00"
        result = parse_ztr(data)
        # Should return empty dict (no chunks) or None if too short
        assert result is None or result == {}

    def test_parse_ztr_invalid_magic(self):
        """Wrong magic bytes returns None"""
        data = b"NOTZTR\r\n\x1a\n" + b"\x01\x00"
        result = parse_ztr(data)
        assert result is None

    def test_parse_ztr_too_short(self):
        """Data too short returns None"""
        result = parse_ztr(b"short")
        assert result is None


# =============================================================================
# History Parser Tests
# =============================================================================


class TestParseHistoryNode:
    def test_parse_history_node_basic(self):
        """Basic history node parsing"""
        # node_index (4 bytes) + sequence_type (1 byte)
        data = struct.pack(">I", 1) + bytes([0])  # type 0 = uncompressed

        # For type 0: seq_length (4 bytes) + sequence
        seq = b"ATCG"
        data += struct.pack(">I", len(seq)) + seq

        result = parse_history_node(data)
        assert result["node_index"] == 1
        assert result["sequence_type"] == 0
        assert result["sequence"] == "ATCG"
        assert result["length"] == 4


# =============================================================================
# Scheme Tests
# =============================================================================


class TestScheme:
    def test_scheme_has_expected_types(self):
        """SCHEME contains expected block types"""
        expected = [0, 1, 5, 6, 7, 8, 10, 11, 16, 17, 18, 21, 29, 30, 32]
        for block_type in expected:
            assert block_type in SCHEME

    def test_scheme_sequence_types(self):
        """Sequence types (0, 21, 32) use parse_sequence"""
        for block_type in [0, 21, 32]:
            _, parser = SCHEME[block_type]
            assert parser == parse_sequence

    def test_scheme_compressed_dna(self):
        """Type 1 uses parse_compressed_dna"""
        _, parser = SCHEME[1]
        assert parser == parse_compressed_dna

    def test_scheme_features(self):
        """Type 10 uses parse_features"""
        _, parser = SCHEME[10]
        assert parser == parse_features


# =============================================================================
# Parse Blocks Integration Test
# =============================================================================


class TestParseBlocks:
    def test_parse_blocks_single(self):
        """Parse single TLV block"""
        # Block type 0 (sequence), length, data
        block_data = bytes([0]) + b"ATCG"  # props + sequence
        stream = BytesIO(
            bytes([0])  # type
            + struct.pack(">I", len(block_data))
            + block_data
        )

        result = parse_blocks(stream)
        assert 0 in result
        assert result[0][0]["sequence"] == "ATCG"

    def test_parse_blocks_multiple(self):
        """Parse multiple TLV blocks"""
        buf = BytesIO()

        # Block 0: sequence
        seq_data = bytes([0]) + b"ATCG"
        buf.write(bytes([0]))
        buf.write(struct.pack(">I", len(seq_data)))
        buf.write(seq_data)

        # Block 6: notes (XML)
        xml_data = b"<Notes><Note>test</Note></Notes>"
        buf.write(bytes([6]))
        buf.write(struct.pack(">I", len(xml_data)))
        buf.write(xml_data)

        buf.seek(0)
        result = parse_blocks(buf)

        assert 0 in result
        assert 6 in result

    def test_parse_blocks_unknown_skipped(self):
        """Unknown block types are skipped"""
        buf = BytesIO()

        # Unknown block type 99
        buf.write(bytes([99]))
        buf.write(struct.pack(">I", 10))
        buf.write(b"0123456789")

        # Known block type 0
        seq_data = bytes([0]) + b"ATCG"
        buf.write(bytes([0]))
        buf.write(struct.pack(">I", len(seq_data)))
        buf.write(seq_data)

        buf.seek(0)
        result = parse_blocks(buf)

        assert 99 not in result
        assert 0 in result

    def test_parse_blocks_empty_stream(self):
        """Empty stream returns empty dict"""
        result = parse_blocks(BytesIO(b""))
        assert result == {}
