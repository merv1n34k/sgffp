"""
Tests for the ZTR chunk codec and trace round-tripping (issue #23)

ZTR chunk data is filtered through a chain of encodings. Decoding only the
outermost layer leaves the rest of the chain to be misread as values, which is
what corrupted traces on write. These tests pin both halves: the decoded values
must be real, and re-encoding must reproduce the original bytes.
"""

import struct

import pytest

from sgffp import ztr
from sgffp.parsers import ZTR_MAGIC, parse_ztr
from sgffp.reader import SgffReader
from sgffp.writer import SgffWriter


def iter_chunks(blob):
    """Yield (type, meta, encoded_data) for each chunk in a ZTR blob."""
    offset = 10
    while offset + 8 <= len(blob):
        chunk_type = blob[offset : offset + 4].decode("ascii").strip()
        meta_len = struct.unpack(">I", blob[offset + 4 : offset + 8])[0]
        meta = blob[offset + 8 : offset + 8 + meta_len]
        offset += 8 + meta_len
        data_len = struct.unpack(">I", blob[offset : offset + 4])[0]
        offset += 4
        yield chunk_type, meta, blob[offset : offset + data_len]
        offset += data_len


def collect_ztr_blobs(path):
    """Every raw ZTR blob the reader reaches, at any nesting depth."""
    import sgffp.parsers as parsers

    blobs = []
    original = parsers.SCHEME[18]

    def capture(data):
        blobs.append(data)
        return original(data)

    parsers.SCHEME[18] = capture
    try:
        SgffReader.from_file(path)
    finally:
        parsers.SCHEME[18] = original
    return blobs


def block_16_bytes(raw):
    """Extract every top-level block 16 payload from raw file bytes."""
    out = []
    offset = 19
    while offset + 5 <= len(raw):
        block_type = raw[offset]
        length = struct.unpack(">I", raw[offset + 1 : offset + 5])[0]
        if block_type == 16:
            out.append(raw[offset + 5 : offset + 5 + length])
        offset += 5 + length
    return out


# =============================================================================
# Codec unit tests
# =============================================================================


class TestZtrFilters:
    """Each filter must invert exactly, for arbitrary payloads."""

    # 16TO8/32TO8 step through 2- and 4-byte words, so the payload they see is
    # always word-aligned in practice: a delta filter pads its header to the
    # alignment it needs. Keep the test payload aligned for the same reason.
    PAYLOAD = bytes([0, 0, 0, 0]) + bytes(range(256)) * 3 + b"\x00" * 40
    assert len(PAYLOAD) % 4 == 0

    @pytest.mark.parametrize("level", [1, 2, 3])
    def test_delta1_roundtrip(self, level):
        chain = [{"format": ztr.DELTA1, "level": level}]
        assert ztr.unpack(ztr.pack(self.PAYLOAD, chain))[0] == self.PAYLOAD

    @pytest.mark.parametrize("level", [1, 2, 3])
    def test_delta2_roundtrip(self, level):
        chain = [{"format": ztr.DELTA2, "level": level}]
        assert ztr.unpack(ztr.pack(self.PAYLOAD, chain))[0] == self.PAYLOAD

    @pytest.mark.parametrize("level", [1, 2, 3])
    def test_delta4_roundtrip(self, level):
        chain = [{"format": ztr.DELTA4, "level": level}]
        assert ztr.unpack(ztr.pack(self.PAYLOAD, chain))[0] == self.PAYLOAD

    def test_zlib_roundtrip(self):
        chain = [{"format": ztr.ZLIB}]
        assert ztr.unpack(ztr.pack(self.PAYLOAD, chain))[0] == self.PAYLOAD

    def test_rle_roundtrip(self):
        chain = [{"format": ztr.RLE, "guard": 0xAA}]
        assert ztr.unpack(ztr.pack(self.PAYLOAD, chain))[0] == self.PAYLOAD

    def test_rle_escapes_literal_guard_byte(self):
        """A guard byte in the input is escaped as 'guard 0', not a run."""
        payload = b"\x00\xaa\xaa\x01\xaa"
        chain = [{"format": ztr.RLE, "guard": 0xAA}]
        assert ztr.unpack(ztr.pack(payload, chain))[0] == payload

    def test_rle_long_run(self):
        """Runs longer than 255 bytes split into multiple encoded runs."""
        payload = b"\x00" + b"\x7f" * 700
        chain = [{"format": ztr.RLE, "guard": 0xAA}]
        assert ztr.unpack(ztr.pack(payload, chain))[0] == payload

    def test_16to8_roundtrip(self):
        chain = [{"format": ztr.SHRINK_16TO8}]
        assert ztr.unpack(ztr.pack(self.PAYLOAD, chain))[0] == self.PAYLOAD

    def test_32to8_roundtrip(self):
        chain = [{"format": ztr.SHRINK_32TO8}]
        assert ztr.unpack(ztr.pack(self.PAYLOAD, chain))[0] == self.PAYLOAD

    def test_follow1_roundtrip(self):
        chain = [{"format": ztr.FOLLOW1}]
        assert ztr.unpack(ztr.pack(self.PAYLOAD, chain))[0] == self.PAYLOAD

    def test_full_chain_roundtrip(self):
        """The 5-layer chain SnapGene uses for SMP4."""
        chain = [
            {"format": ztr.ZLIB},
            {"format": ztr.RLE, "guard": 0x96},
            {"format": ztr.FOLLOW1},
            {"format": ztr.SHRINK_16TO8},
            {"format": ztr.DELTA2, "level": 3},
        ]
        packed = ztr.pack(self.PAYLOAD, chain)
        payload, recovered = ztr.unpack(packed)
        assert payload == self.PAYLOAD
        assert recovered == chain

    def test_raw_payload_needs_no_filters(self):
        assert ztr.unpack(self.PAYLOAD) == (self.PAYLOAD, [])
        assert ztr.pack(self.PAYLOAD, []) == self.PAYLOAD

    def test_unsupported_filter_raises(self):
        """STHUFF and friends are rejected rather than misread as data."""
        with pytest.raises(ztr.ZtrFormatError, match="STHUFF"):
            ztr.unpack(bytes([77, 0]) + b"garbage")

    def test_unsupported_filter_cannot_be_packed(self):
        with pytest.raises(ztr.ZtrFormatError):
            ztr.pack(b"\x00data", [{"format": 79}])


# =============================================================================
# Decoded values
# =============================================================================


class TestTraceDecoding:
    def test_positions_are_monotonic(self, legacy_v13_dna):
        """A chromatogram's peak positions strictly increase along the trace.

        Decoding only the zlib layer produced unordered nonsense here.
        """
        sgff = SgffReader.from_file(legacy_v13_dna)
        for trace in sgff.traces:
            assert trace.positions
            assert all(
                a < b for a, b in zip(trace.positions, trace.positions[1:])
            )

    def test_one_position_per_base(self, legacy_v13_dna):
        sgff = SgffReader.from_file(legacy_v13_dna)
        for trace in sgff.traces:
            assert len(trace.positions) == len(trace.bases)

    def test_four_confidence_values_per_base(self, legacy_v13_dna):
        """CNF4 stores the called-base score plus one per A/C/G/T."""
        sgff = SgffReader.from_file(legacy_v13_dna)
        for trace in sgff.traces:
            assert len(trace.confidence) == 4 * len(trace.bases)

    def test_samples_within_16_bit_range_and_channels_equal(self, legacy_v13_dna):
        sgff = SgffReader.from_file(legacy_v13_dna)
        for trace in sgff.traces:
            samples = trace.samples
            assert samples is not None
            lengths = {len(samples.a), len(samples.c), len(samples.g), len(samples.t)}
            assert len(lengths) == 1
            assert max(samples.a) <= 0xFFFF

    def test_positions_index_into_the_sample_trace(self, legacy_v13_dna):
        """Every base's peak must fall inside the sampled trace."""
        sgff = SgffReader.from_file(legacy_v13_dna)
        for trace in sgff.traces:
            assert trace.positions[-1] < len(trace.samples)

    def test_encoding_records_chain_and_order(self, legacy_v13_dna):
        sgff = SgffReader.from_file(legacy_v13_dna)
        encoding = sgff.traces[0].encoding
        assert encoding["version"] == [1, 2]
        order = [chunk["type"] for chunk in encoding["chunks"]]
        # SnapGene writes samples first; the writer must not re-sort these
        assert order[0] == "SMP4"
        assert set(order) >= {"SMP4", "BASE", "BPOS", "CNF4"}


# =============================================================================
# Byte-exact re-encoding
# =============================================================================


class TestZtrByteIdentity:
    @pytest.mark.parametrize("filename", ["legacy_v13.dna", "test3.dna"])
    def test_every_chunk_repacks_identically(self, data_dir, filename):
        """Replaying a recorded chain must reproduce the source bytes exactly."""
        blobs = collect_ztr_blobs(data_dir / filename)
        assert blobs
        checked = 0
        for blob in blobs:
            for chunk_type, _meta, encoded in iter_chunks(blob):
                payload, chain = ztr.unpack(encoded)
                assert ztr.pack(payload, chain) == encoded, (
                    f"{filename}: {chunk_type} chunk did not repack identically"
                )
                checked += 1
        assert checked > 0

    @pytest.mark.parametrize("filename", ["legacy_v13.dna", "test3.dna"])
    def test_trace_blocks_survive_roundtrip_byte_for_byte(self, data_dir, filename):
        """Issue #23: block 16 came back ~1.9x larger and full of noise."""
        path = data_dir / filename
        original = path.read_bytes()
        sgff = SgffReader.from_file(path)
        rewritten = SgffWriter.to_bytes(sgff, preserve=True)

        before = block_16_bytes(original)
        after = block_16_bytes(rewritten)
        assert before, f"{filename} has no block 16 to compare"
        assert len(before) == len(after)
        for index, (a, b) in enumerate(zip(before, after)):
            assert a == b, f"{filename}: trace container {index} changed on write"

    def test_all_filters_are_exercised_by_the_fixtures(self, data_dir):
        """Guard against a fixture change silently narrowing this coverage."""
        seen = set()
        for blob in collect_ztr_blobs(data_dir / "legacy_v13.dna"):
            for _type, _meta, encoded in iter_chunks(blob):
                seen.update(layer["format"] for layer in ztr.unpack(encoded)[1])
        assert seen == {
            ztr.ZLIB,
            ztr.RLE,
            ztr.DELTA1,
            ztr.DELTA2,
            ztr.DELTA4,
            ztr.SHRINK_16TO8,
            ztr.SHRINK_32TO8,
            ztr.FOLLOW1,
        }


# =============================================================================
# Chunk-level edge cases
# =============================================================================


def build_ztr(chunks):
    """Assemble a ZTR blob from (type, meta, encoded_data) triples."""
    blob = bytearray(ZTR_MAGIC + bytes([1, 2]))
    for chunk_type, meta, data in chunks:
        blob += chunk_type.encode("ascii").ljust(4)[:4]
        blob += struct.pack(">I", len(meta)) + meta
        blob += struct.pack(">I", len(data)) + data
    return bytes(blob)


class TestChunkEdgeCases:
    def _rewrite(self, blob):
        parsed = parse_ztr(blob)
        return SgffWriter(None)._serialize_ztr(parsed)

    def test_empty_text_chunk_is_preserved(self):
        """A TEXT chunk with no pairs is still a chunk; dropping it changes bytes."""
        blob = build_ztr([("TEXT", b"", b"\x00\x00"), ("CLIP", b"", b"\x00" + b"\x00" * 8)])
        assert self._rewrite(blob) == blob

    def test_text_pair_list_terminator_survives(self):
        payload = b"\x00NAME\x00abc\x00LANE\x0033\x00\x00"
        blob = build_ztr([("TEXT", b"", payload)])
        assert self._rewrite(blob) == blob
        assert parse_ztr(blob)["text"] == {"NAME": "abc", "LANE": "33"}

    def test_unsupported_filter_passes_through_untouched(self):
        """An unknown filter keeps its bytes rather than being dropped or misread."""
        exotic = bytes([77, 0]) + b"opaque huffman payload"
        blob = build_ztr([("BASE", b"", b"\x00ACGT"), ("CNF4", b"", exotic)])
        parsed = parse_ztr(blob)
        assert parsed["bases"] == "ACGT"
        assert "confidence" not in parsed
        assert self._rewrite(blob) == blob

    def test_unmodelled_chunk_type_passes_through(self):
        blob = build_ztr([("BASE", b"", b"\x00ACGT"), ("XYZZ", b"", b"\x00payload")])
        assert self._rewrite(blob) == blob

    def test_samp_chunks_carry_channel_metadata(self):
        chunks = [
            ("SAMP", ch.encode() + b"\x00\x00\x00", b"\x00\x00" + struct.pack(">2H", 5, 9))
            for ch in "ACGT"
        ]
        blob = build_ztr(chunks)
        parsed = parse_ztr(blob)
        assert parsed["samples"] == {c: [5, 9] for c in "ACGT"}
        assert self._rewrite(blob) == blob

    def test_multiple_comm_chunks_keep_their_order(self):
        blob = build_ztr([("COMM", b"", b"\x00first"), ("COMM", b"", b"\x00second")])
        assert parse_ztr(blob)["comments"] == ["first", "second"]
        assert self._rewrite(blob) == blob

    def test_trace_built_from_scratch_is_written_unfiltered(self):
        """No recorded encoding means raw chunks in SnapGene's order."""
        data = {
            "bases": "ACGT",
            "positions": [1, 2, 3, 4],
            "confidence": [10] * 16,
            "samples": {ch: [1, 2, 3] for ch in "ACGT"},
        }
        blob = SgffWriter(None)._serialize_ztr(data)
        order = [chunk_type for chunk_type, _m, _d in iter_chunks(blob)]
        assert order == ["SMP4", "BASE", "BPOS", "CNF4"]
        reparsed = parse_ztr(blob)
        assert reparsed["bases"] == "ACGT"
        assert reparsed["positions"] == [1, 2, 3, 4]
        assert reparsed["confidence"] == [10] * 16
        assert reparsed["samples"] == {ch: [1, 2, 3] for ch in "ACGT"}
