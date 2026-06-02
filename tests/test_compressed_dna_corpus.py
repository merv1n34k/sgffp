"""Corpus tests for the compressed-DNA parser/writer.

Two binary corpora live under ``tests/data``:

* ``degenerate_examples/`` — small fixtures where the live block-0 sequence
  equals the embedded block-11 snapshot at history node 0. The decoded
  snapshot must match the live sequence verbatim.
* ``format_31_examples/<case>/`` — for each case the assembly product
  ``out.dna`` contains a history snapshot of the original ``vector.dna``.
  The decoded snapshot must equal ``vector.dna``'s live sequence.
"""

from pathlib import Path

import pytest

from sgffp.reader import SgffReader
from sgffp.writer import SgffWriter

DATA = Path(__file__).parent / "data"
DEGEN = DATA / "degenerate_examples"
FMT31 = DATA / "format_31_examples"


def _vector_sequence(case_dir: Path) -> str:
    return SgffReader.from_file(case_dir / "vector.dna").blocks[0][0]["sequence"]


def _snapshot_sequence_for_vector(out_dna: Path, expected_length: int) -> str:
    sgff = SgffReader.from_file(out_dna)
    for node in sgff.blocks.get(11, []):
        if node.get("sequence_type") == 1 and node.get("length") == expected_length:
            return node["sequence"]
    raise AssertionError(f"no matching block-11 snapshot in {out_dna}")


@pytest.mark.parametrize("dna_path", sorted(DEGEN.glob("*.dna")))
def test_degenerate_examples_roundtrip(dna_path):
    """Block-0 live sequence equals the first block-11 snapshot."""
    sgff = SgffReader.from_file(dna_path)
    live = sgff.blocks[0][0]["sequence"]
    snapshot = sgff.history.get_sequence_at(0)
    assert snapshot == live, f"{dna_path.name}: snapshot != live"


@pytest.mark.parametrize(
    "case",
    [
        "longB",
        "longN",
        "longW",
        "lowercase",
        "mixed",
        "singleN-beginning",
        "singleN-end",
    ],
)
def test_format_31_decodes_vector(case):
    """Assembly output's block-11 snapshot equals the source vector sequence."""
    case_dir = FMT31 / case
    vector_seq = _vector_sequence(case_dir)
    snapshot_seq = _snapshot_sequence_for_vector(
        case_dir / "out.dna", len(vector_seq)
    )
    assert snapshot_seq == vector_seq


@pytest.mark.parametrize("dna_path", sorted(DEGEN.glob("*.dna")))
def test_degenerate_examples_full_roundtrip(dna_path):
    """Read → write → read preserves the live sequence and every block-11 snapshot."""
    original = SgffReader.from_file(dna_path)
    written = SgffWriter.to_bytes(original)
    restored = SgffReader.from_bytes(written)

    assert restored.blocks[0][0]["sequence"] == original.blocks[0][0]["sequence"]

    orig_snaps = [
        n["sequence"]
        for n in original.blocks.get(11, [])
        if n.get("sequence_type") == 1
    ]
    rest_snaps = [
        n["sequence"]
        for n in restored.blocks.get(11, [])
        if n.get("sequence_type") == 1
    ]
    assert orig_snaps == rest_snaps


@pytest.mark.parametrize(
    "case",
    ["longB", "longN", "longW", "lowercase", "mixed", "singleN-beginning", "singleN-end"],
)
def test_format_31_full_roundtrip(case):
    """Format-31 cases round-trip cleanly through our writer."""
    src = FMT31 / case / "out.dna"
    original = SgffReader.from_file(src)
    written = SgffWriter.to_bytes(original)
    restored = SgffReader.from_bytes(written)

    orig_snaps = [
        n["sequence"]
        for n in original.blocks.get(11, [])
        if n.get("sequence_type") == 1
    ]
    rest_snaps = [
        n["sequence"]
        for n in restored.blocks.get(11, [])
        if n.get("sequence_type") == 1
    ]
    assert orig_snaps == rest_snaps


def test_format_31_double_roundtrip_stable():
    """Output stabilises after one read/write cycle on the hardest case."""
    src = FMT31 / "mixed" / "out.dna"
    written1 = SgffWriter.to_bytes(SgffReader.from_file(src))
    written2 = SgffWriter.to_bytes(SgffReader.from_bytes(written1))
    assert written1 == written2
