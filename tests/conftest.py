"""
Shared fixtures for SGFFP tests
"""

import pytest
from pathlib import Path


@pytest.fixture
def data_dir():
    """Path to test data directory"""
    return Path(__file__).parent / "data"


@pytest.fixture
def test_dna(data_dir):
    """Path to test.dna file"""
    return data_dir / "test.dna"


@pytest.fixture
def test2_dna(data_dir):
    """Path to test2.dna file"""
    return data_dir / "test2.dna"


@pytest.fixture
def test3_dna(data_dir):
    """Path to test3.dna file"""
    return data_dir / "test3.dna"


@pytest.fixture
def test_rna(data_dir):
    """Path to test.rna file"""
    return data_dir / "test.rna"


@pytest.fixture
def test_prot(data_dir):
    """Path to test.prot file"""
    return data_dir / "test.prot"


@pytest.fixture
def all_test_files(test_dna, test2_dna, test3_dna, test_rna, test_prot):
    """All test files for parametrized tests"""
    return [test_dna, test2_dna, test3_dna, test_rna, test_prot]
