"""
Tests for CLI commands
"""

import json
import subprocess
import sys
import tempfile
from io import BytesIO
from pathlib import Path
from unittest.mock import patch

import pytest

from sgffp.cli import cmd_parse, cmd_info, cmd_tree, cmd_check, cmd_filter, main


class MockArgs:
    """Mock argparse namespace"""

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


# =============================================================================
# Parse Command Tests
# =============================================================================


class TestParseCommand:
    def test_parse_stdout(self, test_dna, capsys):
        """JSON to stdout"""
        args = MockArgs(input=str(test_dna), output=None)
        cmd_parse(args)

        captured = capsys.readouterr()
        output = json.loads(captured.out)

        assert "cookie" in output
        assert "blocks" in output

    def test_parse_file_output(self, test_dna):
        """JSON to file with -o"""
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            output_path = f.name

        try:
            args = MockArgs(input=str(test_dna), output=output_path)
            cmd_parse(args)

            assert Path(output_path).exists()

            with open(output_path) as f:
                output = json.load(f)

            assert "cookie" in output
            assert "blocks" in output
        finally:
            Path(output_path).unlink(missing_ok=True)

    def test_parse_has_cookie_fields(self, test_dna, capsys):
        """Cookie contains expected fields"""
        args = MockArgs(input=str(test_dna), output=None)
        cmd_parse(args)

        captured = capsys.readouterr()
        output = json.loads(captured.out)

        cookie = output["cookie"]
        assert "type_of_sequence" in cookie
        assert "export_version" in cookie
        assert "import_version" in cookie


# =============================================================================
# Info Command Tests
# =============================================================================


class TestInfoCommand:
    def test_info_output(self, test_dna, capsys):
        """Metadata display"""
        args = MockArgs(input=str(test_dna))
        cmd_info(args)

        captured = capsys.readouterr()
        assert "File:" in captured.out
        assert "Type:" in captured.out
        assert "export v" in captured.out
        assert "import v" in captured.out
        assert "Sequence:" in captured.out

    def test_info_block_counts(self, test_dna, capsys):
        """Block counts shown"""
        args = MockArgs(input=str(test_dna))
        cmd_info(args)

        captured = capsys.readouterr()
        assert "Blocks:" in captured.out

    def test_info_verbose(self, test3_dna, capsys):
        """Verbose info shows detailed sections"""
        args = MockArgs(input=str(test3_dna), verbose=True)
        cmd_info(args)

        captured = capsys.readouterr()
        assert "--- Notes ---" in captured.out
        assert "--- Primers" in captured.out
        assert "--- History ---" in captured.out

    def test_info_verbose_minimal(self, test_dna, capsys):
        """Verbose on minimal file doesn't crash"""
        args = MockArgs(input=str(test_dna), verbose=True)
        cmd_info(args)

        captured = capsys.readouterr()
        assert "File:" in captured.out


# =============================================================================
# Check Command Tests
# =============================================================================


class TestCheckCommand:
    def test_check_lists_blocks(self, test_dna, capsys):
        """Check with -l flag lists blocks"""
        args = MockArgs(input=str(test_dna), list=True, dump=False)
        cmd_check(args)

        captured = capsys.readouterr()
        # Should list block types with counts
        # Format: "  0:  1"
        assert ":" in captured.out

    def test_check_without_flags(self, test_dna, capsys):
        """Check without flags runs without error"""
        args = MockArgs(input=str(test_dna), list=False, dump=False)
        cmd_check(args)
        # Should not raise an exception


# =============================================================================
# Tree Command Tests
# =============================================================================


class TestTreeCommand:
    def test_tree_output(self, test3_dna, capsys):
        """Tree displays timeline with operations"""
        args = MockArgs(input=str(test3_dna), verbose=False)
        cmd_tree(args)

        captured = capsys.readouterr()
        assert "test.rna" in captured.out
        assert "makeDna" in captured.out
        assert "amplifyFragment" in captured.out
        assert "test_ampf.dna" in captured.out

    def test_tree_verbose(self, test3_dna, capsys):
        """Verbose tree shows operation details"""
        args = MockArgs(input=str(test3_dna), verbose=True)
        cmd_tree(args)

        captured = capsys.readouterr()
        assert "oligo:" in captured.out

    def test_tree_no_history(self, test_dna, capsys):
        """Tree on file without history exits with error"""
        # Create a minimal file with only block 0
        with tempfile.NamedTemporaryFile(suffix=".dna", delete=False) as f:
            output_path = f.name

        try:
            from sgffp.cli import cmd_filter

            filter_args = MockArgs(input=str(test_dna), keep="0", output=output_path)
            cmd_filter(filter_args)
            capsys.readouterr()  # clear filter output

            args = MockArgs(input=output_path, verbose=False)
            with pytest.raises(SystemExit) as exc_info:
                cmd_tree(args)
            assert exc_info.value.code == 1
        finally:
            Path(output_path).unlink(missing_ok=True)


# =============================================================================
# Filter Command Tests
# =============================================================================


class TestFilterCommand:
    def test_filter_single_type(self, test_dna, capsys):
        """Filter to single block type"""
        with tempfile.NamedTemporaryFile(suffix=".dna", delete=False) as f:
            output_path = f.name

        try:
            args = MockArgs(input=str(test_dna), keep="0", output=output_path)
            cmd_filter(args)

            captured = capsys.readouterr()
            assert "Filtered file written to" in captured.out

            # Verify output file exists and is valid
            assert Path(output_path).exists()
            assert Path(output_path).stat().st_size > 0
        finally:
            Path(output_path).unlink(missing_ok=True)

    def test_filter_multiple_types(self, test_dna, capsys):
        """Filter to multiple block types"""
        with tempfile.NamedTemporaryFile(suffix=".dna", delete=False) as f:
            output_path = f.name

        try:
            args = MockArgs(input=str(test_dna), keep="0,6,10", output=output_path)
            cmd_filter(args)

            assert Path(output_path).exists()
        finally:
            Path(output_path).unlink(missing_ok=True)


# =============================================================================
# Stdin Tests
# =============================================================================


class TestStdin:
    def test_info_stdin(self, test_dna, capsys):
        """Info reads from stdin when input is -"""
        data = test_dna.read_bytes()
        with patch("sys.stdin", **{"buffer": BytesIO(data)}):
            args = MockArgs(input="-")
            cmd_info(args)

        captured = capsys.readouterr()
        assert "File: <stdin>" in captured.out
        assert "Type:" in captured.out

    def test_parse_stdin(self, test_dna, capsys):
        """Parse reads from stdin when input is -"""
        data = test_dna.read_bytes()
        with patch("sys.stdin", **{"buffer": BytesIO(data)}):
            args = MockArgs(input="-", output=None)
            cmd_parse(args)

        captured = capsys.readouterr()
        output = json.loads(captured.out)
        assert "cookie" in output
        assert "blocks" in output

    def test_check_stdin(self, test_dna, capsys):
        """Check reads from stdin when input is -"""
        data = test_dna.read_bytes()
        with patch("sys.stdin", **{"buffer": BytesIO(data)}):
            args = MockArgs(input="-", list=True, dump=False)
            cmd_check(args)

        captured = capsys.readouterr()
        assert ":" in captured.out


# =============================================================================
# Main Function Tests
# =============================================================================


class TestMain:
    def test_main_no_args(self, capsys):
        """No arguments shows help and exits"""
        with patch("sys.argv", ["sff"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 1

    def test_main_info_command(self, test_dna, capsys):
        """Main with info command"""
        with patch("sys.argv", ["sff", "info", str(test_dna)]):
            main()

        captured = capsys.readouterr()
        assert "File:" in captured.out


# =============================================================================
# Subprocess Integration Tests
# =============================================================================


class TestSubprocess:
    def test_cli_module_execution(self, test_dna):
        """Run CLI as module"""
        result = subprocess.run(
            [sys.executable, "-m", "sgffp.cli", "info", str(test_dna)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "File:" in result.stdout

    def test_cli_parse_json(self, test_dna):
        """Parse command outputs valid JSON"""
        result = subprocess.run(
            [sys.executable, "-m", "sgffp.cli", "parse", str(test_dna)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0

        # Output should be valid JSON
        output = json.loads(result.stdout)
        assert "cookie" in output
        assert "blocks" in output

    def test_cli_check_list(self, test_dna):
        """Check command with -l flag"""
        result = subprocess.run(
            [sys.executable, "-m", "sgffp.cli", "check", "-l", str(test_dna)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
