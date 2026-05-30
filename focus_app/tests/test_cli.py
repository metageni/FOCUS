#!/usr/bin/env python3
"""Tests for focus_app.cli and the profile() pipeline function.

Uses Click's CliRunner for all main() invocations.
"""

import pytest

from pathlib import Path
from unittest.mock import patch, MagicMock

import numpy as np

from click.testing import CliRunner

from focus_app.cli import (
    _get_jellyfish_version,
    _collect_query_files,
    main,
)
from focus_app.focus import profile

DATA = Path(__file__).parent / "data"
DB_SMALL = DATA / "k6_small_sample.txt"
FASTA = DATA / "mock_sample.fasta"
FASTQ = DATA / "mock_sample.fastq"

runner = CliRunner()


# ---------------------------------------------------------------------------
# _get_jellyfish_version
# ---------------------------------------------------------------------------

def test_get_jellyfish_version_none_when_no_path():
    """Returns None when jellyfish binary path is falsy."""
    assert _get_jellyfish_version(None) is None
    assert _get_jellyfish_version("") is None


def test_get_jellyfish_version_returns_major():
    """Returns the major version digit string."""
    assert _get_jellyfish_version("/opt/homebrew/bin/jellyfish") == "2"


# ---------------------------------------------------------------------------
# _collect_query_files
# ---------------------------------------------------------------------------

def test_collect_query_files_single_file():
    """A single valid file path is returned as a one-element list."""
    result = _collect_query_files([str(FASTA)])
    assert len(result) == 1 and result[0] == FASTA


def test_collect_query_files_directory(tmp_path):
    """Files inside a directory are discovered and filtered."""
    (tmp_path / "a.fasta").write_text(">r\nAAAA\n")
    (tmp_path / "b.txt").write_text("ignore")
    result = _collect_query_files([str(tmp_path)])
    assert len(result) == 1 and result[0].name == "a.fasta"


def test_collect_query_files_mixed(tmp_path):
    """Mix of file and directory args are merged and filtered."""
    (tmp_path / "c.fastq").write_text("@r\nAAAA\n+\nIIII\n")
    result = _collect_query_files([str(FASTA), str(tmp_path)])
    names = [p.name for p in result]
    assert "mock_sample.fasta" in names and "c.fastq" in names


def test_collect_query_files_no_valid_files(tmp_path):
    """Returns empty list when no valid sequence files exist."""
    (tmp_path / "readme.txt").write_text("nothing")
    assert _collect_query_files([str(tmp_path)]) == []


# ---------------------------------------------------------------------------
# main() — Click CLI validation exit paths
# ---------------------------------------------------------------------------

def test_main_exits_no_query_files(tmp_path):
    """main() exits non-zero when query dir has no sequence files."""
    (tmp_path / "junk.txt").write_text("x")
    result = runner.invoke(main, ["-q", str(tmp_path), "-o", str(tmp_path / "out")])
    assert result.exit_code != 0


def test_main_exits_db_not_found(tmp_path):
    """main() exits non-zero when the database path does not exist."""
    result = runner.invoke(main, ["-q", str(FASTA), "-o", str(tmp_path / "out"),
                                  "-b", str(tmp_path)])
    assert result.exit_code != 0


def test_main_exits_jellyfish_missing(tmp_path):
    """main() exits non-zero when jellyfish is not on PATH."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    (db_dir / "k6").write_text("")

    with patch("focus_app.cli.which", return_value=None):
        result = runner.invoke(main, ["-q", str(FASTA), "-o", str(tmp_path / "out"),
                                      "-b", str(tmp_path)])
    assert result.exit_code != 0


def test_main_exits_jellyfish_wrong_version(tmp_path):
    """main() exits non-zero when jellyfish major version is not 2."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    (db_dir / "k6").write_text("")

    with patch("focus_app.cli.which", return_value="/usr/bin/jellyfish"), \
         patch("focus_app.cli._get_jellyfish_version", return_value="1"):
        result = runner.invoke(main, ["-q", str(FASTA), "-o", str(tmp_path / "out"),
                                      "-b", str(tmp_path)])
    assert result.exit_code != 0


def test_main_exits_invalid_kmer_size(tmp_path):
    """Click rejects k-mer size other than 6 or 7 (exit code 2)."""
    result = runner.invoke(main, ["-q", str(FASTA), "-o", str(tmp_path / "out"), "-k", "5"])
    assert result.exit_code == 2


def test_main_runs_and_writes_outputs(tmp_path):
    """main() completes successfully and writes all expected CSV files."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    import shutil
    shutil.copy(DB_SMALL, db_dir / "k6")

    fake_results = {
        "Bacteria\tFirmicutes\tBacilli\tLactobacillales\tStreptococcaceae"
        "\tStreptococcus\tStreptococcus_suis\tStreptococcus_suis_uid1": [0.8],
    }

    with patch("focus_app.cli.which", return_value="/usr/bin/jellyfish"), \
         patch("focus_app.cli._get_jellyfish_version", return_value="2"), \
         patch("focus_app.cli.profile", return_value=(fake_results, list(fake_results))):
        result = runner.invoke(main, ["-q", str(FASTA), "-o", str(tmp_path / "out"),
                                      "-b", str(tmp_path)])

    assert result.exit_code == 0
    out = tmp_path / "out"
    assert (out / "output_All_levels.csv").exists()
    assert (out / "output_Kingdom_tabular.csv").exists()
    assert (out / "output_Strain_tabular.csv").exists()


def test_main_list_output_returns_rows(tmp_path):
    """main() returns a list of lists when --list_output is set."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    import shutil
    shutil.copy(DB_SMALL, db_dir / "k6")

    fake_results = {
        "Bacteria\tFirmicutes\tBacilli\tLactobacillales\tStreptococcaceae"
        "\tStreptococcus\tStreptococcus_suis\tStreptococcus_suis_uid1": [0.8],
    }

    with patch("focus_app.cli.which", return_value="/usr/bin/jellyfish"), \
         patch("focus_app.cli._get_jellyfish_version", return_value="2"), \
         patch("focus_app.cli.profile", return_value=(fake_results, list(fake_results))):
        result = runner.invoke(main, ["-q", str(FASTA), "-o", str(tmp_path / "out"),
                                      "-b", str(tmp_path), "--list_output"])

    assert result.exit_code == 0


def test_main_log_file(tmp_path):
    """main() passes the log filename to logging.basicConfig when -l is provided."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    import shutil
    shutil.copy(DB_SMALL, db_dir / "k6")
    log_file = tmp_path / "run.log"

    fake_results = {
        "Bacteria\tFirmicutes\tBacilli\tLactobacillales\tStreptococcaceae"
        "\tStreptococcus\tStreptococcus_suis\tStreptococcus_suis_uid1": [0.8],
    }
    captured = {}

    with patch("focus_app.cli.which", return_value="/usr/bin/jellyfish"), \
         patch("focus_app.cli._get_jellyfish_version", return_value="2"), \
         patch("focus_app.cli.profile", return_value=(fake_results, list(fake_results))), \
         patch("focus_app.cli.logging.basicConfig", side_effect=lambda **kw: captured.update(kw)):
        runner.invoke(main, ["-q", str(FASTA), "-o", str(tmp_path / "out"),
                             "-b", str(tmp_path), "-l", str(log_file)])

    assert captured.get("filename") == str(log_file)


# ---------------------------------------------------------------------------
# profile()
# ---------------------------------------------------------------------------

def test_profile_returns_correct_structure():
    """profile() returns a results dict and organism list with correct shapes."""
    with patch("focus_app.focus.count_kmers", return_value=[990, 1439, 1320]):
        results, organisms = profile([FASTA], DB_SMALL, kmer_size="6", threads="1")
    assert isinstance(results, dict)
    assert len(organisms) == 2
    assert all(len(v) == 1 for v in results.values())


def test_profile_abundances_sum_to_one():
    """Abundances across all organisms for a single query must sum to ~1."""
    with patch("focus_app.focus.count_kmers", return_value=[990, 1439, 1320]):
        results, organisms = profile([FASTA], DB_SMALL, kmer_size="6", threads="1")
    np.testing.assert_almost_equal(
        sum(results[org][0] for org in organisms), 1.0, decimal=6)


def test_profile_multiple_queries():
    """profile() handles multiple query files, producing one column per file."""
    with patch("focus_app.focus.count_kmers", return_value=[548, 753, 661]):
        results, organisms = profile([FASTA, FASTQ], DB_SMALL, kmer_size="6", threads="1")
    assert all(len(v) == 2 for v in results.values())


# ---------------------------------------------------------------------------
# count_kmers error branches (mocked)
# ---------------------------------------------------------------------------

def test_count_kmers_dump_missing_raises():
    """Raises Exception when jellyfish dump produces no output file."""
    from focus_app.focus import count_kmers

    def fake_exists(self):
        return "kmer_counting" in str(self)

    with patch("focus_app.focus.random.random", return_value=0.0), \
         patch("focus_app.focus.os.system"), \
         patch("focus_app.focus.Path.exists", fake_exists), \
         patch("focus_app.focus.os.remove"):
        with pytest.raises(Exception, match="dump step failed"):
            count_kmers(FASTA, "6", "1", ["AAAAAA"])


def test_count_kmers_empty_dump_raises():
    """Raises Exception when jellyfish dump file exists but is empty."""
    from focus_app.focus import count_kmers

    def fake_exists(self):
        return True

    def fake_stat(self):
        m = MagicMock()
        m.st_size = 0
        return m

    with patch("focus_app.focus.random.random", return_value=0.0), \
         patch("focus_app.focus.os.system"), \
         patch("focus_app.focus.Path.exists", fake_exists), \
         patch("focus_app.focus.Path.stat", fake_stat), \
         patch("focus_app.focus.os.remove"):
        with pytest.raises(Exception, match="no k-mer counts"):
            count_kmers(FASTA, "6", "1", ["AAAAAA"])
