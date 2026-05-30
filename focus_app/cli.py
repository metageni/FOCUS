#!/usr/bin/env python3
"""CLI entry point for focus_app — built with Click, separate from core logic."""

import os
import logging

from pathlib import Path
from shutil import which

import click

from focus_app import __version__
from focus_app.focus import (
    TAXONOMY_LEVELS,
    is_wanted_file,
    profile,
    aggregate_level,
    write_results,
    refine_results,
)

LOGGER_FORMAT = "[%(asctime)s - %(levelname)s] %(message)s"


def _get_jellyfish_version(jf_path):
    """Return the major version string of the installed Jellyfish, or None.

    Args:
        jf_path (str or None): Path to the jellyfish binary.

    Returns:
        str or None: Major version digit, e.g. "2", or None if not found.
    """
    if not jf_path:
        return None
    return os.popen("jellyfish count --version").read().split(".")[0]


def _collect_query_files(query_args):
    """Expand query arguments (files and/or directories) into a flat file list.

    Args:
        query_args (tuple[str]): Paths provided via -q flags.

    Returns:
        list[Path]: Filtered, sorted list of valid sequence files.
    """
    raw = []
    for q in query_args:
        p = Path(q)
        if p.is_dir():
            raw += [p / f for f in os.listdir(p)]
        elif p.is_file():
            raw.append(p)
    return is_wanted_file(raw)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, "-v", "--version", prog_name="focus")
@click.option("-q", "--query", multiple=True, required=True,
              help="FAST(A/Q) file or directory (repeatable).")
@click.option("-o", "--output_directory", required=True,
              help="Directory for output CSV files.")
@click.option("-k", "--kmer_size", default="6", show_default=True,
              type=click.Choice(["6", "7"]), help="K-mer size.")
@click.option("-b", "--alternate_directory", default="",
              help="Alternate directory containing the db/ folder.")
@click.option("-p", "--output_prefix", default="output", show_default=True,
              help="Prefix for output files.")
@click.option("-t", "--threads", default="4", show_default=True,
              help="Jellyfish thread count.")
@click.option("--list_output", is_flag=True,
              help="Return results as a list of lists.")
@click.option("-l", "--log", default=None,
              help="Log file path (default: STDOUT).")
def main(query, output_directory, kmer_size, alternate_directory,
         output_prefix, threads, list_output, log):
    """FOCUS: An Agile Profiler for Metagenomic Data.

    Example:

        focus -q samples/ -o results/
    """
    log_kwargs = dict(format=LOGGER_FORMAT, level=logging.INFO)
    if log:
        log_kwargs["filename"] = log
    logging.basicConfig(**log_kwargs)
    logger = logging.getLogger(__name__)

    logger.info("focus version %s", __version__)

    output_dir = Path(output_directory)
    output_dir.mkdir(parents=True, exist_ok=True)

    work_dir = Path(alternate_directory) if alternate_directory else Path(__file__).parent
    db_path = work_dir / "db" / "k{}".format(kmer_size)

    query_files = _collect_query_files(query)
    if not query_files:
        logger.critical("No FASTA/FASTQ files found in: %s", list(query))
        raise SystemExit(1)

    if not db_path.exists():
        logger.critical("Database not found: %s — extract db.zip first.", db_path)
        raise SystemExit(1)

    jf_path = which("jellyfish")
    jf_version = _get_jellyfish_version(jf_path)
    if not jf_path:
        logger.critical("Jellyfish is not installed.")
        raise SystemExit(1)
    if jf_version != "2":
        logger.critical("Jellyfish 2.x required; found version %s.", jf_version)
        raise SystemExit(1)

    results, _ = profile(query_files, db_path, kmer_size, threads)

    write_results(results, output_dir / "{}_All_levels.csv".format(output_prefix),
                  query_files, TAXONOMY_LEVELS)

    for pos, level in enumerate(TAXONOMY_LEVELS):
        level_results = aggregate_level(results, pos)
        write_results(level_results,
                      output_dir / "{}_{}_tabular.csv".format(output_prefix, level),
                      query_files, [level])

    logger.info("Results written to %s", output_dir)

    if list_output:
        return refine_results(results, query_files, TAXONOMY_LEVELS)


if __name__ == "__main__":
    main()
