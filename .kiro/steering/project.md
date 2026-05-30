# FOCUS — Project Overview

## What this project is

FOCUS (version 2.0.0) is a rewrite of the original FOCUS metagenomics profiler (v1.6).
It is an alignment-free tool that estimates the relative abundance of organisms in a
metagenome using k-mer counting (Jellyfish) and non-negative least squares (NNLS).

## Package structure

```
focus/                          ← repo root
├── focus_app/                  ← Python package (named focus_app to match original)
│   ├── __init__.py             ← version = "2.0.0"
│   ├── focus.py                ← ALL core logic (pure functions, no CLI)
│   ├── cli.py                  ← CLI only, built with Click
│   └── tests/
│       ├── data/               ← test fixtures (FASTA, FASTQ, small DB files)
│       ├── test_focus.py       ← unit tests for focus.py
│       └── test_cli.py         ← unit tests for cli.py (uses Click CliRunner)
├── recipe/
│   └── meta.yaml               ← Bioconda conda recipe
├── pyproject.toml              ← build config, deps, entry points
├── Dockerfile                  ← miniconda3-based, local use
├── Docker/
│   ├── Dockerfile              ← ubuntu:24.04 + miniforge, HPC/Singularity
│   └── README.md
├── README.md
├── CHANGELOG
├── COPYING
├── logo/
│   └── focus_small_logo.png
└── MANIFEST.in
```

## Entry point

```
focus = focus_app.cli:main
```

The `focus` command is installed via `pyproject.toml`. Run with:
```bash
focus -q samples/ -o results/
```

## Database

The k-mer database lives at `focus_app/db/` (extracted from `db.zip`).
It is NOT committed to git (too large). Users must extract it:
```bash
cd focus_app && unzip db.zip
```

## Key design decisions

- `focus.py` is pure logic — no argparse, no sys.exit, no logging.basicConfig
- `cli.py` is CLI only — all Click decorators, validation, and logging setup live here
- `profile()` in `focus.py` is the clean programmatic entry point
- Tests use `click.testing.CliRunner` for CLI tests, not sys.argv patching
