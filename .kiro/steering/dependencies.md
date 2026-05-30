# Dependencies & Tooling

## Runtime dependencies
- `click>=8.0` — CLI framework (replaced argparse)
- `numpy>=1.12.1` — array math
- `scipy>=0.19.0` — NNLS solver (`scipy.optimize.nnls`)
- `jellyfish` (external binary, 2.x) — k-mer counting, NOT a pip package

## Dev dependencies
- `pytest>=7.0`
- `pytest-cov>=4.0`

## Installation
```bash
# recommended (dev)
uv pip install -e ".[dev]"

# standard
pip install -e .
```
uv is a developer convenience tool only — NOT a runtime or build dependency.
The package installs with plain pip/conda.

## Running tests
```bash
pytest                                          # all tests
pytest --cov=focus_app --cov-report=term-missing  # with coverage
```
Current coverage: 99%. 40 tests total.

## Test files
- `test_focus.py` — unit tests for every function in focus.py
- `test_cli.py` — CLI tests using `click.testing.CliRunner`

## Bioconda
- Recipe: `recipe/meta.yaml`
- Package: `noarch: python`
- Entry point: `focus = focus_app.cli:main`
- `jellyfish` listed as a conda run dependency (available in bioconda channel)
- Before submitting: fill in `sha256` and real GitHub URL in `meta.yaml`

## Version
Current: `2.0.0` (major bump from FOCUS 1.6 original).
Defined in `focus_app/__init__.py` and `pyproject.toml`.
