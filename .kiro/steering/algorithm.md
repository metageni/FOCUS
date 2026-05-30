# Algorithm & Science

## What FOCUS does
Profiles metagenomes without assembly or read mapping.

## Pipeline (profile() in focus.py)
1. Load k-mer database (`load_database`) — tab-separated file, 8 taxonomy cols + k-mer counts
2. For each query FASTA/FASTQ:
   a. Count k-mers with Jellyfish (`count_kmers`) — calls `jellyfish count` + `jellyfish dump`
   b. Normalise counts (`normalise`) — divide by sum, raises RuntimeWarning if all-zero
   c. Run NNLS (`run_nnls`) — `scipy.optimize.nnls`, result normalised to sum=1
3. Aggregate by taxonomy level (`aggregate_level`) — positions 0–7 = Kingdom→Strain
4. Write CSV outputs (`write_results`) — abundances as percentages (×100)

## K-mer sizes
- `6` (default) or `7`
- Database files: `focus_app/db/k6` and `focus_app/db/k7`

## Taxonomy levels (in order)
Kingdom, Phylum, Class, Order, Family, Genus, Species, Strain

## Output files
- `{prefix}_All_levels.csv` — all 8 taxonomy columns
- `{prefix}_{Level}_tabular.csv` — one per level (8 files)
- Abundances are percentages (0–100)
- Compatible with STAMP for statistical analysis

## Database format
Tab-separated. Header: 8 taxonomy names + k-mer strings.
Each row: 8 taxonomy values + raw k-mer counts (normalised on load).

## Jellyfish requirement
Must be version 2.x. Checked at startup via `jellyfish count --version`.
On macOS: install via `brew install jellyfish` or `conda install -c bioconda jellyfish`.
