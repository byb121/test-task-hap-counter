# hap-counter

Computes per-haplotype REF/ALT read support for biallelic SNVs, given a haplotagged BAM
(reads carrying an `HP` tag of `1` or `2`) and a phased VCF.

## Installation

Requires Python >= 3.9. This repo uses [`uv`](https://docs.astral.sh/uv/) to manage the
`.env` virtual environment.

```bash
uv pip install -e .
```

This installs `hap-counter` as a console script into `.env/bin/`, along with its
dependencies (`pysam`, `click`).

## Usage

```bash
hap-counter --bam alignments.bam --vcf variants.vcf.gz --output support.tsv
```

Options:

| Option | Required | Description |
| --- | --- | --- |
| `--bam` / `-b` | yes | Coordinate-sorted, indexed BAM. Primary alignments may carry an `HP` tag (1 or 2). |
| `--vcf` / `-v` | yes | Phased VCF, optionally bgzipped/tabix-indexed. |
| `--output` / `-o` | yes | Path to write the output TSV to. |

### Example

```bash
hap-counter \
  --bam test_data/giab_2023.05.hg002.haplotagged.chr16_28000000_29000000.processed.30x.bam \
  --vcf test_data/giab_2023.05.hg002.wf_snp.chr16_28000000_29000000.vcf.gz \
  --output support.tsv
```

### Output

One row per biallelic SNV found in the VCF (indels, MNPs, multiallelic sites, and
non-`PASS` records are skipped). Columns:

| Column | Description |
| --- | --- |
| `chrom`, `pos` | Variant coordinates |
| `ref`, `alt` | REF/ALT alleles |
| `h1_REF`, `h1_ALT` | Primary alignments with `HP==1` supporting REF / ALT |
| `h2_REF`, `h2_ALT` | Primary alignments with `HP==2` supporting REF / ALT |
| `h1_other`, `h2_other` | Covering reads per haplotype whose base matched neither REF nor ALT |
| `n_untagged` | Covering reads with no (or an invalid) `HP` tag |

Only primary alignments are counted — secondary, supplementary, and QC-failed alignments
are excluded. Genotype/phasing in the VCF is not used for the counts; haplotype
assignment comes entirely from the BAM's `HP` tag.

## Performance

`pysam`/`htslib` pays a fixed per-call setup cost each time `bam.pileup()` is invoked,
even for a single-column, `truncate=True` query. An earlier version of this tool called
`pileup()` once per SNV site and took **~25s wall clock (~24ms/site)** on the bundled
`test_data/` (939 SNV sites, 35 MB BAM, chr16:28–29 Mb) — at that rate, a few million
SNVs (a typical whole-genome call set) would take on the order of many hours.

To fix this, `count_haplotype_support_batches` (`src/hap_counter/counting.py`) now
groups sites into batches of up to `DEFAULT_BATCH_SIZE` (150) **same-chromosome** sites
and issues one `pileup()` call per batch — spanning `min(pos)..max(pos)` of the batch —
via `fetch_primary_pileup_reads_for_positions` (`src/hap_counter/bam_io.py`), instead of
one call per site. This amortizes the fixed per-call cost across many sites at once.

On the same `test_data/`, the CLI now takes **~4s wall clock (~4ms/site)**, a ~6x
speedup. Benchmarking batch sizes from 1 to 939 (all sites in one call) showed wall time
drops steeply up to ~100 sites/batch and then flattens — going all the way to a single
batch only saves another ~10-15%, within run-to-run noise. `DEFAULT_BATCH_SIZE = 150`
was chosen from that plateau: it captures most of the achievable speedup while keeping
each `pileup()` call's genomic span (and read buffer) bounded and predictable, since
real VCFs can have far sparser/more irregular variant spacing than this test region
(gaps here reach 21.6kb; an unbounded single batch could span megabases elsewhere).

## Developer setup

Install with the `test` extra to also pull in `pytest`:

```bash
uv pip install -e ".[test]"
```

Run the test suite:

```bash
pytest
```

This runs unit tests for VCF filtering (`tests/test_vcf_io.py`) and BAM allele/haplotype
classification (`tests/test_bam_io.py`), plus an end-to-end test
(`tests/test_cli_e2e.py`) that runs the CLI against the bundled `test_data/` and checks
the row count and a couple of known site counts.
