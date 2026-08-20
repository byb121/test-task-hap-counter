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
dependencies (`pysam`, `click`, `matplotlib`).

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
| `--sample` / `-s` | no | VCF sample to genotype against. Required if the VCF has more than one sample (the CLI errors out listing the available sample names); defaults to the VCF's only sample otherwise. |
| `--min-genotyping-hap-reads` / `-m` | no | Minimum REF+ALT reads a haplotype needs at a site before `genotyping_from_input_bam` attempts a call. Default `10`. |
| `--genotype-allele-fraction` / `-f` | no | Fraction of a haplotype's REF+ALT reads that must support one allele for `genotyping_from_input_bam` to call it. Default `0.8`. |
| `--discrepancy-plot` / `-d` | no | If set, also computes per-(chromosome, haplotype) discrepant-vote statistics across confidently-genotyped SNVs, prints the aggregate summary to stdout, and saves a histogram of discrepancy fractions to this path. |

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
| `genotyping_from_input_vcf` | The selected sample's VCF genotype (`GT`) at this site, passed through verbatim — `\|`-separated if the record is phased, `/`-separated if not, `.` for any missing allele |
| `genotyping_from_input_bam` | Per-haplotype genotype call derived from read support: for each haplotype, if its REF+ALT reads (excluding `other`) total at least `--min-genotyping-hap-reads` and one allele reaches `--genotype-allele-fraction` of them, that haplotype is called `0` (REF) or `1` (ALT); otherwise `.` (no call). Formatted `H1\|H2`, e.g. `0\|1` |

Only primary alignments are counted — secondary, supplementary, and QC-failed alignments
are excluded. Haplotype assignment for the counting columns comes entirely from the
BAM's `HP` tag, not from VCF phase.

Reading the VCF's genotype also compares the BAM's sample name (from its `@RG`/`SM`
tags) against the VCF sample being used, and prints a warning if they don't match — or
if the BAM has no sample name in its header at all (as with the bundled `test_data/`
BAM, which has no `@RG` lines).

### Discrepancy analysis (`--discrepancy-plot`)

When `--discrepancy-plot PATH` is given, for every haplotype that got a confident call
in `genotyping_from_input_bam` (i.e. not `.`), the reads that voted for the *other*
allele are counted as "discrepant votes". These are aggregated per
`(chromosome, haplotype)` and:

* printed to stdout as a summary (total discrepant votes, total votes, site count, and
  overall discrepancy rate per group), and
* plotted as a histogram of per-site discrepancy fractions (`PATH`), with one overlaid
  series per `(chromosome, haplotype)` group and the same summary text embedded below
  the plot.

Sites where a haplotype got no confident call are excluded from this analysis — it
measures disagreement with an actual genotype call, not general site ambiguity.

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

This runs unit tests for VCF filtering/genotyping (`tests/test_vcf_io.py`), BAM
allele/haplotype/sample-name handling (`tests/test_bam_io.py`), per-haplotype genotype
calling (`tests/test_counting.py`), discrepant-vote aggregation
(`tests/test_discrepancy.py`), the histogram plot (`tests/test_plotting.py`), and CLI
helpers (`tests/test_cli.py`), plus end-to-end tests (`tests/test_cli_e2e.py`) that run
the CLI against the bundled `test_data/` and check row counts, known site counts, the
new genotyping columns, the sample-name-mismatch warning, and `--discrepancy-plot`.
