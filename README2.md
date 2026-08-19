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
| `h1_ALT`, `h1_REF` | Primary alignments with `HP==1` supporting ALT / REF |
| `h2_ALT`, `h2_REF` | Primary alignments with `HP==2` supporting ALT / REF |
| `h1_other`, `h2_other` | Covering reads per haplotype whose base matched neither REF nor ALT |
| `n_untagged` | Covering reads with no (or an invalid) `HP` tag |

Only primary alignments are counted — secondary, supplementary, and QC-failed alignments
are excluded. Genotype/phasing in the VCF is not used for the counts; haplotype
assignment comes entirely from the BAM's `HP` tag.

## Performance

On the bundled `test_data/` (939 SNV sites, 35 MB BAM, chr16:28–29 Mb), the CLI takes
**~25s wall clock (~24ms per SNV site)**. Profiling shows this time is almost entirely
inside `bam.pileup()` itself (our Python-level filtering/classification logic accounts
for well under 10% of it) — `pysam`/`htslib` pays a fixed per-call setup cost each time
`pileup()` is invoked, even for a single-column, `truncate=True` query, and we currently
call it once per SNV site (`fetch_primary_pileup_reads` in `src/hap_counter/bam_io.py`).

This is fine at the scale of this test region, but it does not scale to genome-wide VCFs:
at ~24ms/site, a few million SNVs (a typical whole-genome call set) would take on the
order of many hours.

**Potential improvement (next iteration):** replace the per-site `pileup()` calls with a
single coordinated linear pass — walk the coordinate-sorted BAM and the coordinate-sorted
VCF together in lockstep (e.g. advance a read-window buffer as variant positions are
consumed), rather than doing a fresh indexed pileup lookup per variant. This changes the
BAM-reading function from "coordinates in, reads out per call" to a streaming/batched
shape, so it wasn't pursued in this iteration in order to keep that function's interface
as originally specified.

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
