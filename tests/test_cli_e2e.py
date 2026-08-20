import csv
import textwrap
from pathlib import Path

import pytest
from click.testing import CliRunner

from hap_counter.cli import main

TEST_DATA = Path(__file__).resolve().parent.parent / "test_data"
BAM_PATH = TEST_DATA / "giab_2023.05.hg002.haplotagged.chr16_28000000_29000000.processed.30x.bam"
VCF_PATH = TEST_DATA / "giab_2023.05.hg002.wf_snp.chr16_28000000_29000000.vcf.gz"

TWO_SAMPLE_VCF_HEADER = textwrap.dedent(
    """\
    ##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##contig=<ID=chr16,length=90338345>
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_a\tsample_b
    """
)


def run_cli(tmp_path, extra_args=()):
    output_path = tmp_path / "support.tsv"
    runner = CliRunner()
    args = ["--bam", str(BAM_PATH), "--vcf", str(VCF_PATH), "--output", str(output_path)]
    args.extend(extra_args)
    with pytest.warns(UserWarning, match="no sample name"):
        result = runner.invoke(main, args)
    assert result.exit_code == 0, result.output
    with open(output_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return result, rows


def test_output_has_expected_header_and_row_count(tmp_path):
    _, rows = run_cli(tmp_path)
    # 939 biallelic PASS SNVs in the test VCF (indels/multiallelic/non-PASS excluded).
    assert len(rows) == 939


def test_known_site_counts(tmp_path):
    _, rows = run_cli(tmp_path)
    by_pos = {row["pos"]: row for row in rows}

    # Cross-checked independently against a standalone pysam pileup script
    # using pysam's default min_base_quality (13) and flag_filter.
    row = by_pos["28001381"]
    assert row["chrom"] == "chr16"
    assert row["ref"] == "G"
    assert row["alt"] == "A"
    assert row["h1_ALT"] == "4"
    assert row["h1_REF"] == "11"
    assert row["h2_ALT"] == "11"
    assert row["h2_REF"] == "3"

    row = by_pos["28002344"]
    assert row["h1_ALT"] == "4"
    assert row["h1_REF"] == "10"
    assert row["h2_ALT"] == "9"
    assert row["h2_REF"] == "3"
    assert row["h2_other"] == "0"


def test_output_has_genotyping_columns(tmp_path):
    _, rows = run_cli(tmp_path)
    assert "genotyping_from_input_vcf" in rows[0]
    assert "genotyping_from_input_bam" in rows[0]

    row = {r["pos"]: r for r in rows}["28001381"]
    # h1_REF=11,h1_ALT=4 (frac 11/15=0.733) and h2_REF=3,h2_ALT=11 (frac 11/14=0.786):
    # neither haplotype clears the default 0.8 threshold.
    assert row["genotyping_from_input_bam"] == ".|."
    # VCF genotype is passed through verbatim (phased/unphased as encoded).
    assert row["genotyping_from_input_vcf"] in {"1|0", "0|1", "0/1", "1/0"}


def test_lowering_allele_fraction_threshold_produces_a_call(tmp_path):
    _, rows = run_cli(tmp_path, extra_args=["--genotype-allele-fraction", "0.7"])
    row = {r["pos"]: r for r in rows}["28001381"]
    # h1 REF frac 11/15=0.733 >= 0.7 -> "0"; h2 ALT frac 11/14=0.786 >= 0.7 -> "1".
    assert row["genotyping_from_input_bam"] == "0|1"


def test_raising_min_reads_forces_every_site_to_no_call(tmp_path):
    _, rows = run_cli(tmp_path, extra_args=["--min-genotyping-hap-reads", "1000000"])
    assert all(row["genotyping_from_input_bam"] == ".|." for row in rows)


def test_explicit_sample_matches_default_single_sample(tmp_path):
    _, default_rows = run_cli(tmp_path)
    _, explicit_rows = run_cli(tmp_path, extra_args=["--sample", "hg002"])
    assert [r["genotyping_from_input_vcf"] for r in default_rows] == [
        r["genotyping_from_input_vcf"] for r in explicit_rows
    ]


def test_unknown_sample_errors_out_listing_available_samples(tmp_path):
    output_path = tmp_path / "support.tsv"
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "--bam", str(BAM_PATH),
            "--vcf", str(VCF_PATH),
            "--output", str(output_path),
            "--sample", "nonexistent_sample",
        ],
    )
    assert result.exit_code != 0
    assert "nonexistent_sample" in result.output
    assert "hg002" in result.output


def test_multi_sample_vcf_without_sample_flag_errors_out(tmp_path):
    vcf_path = tmp_path / "multi_sample.vcf"
    vcf_path.write_text(
        TWO_SAMPLE_VCF_HEADER + "chr16\t28001381\t.\tG\tA\t30\tPASS\t.\tGT\t1|0\t0|1\n"
    )
    output_path = tmp_path / "support.tsv"
    runner = CliRunner()
    result = runner.invoke(
        main,
        ["--bam", str(BAM_PATH), "--vcf", str(vcf_path), "--output", str(output_path)],
    )
    assert result.exit_code != 0
    assert "sample_a" in result.output
    assert "sample_b" in result.output
    assert "--sample" in result.output


def test_discrepancy_plot_writes_file_and_prints_summary(tmp_path):
    plot_path = tmp_path / "discrepancy.png"
    result, _ = run_cli(tmp_path, extra_args=["--discrepancy-plot", str(plot_path)])
    assert plot_path.exists()
    assert plot_path.stat().st_size > 0
    assert "discrepant" in result.output.lower()


def test_omitting_discrepancy_plot_prints_no_summary(tmp_path):
    result, _ = run_cli(tmp_path)
    assert "discrepant" not in result.output.lower()
