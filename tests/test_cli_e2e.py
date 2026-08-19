import csv
from pathlib import Path

from click.testing import CliRunner

from hap_counter.cli import main

TEST_DATA = Path(__file__).resolve().parent.parent / "test_data"
BAM_PATH = TEST_DATA / "giab_2023.05.hg002.haplotagged.chr16_28000000_29000000.processed.30x.bam"
VCF_PATH = TEST_DATA / "giab_2023.05.hg002.wf_snp.chr16_28000000_29000000.vcf.gz"


def run_cli(tmp_path):
    output_path = tmp_path / "support.tsv"
    runner = CliRunner()
    result = runner.invoke(
        main,
        ["--bam", str(BAM_PATH), "--vcf", str(VCF_PATH), "--output", str(output_path)],
    )
    assert result.exit_code == 0, result.output
    with open(output_path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_output_has_expected_header_and_row_count(tmp_path):
    rows = run_cli(tmp_path)
    # 939 biallelic PASS SNVs in the test VCF (indels/multiallelic/non-PASS excluded).
    assert len(rows) == 939


def test_known_site_counts(tmp_path):
    rows = run_cli(tmp_path)
    by_pos = {row["pos"]: row for row in rows}

    # Cross-checked independently against a standalone pysam pileup script
    # using the same min_base_quality=0/flag_filter=0 parameters.
    row = by_pos["28001381"]
    assert row["chrom"] == "chr16"
    assert row["ref"] == "G"
    assert row["alt"] == "A"
    assert row["h1_ALT"] == "4"
    assert row["h1_REF"] == "11"
    assert row["h2_ALT"] == "12"
    assert row["h2_REF"] == "3"

    row = by_pos["28002344"]
    assert row["h1_ALT"] == "5"
    assert row["h1_REF"] == "11"
    assert row["h2_ALT"] == "12"
    assert row["h2_REF"] == "3"
    assert row["h2_other"] == "1"
