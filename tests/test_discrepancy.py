from hap_counter.discrepancy import DiscrepancyAccumulator


def test_add_skips_no_call_sites():
    acc = DiscrepancyAccumulator()
    acc.add(chrom="chr1", hap=1, ref_count=5, alt_count=5, call=".")
    assert acc.fractions_by_group() == {}


def test_add_records_fraction_for_confident_call():
    acc = DiscrepancyAccumulator()
    acc.add(chrom="chr1", hap=1, ref_count=8, alt_count=2, call="0")
    assert acc.fractions_by_group() == {("chr1", 1): [0.2]}


def test_fractions_by_group_groups_by_chrom_and_hap():
    acc = DiscrepancyAccumulator()
    acc.add(chrom="chr1", hap=1, ref_count=9, alt_count=1, call="0")
    acc.add(chrom="chr1", hap=2, ref_count=1, alt_count=9, call="1")
    acc.add(chrom="chr2", hap=1, ref_count=10, alt_count=0, call="0")

    groups = acc.fractions_by_group()
    assert groups[("chr1", 1)] == [0.1]
    assert groups[("chr1", 2)] == [0.1]
    assert groups[("chr2", 1)] == [0.0]


def test_summary_totals_and_overall_rate():
    acc = DiscrepancyAccumulator()
    acc.add(chrom="chr1", hap=1, ref_count=9, alt_count=1, call="0")
    acc.add(chrom="chr1", hap=1, ref_count=7, alt_count=3, call="0")

    summaries = {(s.chrom, s.hap): s for s in acc.summary()}
    group = summaries[("chr1", 1)]
    assert group.n_sites == 2
    assert group.discrepant_votes == 1 + 3
    assert group.total_votes == 10 + 10
    assert group.overall_rate == (1 + 3) / (10 + 10)


def test_summary_excludes_no_call_sites():
    acc = DiscrepancyAccumulator()
    acc.add(chrom="chr1", hap=1, ref_count=5, alt_count=5, call=".")
    assert acc.summary() == []


def test_format_summary_text_includes_group_numbers():
    acc = DiscrepancyAccumulator()
    acc.add(chrom="chr1", hap=1, ref_count=9, alt_count=1, call="0")

    text = acc.format_summary_text()
    assert "chr1" in text
    assert "1" in text  # discrepant vote count and/or haplotype number
    assert "10" in text  # total votes


def test_format_summary_text_for_no_data():
    acc = DiscrepancyAccumulator()
    text = acc.format_summary_text()
    assert isinstance(text, str)
    assert text != ""
