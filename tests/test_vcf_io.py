import textwrap

import pytest

from hap_counter.vcf_io import (
    SnvSite,
    VcfSampleSelectionError,
    format_genotype,
    get_vcf_sample_names,
    read_biallelic_snvs,
    resolve_vcf_sample,
)

VCF_HEADER = textwrap.dedent(
    """\
    ##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##FILTER=<ID=LowQual,Description="Low quality variant">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##contig=<ID=chr1,length=1000000>
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
    """
)

TWO_SAMPLE_VCF_HEADER = textwrap.dedent(
    """\
    ##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##contig=<ID=chr1,length=1000000>
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_a\tsample_b
    """
)


def write_vcf(tmp_path, records, header=VCF_HEADER, name="variants.vcf"):
    path = tmp_path / name
    path.write_text(header + "\n".join(records) + ("\n" if records else ""))
    return str(path)


# --- format_genotype ---


def test_format_genotype_phased():
    assert format_genotype((0, 1), True) == "0|1"


def test_format_genotype_unphased():
    assert format_genotype((1, 0), False) == "1/0"


def test_format_genotype_hom_ref():
    assert format_genotype((0, 0), True) == "0|0"


def test_format_genotype_hom_alt():
    assert format_genotype((1, 1), False) == "1/1"


def test_format_genotype_missing_allele():
    assert format_genotype((None, 1), True) == ".|1"


def test_format_genotype_fully_missing():
    assert format_genotype((None, None), False) == "./."


def test_format_genotype_haploid():
    assert format_genotype((1,), True) == "1"


# --- resolve_vcf_sample ---


def test_resolve_vcf_sample_defaults_to_only_sample(tmp_path):
    vcf_path = write_vcf(tmp_path, ["chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t1|0"])
    assert resolve_vcf_sample(vcf_path, None) == "sample"


def test_resolve_vcf_sample_accepts_explicit_matching_sample(tmp_path):
    vcf_path = write_vcf(tmp_path, ["chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t1|0"])
    assert resolve_vcf_sample(vcf_path, "sample") == "sample"


def test_resolve_vcf_sample_raises_for_unknown_explicit_sample(tmp_path):
    vcf_path = write_vcf(tmp_path, ["chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t1|0"])
    with pytest.raises(VcfSampleSelectionError) as exc_info:
        resolve_vcf_sample(vcf_path, "nope")
    message = str(exc_info.value)
    assert "nope" in message
    assert "sample" in message


def test_resolve_vcf_sample_raises_when_multiple_samples_and_none_given(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t1|0\t0|1"],
        header=TWO_SAMPLE_VCF_HEADER,
    )
    with pytest.raises(VcfSampleSelectionError) as exc_info:
        resolve_vcf_sample(vcf_path, None)
    message = str(exc_info.value)
    assert "sample_a" in message
    assert "sample_b" in message


def test_resolve_vcf_sample_picks_explicit_sample_among_multiple(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t1|0\t0|1"],
        header=TWO_SAMPLE_VCF_HEADER,
    )
    assert resolve_vcf_sample(vcf_path, "sample_b") == "sample_b"


# --- get_vcf_sample_names ---


def test_get_vcf_sample_names_single(tmp_path):
    vcf_path = write_vcf(tmp_path, ["chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t1|0"])
    assert get_vcf_sample_names(vcf_path) == ["sample"]


def test_get_vcf_sample_names_multiple(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t1|0\t0|1"],
        header=TWO_SAMPLE_VCF_HEADER,
    )
    assert get_vcf_sample_names(vcf_path) == ["sample_a", "sample_b"]


# --- read_biallelic_snvs (genotype field) ---


def test_keeps_pass_biallelic_snv(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t1|0"],
    )
    assert list(read_biallelic_snvs(vcf_path, "sample")) == [
        SnvSite("chr1", 100, "G", "A", "1|0")
    ]


def test_skips_non_pass_filter(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tG\tA\t30\tLowQual\t.\tGT\t1|0"],
    )
    assert list(read_biallelic_snvs(vcf_path, "sample")) == []


def test_skips_indel(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tCAGAAT\tC\t30\tPASS\t.\tGT\t0/1"],
    )
    assert list(read_biallelic_snvs(vcf_path, "sample")) == []


def test_skips_multiallelic(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tG\tA,T\t30\tPASS\t.\tGT\t1/2"],
    )
    assert list(read_biallelic_snvs(vcf_path, "sample")) == []


def test_skips_symbolic_allele(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tG\t*\t30\tPASS\t.\tGT\t0/1"],
    )
    assert list(read_biallelic_snvs(vcf_path, "sample")) == []


def test_keeps_unphased_and_hom_alt_sites(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        [
            "chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t0/1",
            "chr1\t200\t.\tC\tT\t30\tPASS\t.\tGT\t1/1",
        ],
    )
    assert list(read_biallelic_snvs(vcf_path, "sample")) == [
        SnvSite("chr1", 100, "G", "A", "0/1"),
        SnvSite("chr1", 200, "C", "T", "1/1"),
    ]


def test_keeps_only_pass_biallelic_snvs_from_mixed_records(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        [
            # PASS SNVs (phased, unphased, hom-alt) - 4 total, kept
            "chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t1|0",
            "chr1\t200\t.\tC\tT\t30\tPASS\t.\tGT\t0/1",
            "chr1\t300\t.\tA\tG\t30\tPASS\t.\tGT\t1/1",
            "chr1\t400\t.\tT\tC\t30\tPASS\t.\tGT\t1|0",
            # non-PASS SNV, skipped
            "chr1\t500\t.\tG\tA\t30\tLowQual\t.\tGT\t0/1",
            # PASS indels - 3 total, skipped
            "chr1\t600\t.\tCAGAAT\tC\t30\tPASS\t.\tGT\t0/1",
            "chr1\t700\t.\tC\tCAG\t30\tPASS\t.\tGT\t0/1",
            "chr1\t750\t.\tGAT\tG\t30\tPASS\t.\tGT\t1/1",
            # non-PASS indel, skipped
            "chr1\t800\t.\tA\tATG\t30\tLowQual\t.\tGT\t0/1",
            # PASS MNVs - 3 total, skipped
            "chr1\t900\t.\tAC\tGT\t30\tPASS\t.\tGT\t0/1",
            "chr1\t950\t.\tTAG\tCGA\t30\tPASS\t.\tGT\t1/1",
            "chr1\t1000\t.\tGGA\tCCT\t30\tPASS\t.\tGT\t0/1",
            # non-PASS MNV, skipped
            "chr1\t1050\t.\tAC\tGT\t30\tLowQual\t.\tGT\t0/1",
            # PASS multiallelic SNV, skipped
            "chr1\t1100\t.\tG\tA,T\t30\tPASS\t.\tGT\t1/2",
            # PASS symbolic allele, skipped
            "chr1\t1150\t.\tG\t*\t30\tPASS\t.\tGT\t0/1",
        ],
    )
    assert list(read_biallelic_snvs(vcf_path, "sample")) == [
        SnvSite("chr1", 100, "G", "A", "1|0"),
        SnvSite("chr1", 200, "C", "T", "0/1"),
        SnvSite("chr1", 300, "A", "G", "1/1"),
        SnvSite("chr1", 400, "T", "C", "1|0"),
    ]


def test_read_biallelic_snvs_uses_requested_sample_among_multiple(tmp_path):
    vcf_path = write_vcf(
        tmp_path,
        ["chr1\t100\t.\tG\tA\t30\tPASS\t.\tGT\t1|0\t0|1"],
        header=TWO_SAMPLE_VCF_HEADER,
    )
    assert list(read_biallelic_snvs(vcf_path, "sample_b")) == [
        SnvSite("chr1", 100, "G", "A", "0|1")
    ]
