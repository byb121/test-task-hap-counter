import pysam
import pytest

from hap_counter.bam_io import (classify_allele,
        fetch_primary_pileup_reads_for_positions, get_bam_sample_names, get_haplotype)

CHROM = "chr1"
VARIANT_POS = 10  # 1-based; 0-based reference position is 9
REF = "G"
ALT = "A"


def _make_read(header, name, flag, reference_start, cigar, seq, hp=None):
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.flag = flag
    read.reference_id = 0
    read.reference_start = reference_start
    read.mapping_quality = 60
    read.cigarstring = cigar
    read.query_sequence = seq
    read.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    if hp is not None:
        read.set_tag("HP", hp, value_type="i")
    return read


@pytest.fixture
def synthetic_bam(tmp_path):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": 1000}]}
    )

    # Reads deliberately vary in length (5-20bp) and reference_start so pileup
    # logic is exercised against reads that reach the variant position from
    # different offsets, not just uniform 6bp reads starting at the variant.
    reads = [
        # HP=1 supporting ALT (length 11, starts 7bp before the variant)
        _make_read(header, "h1_alt", 0, 2, "11M", "A" * 11, hp=1),
        # HP=1 supporting REF (length 6, starts at the variant)
        _make_read(header, "h1_ref", 0, 9, "6M", "GAAAAA", hp=1),
        # HP=2 supporting ALT (length 20, starts 9bp before the variant)
        _make_read(header, "h2_alt", 0, 0, "20M", "A" * 20, hp=2),
        # HP=2 with a third allele (neither REF nor ALT); length 9
        _make_read(header, "h2_other", 0, 5, "9M", "AAAATAAAA", hp=2),
        # No HP tag at all; length 16
        _make_read(header, "untagged", 0, 1, "16M", "A" * 16, hp=None),
        # Deletion spanning the variant position; length 7
        _make_read(header, "deletion", 0, 7, "2M1D5M", "A" * 7, hp=1),
        # Secondary alignment - must be excluded even though it looks like ALT support
        _make_read(header, "secondary", 0x100, 6, "5M", "A" * 5, hp=1),
        # Supplementary alignment - must be excluded; length 18
        _make_read(header, "supplementary", 0x800, 3, "18M", "A" * 18, hp=1),
        # QC-failed alignment - must be excluded; length 8
        _make_read(header, "qcfail", 0x200, 4, "8M", "A" * 8, hp=1),
    ]

    unsorted_path = tmp_path / "unsorted.bam"
    with pysam.AlignmentFile(str(unsorted_path), "wb", header=header) as out:
        for read in reads:
            out.write(read)

    sorted_path = tmp_path / "reads.bam"
    pysam.sort("-o", str(sorted_path), str(unsorted_path))
    pysam.index(str(sorted_path))

    with pysam.AlignmentFile(str(sorted_path), "rb") as bam:
        yield bam


def test_fetch_primary_pileup_reads_for_positions_excludes_non_primary(synthetic_bam):
    result = fetch_primary_pileup_reads_for_positions(synthetic_bam, CHROM, [VARIANT_POS])
    names = {pr.alignment.query_name for pr in result[VARIANT_POS]}
    assert names == {"h1_alt", "h1_ref", "h2_alt", "h2_other", "untagged", "deletion"}


def test_fetch_primary_pileup_reads_for_positions_omits_uncovered_positions(synthetic_bam):
    result = fetch_primary_pileup_reads_for_positions(synthetic_bam, CHROM, [VARIANT_POS, 999])
    assert 999 not in result


def test_classify_allele_and_get_haplotype(synthetic_bam):
    reads = {
        pr.alignment.query_name: pr
        for pr in fetch_primary_pileup_reads_for_positions(synthetic_bam, CHROM, [VARIANT_POS])[
            VARIANT_POS
        ]
    }

    assert classify_allele(reads["h1_alt"], REF, ALT) == "ALT"
    assert get_haplotype(reads["h1_alt"]) == 1

    assert classify_allele(reads["h1_ref"], REF, ALT) == "REF"
    assert get_haplotype(reads["h1_ref"]) == 1

    assert classify_allele(reads["h2_alt"], REF, ALT) == "ALT"
    assert get_haplotype(reads["h2_alt"]) == 2

    assert classify_allele(reads["h2_other"], REF, ALT) is None
    assert get_haplotype(reads["h2_other"]) == 2

    assert get_haplotype(reads["untagged"]) is None

    assert classify_allele(reads["deletion"], REF, ALT) is None


# --- get_bam_sample_names ---


def test_get_bam_sample_names_single_rg():
    header = {"RG": [{"ID": "rg1", "SM": "sampleA"}]}
    assert get_bam_sample_names(header) == ["sampleA"]


def test_get_bam_sample_names_dedups_same_sample_across_rgs():
    header = {"RG": [{"ID": "rg1", "SM": "sampleA"}, {"ID": "rg2", "SM": "sampleA"}]}
    assert get_bam_sample_names(header) == ["sampleA"]


def test_get_bam_sample_names_multiple_distinct_samples():
    header = {"RG": [{"ID": "rg1", "SM": "sampleA"}, {"ID": "rg2", "SM": "sampleB"}]}
    assert get_bam_sample_names(header) == ["sampleA", "sampleB"]


def test_get_bam_sample_names_no_rg_at_all():
    header = {"HD": {"VN": "1.6"}}
    assert get_bam_sample_names(header) == []


def test_get_bam_sample_names_rg_without_sm():
    header = {"RG": [{"ID": "rg1"}]}
    assert get_bam_sample_names(header) == []
