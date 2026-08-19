import pysam
import pytest

from hap_counter.bam_io import classify_allele, fetch_primary_pileup_reads, get_haplotype

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

    reads = [
        # HP=1 supporting ALT
        _make_read(header, "h1_alt", 0, 9, "6M", "AAAAAA", hp=1),
        # HP=1 supporting REF
        _make_read(header, "h1_ref", 0, 9, "6M", "GAAAAA", hp=1),
        # HP=2 supporting ALT
        _make_read(header, "h2_alt", 0, 9, "6M", "AAAAAA", hp=2),
        # HP=2 with a third allele (neither REF nor ALT)
        _make_read(header, "h2_other", 0, 9, "6M", "TAAAAA", hp=2),
        # No HP tag at all
        _make_read(header, "untagged", 0, 9, "6M", "AAAAAA", hp=None),
        # Deletion spanning the variant position
        _make_read(header, "deletion", 0, 8, "1M1D4M", "AAAAA", hp=1),
        # Secondary alignment - must be excluded even though it looks like ALT support
        _make_read(header, "secondary", 0x100, 9, "6M", "AAAAAA", hp=1),
        # Supplementary alignment - must be excluded
        _make_read(header, "supplementary", 0x800, 9, "6M", "AAAAAA", hp=1),
        # QC-failed alignment - must be excluded
        _make_read(header, "qcfail", 0x200, 9, "6M", "AAAAAA", hp=1),
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


def test_fetch_primary_pileup_reads_excludes_non_primary(synthetic_bam):
    reads = fetch_primary_pileup_reads(synthetic_bam, CHROM, VARIANT_POS)
    names = {pr.alignment.query_name for pr in reads}
    assert names == {"h1_alt", "h1_ref", "h2_alt", "h2_other", "untagged", "deletion"}


def test_classify_allele_and_get_haplotype(synthetic_bam):
    reads = {
        pr.alignment.query_name: pr
        for pr in fetch_primary_pileup_reads(synthetic_bam, CHROM, VARIANT_POS)
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
