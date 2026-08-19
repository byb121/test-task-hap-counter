import pysam
import pytest

from hap_counter.counting import count_haplotype_support_batches, _chunk_sites_by_chrom
from hap_counter.vcf_io import SnvSite

CHROM = "chr1"
POS1 = 10  # 1-based; REF=G ALT=A
POS2 = 15  # 1-based; REF=C ALT=T


def _make_read(header, name, reference_start, cigar, seq, hp=None):
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.flag = 0
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
def two_site_bam(tmp_path):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": 1000}]}
    )
    reads = [
        # Spans both sites: ALT at POS1 (offset 0), ALT at POS2 (offset 5). HP=1.
        _make_read(header, "both_alt", 9, "10M", "AAAAATAAAA", hp=1),
        # Spans both sites: REF at POS1, REF at POS2. HP=2.
        _make_read(header, "both_ref", 9, "10M", "GAAAACAAAA", hp=2),
        # Covers only POS1 (ends before POS2): ALT. HP=1.
        _make_read(header, "only_pos1", 9, "5M", "AAAAA", hp=1),
        # Covers only POS2 (starts after POS1): ALT, no HP tag.
        _make_read(header, "only_pos2_untagged", 14, "5M", "TAAAA", hp=None),
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


def test_count_haplotype_support_batches_gives_expected_counts(two_site_bam):
    sites = [
        SnvSite(CHROM, POS1, "G", "A"),
        SnvSite(CHROM, POS2, "C", "T"),
    ]

    batched_results = dict(count_haplotype_support_batches(two_site_bam, sites, batch_size=10))

    assert batched_results[sites[0]] == {
        "h1_ALT": 2,
        "h1_REF": 0,
        "h2_ALT": 0,
        "h2_REF": 1,
        "h1_other": 0,
        "h2_other": 0,
        "n_untagged": 0,
    }
    assert batched_results[sites[1]] == {
        "h1_ALT": 1,
        "h1_REF": 0,
        "h2_ALT": 0,
        "h2_REF": 1,
        "h1_other": 0,
        "h2_other": 0,
        "n_untagged": 1,
    }


def test_chunk_sites_by_chrom_splits_on_chrom_change():
    sites = [
        SnvSite("chr1", 1, "A", "G"),
        SnvSite("chr1", 2, "A", "G"),
        SnvSite("chr2", 3, "A", "G"),
        SnvSite("chr2", 4, "A", "G"),
    ]

    batches = list(_chunk_sites_by_chrom(sites, batch_size=10))

    assert [[s.pos for s in batch] for batch in batches] == [[1, 2], [3, 4]]
    assert all(len({s.chrom for s in batch}) == 1 for batch in batches)


def test_chunk_sites_by_chrom_splits_on_batch_size():
    sites = [SnvSite("chr1", i, "A", "G") for i in range(1, 6)]

    batches = list(_chunk_sites_by_chrom(sites, batch_size=2))

    assert [len(batch) for batch in batches] == [2, 2, 1]
