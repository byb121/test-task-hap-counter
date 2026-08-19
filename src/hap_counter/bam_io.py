"""Reading BAM alignments at specific coordinates and classifying alleles."""

from typing import List, Optional

import pysam


def fetch_primary_pileup_reads(
    bam: pysam.AlignmentFile, chrom: str, pos: int
) -> List[pysam.PileupRead]:
    """Return primary-alignment pileup reads covering a 1-based position.

    Excludes secondary, supplementary, and QC-failed alignments. Uses
    min_base_quality=0 and flag_filter=0 to disable pysam's own hidden
    filtering defaults (a base-quality threshold, and a flag filter that
    excludes secondary/qcfail/duplicate but not supplementary), so all
    filtering is explicit here.
    """
    start = pos - 1
    end = pos
    reads: List[pysam.PileupRead] = []
    for column in bam.pileup(
        chrom,
        start,
        end,
        truncate=True,
        min_base_quality=0,
        flag_filter=0,
        ignore_overlaps=False,
    ):
        if column.pos != start:
            continue
        for pileup_read in column.pileups:
            alignment = pileup_read.alignment
            if alignment.is_secondary or alignment.is_supplementary or alignment.is_qcfail:
                continue
            reads.append(pileup_read)
    return reads


def classify_allele(pileup_read: pysam.PileupRead, ref: str, alt: str) -> Optional[str]:
    """Classify the base a read carries at the pileup column as REF, ALT, or None.

    Returns None when the read has a deletion/reference-skip at the site,
    or its base matches neither the REF nor the ALT allele.
    """
    if pileup_read.is_del or pileup_read.is_refskip:
        return None

    alignment = pileup_read.alignment
    base = alignment.query_sequence[pileup_read.query_position].upper()
    if base == ref:
        return "REF"
    if base == alt:
        return "ALT"
    return None


def get_haplotype(pileup_read: pysam.PileupRead) -> Optional[int]:
    """Return 1 or 2 if the read has a valid HP tag, else None."""
    alignment = pileup_read.alignment
    if not alignment.has_tag("HP"):
        return None
    hp = alignment.get_tag("HP")
    if hp in (1, 2):
        return hp
    return None
