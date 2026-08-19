"""Reading BAM alignments at specific coordinates and classifying alleles."""

from typing import Dict, List, Optional

import pysam


def fetch_primary_pileup_reads_for_positions(
    bam: pysam.AlignmentFile, chrom: str, positions: List[int]
) -> Dict[int, List[pysam.PileupRead]]:
    """Return primary-alignment pileup reads at several 1-based positions on one
    chromosome, via a single pileup() call spanning min(positions)..max(positions).

    Amortizes htslib's per-call pileup setup cost across a batch of nearby
    positions instead of paying it once per position. A position absent from
    the returned dict had zero covering reads. Uses pysam's default
    min_base_quality (13): a read whose base quality at a given position
    falls below that threshold is silently excluded from that position's
    pileup entirely (it won't be counted as REF/ALT/other/untagged there).
    Excludes secondary, supplementary, QC-failed, and duplicate alignments:
    pysam's default flag_filter (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL |
    BAM_FDUP) excludes unmapped/secondary/QC-failed/duplicate reads;
    supplementary alignments aren't covered by that default, so they're
    filtered explicitly here.
    """
    if not positions:
        return {}

    start = min(positions) - 1
    end = max(positions)
    wanted = set(positions)
    reads_by_pos: Dict[int, List[pysam.PileupRead]] = {}
    for column in bam.pileup(
        chrom,
        start,
        end,
        truncate=True,
        ignore_overlaps=False,
    ):
        pos = column.pos + 1
        if pos not in wanted:
            continue
        reads_by_pos[pos] = [
            pileup_read
            for pileup_read in column.pileups
            if not pileup_read.alignment.is_supplementary
        ]
    return reads_by_pos


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
