"""Per-SNV haplotype support counting."""

from typing import Dict

import pysam

from hap_counter.bam_io import classify_allele, fetch_primary_pileup_reads, get_haplotype

COUNT_FIELDS = (
    "h1_ALT",
    "h1_REF",
    "h2_ALT",
    "h2_REF",
    "h1_other",
    "h2_other",
    "n_untagged",
)


def count_haplotype_support(
    bam: pysam.AlignmentFile, chrom: str, pos: int, ref: str, alt: str
) -> Dict[str, int]:
    """Count primary-alignment REF/ALT support per haplotype at one SNV site."""
    counts = {field: 0 for field in COUNT_FIELDS}

    for pileup_read in fetch_primary_pileup_reads(bam, chrom, pos):
        haplotype = get_haplotype(pileup_read)
        if haplotype is None:
            counts["n_untagged"] += 1
            continue

        allele = classify_allele(pileup_read, ref, alt)
        if allele is None:
            counts[f"h{haplotype}_other"] += 1
        else:
            counts[f"h{haplotype}_{allele}"] += 1

    return counts
