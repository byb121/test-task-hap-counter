"""Per-SNV haplotype support counting."""

from typing import Dict, Iterable, Iterator, List, Tuple

import pysam

from hap_counter.bam_io import (
    classify_allele,
    fetch_primary_pileup_reads_for_positions,
    get_haplotype,
)
from hap_counter.vcf_io import SnvSite

COUNT_FIELDS = (
    "h1_REF",
    "h1_ALT",
    "h2_REF",
    "h2_ALT",
    "h1_other",
    "h2_other",
    "n_untagged",
)

# Sweet spot found by benchmarking on test_data/: wall time drops steeply up
# to ~100 sites/batch, then flattens (see README2's Performance section).
# Bounded rather than "everything in one batch" so a single pileup() call's
# genomic span/read buffer stays predictable on real, sparser VCFs.
DEFAULT_BATCH_SIZE = 150


def _count_from_reads(reads, ref: str, alt: str) -> Dict[str, int]:
    counts = {field: 0 for field in COUNT_FIELDS}
    for pileup_read in reads:
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


def _chunk_sites_by_chrom(
    sites: Iterable[SnvSite], batch_size: int
) -> Iterator[List[SnvSite]]:
    """Group sites into batches of at most batch_size, never spanning chromosomes."""
    batch: List[SnvSite] = []
    for site in sites:
        if batch and (site.chrom != batch[0].chrom or len(batch) >= batch_size):
            yield batch
            batch = []
        batch.append(site)
    if batch:
        yield batch


def count_haplotype_support_batches(
    bam: pysam.AlignmentFile,
    sites: Iterable[SnvSite],
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> Iterator[Tuple[SnvSite, Dict[str, int]]]:
    """Yield (site, counts) for each site, computed via one pileup() call per
    batch of up to batch_size same-chromosome sites instead of one call per site.

    Sites must be sorted by position within each chromosome (as VCF/read_biallelic_snvs
    yields them); batches never span a chromosome boundary.
    """
    for batch in _chunk_sites_by_chrom(sites, batch_size):
        reads_by_pos = fetch_primary_pileup_reads_for_positions(
            bam, batch[0].chrom, [site.pos for site in batch]
        )
        for site in batch:
            yield site, _count_from_reads(reads_by_pos.get(site.pos, []), site.ref, site.alt)
