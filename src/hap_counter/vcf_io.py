"""Reading and filtering of VCF records."""

from typing import Iterator, NamedTuple

import pysam

_VALID_BASES = frozenset("ACGT")


class SnvSite(NamedTuple):
    chrom: str
    pos: int
    ref: str
    alt: str


def read_biallelic_snvs(vcf_path: str) -> Iterator[SnvSite]:
    """Yield biallelic SNV sites from a (optionally gzipped) VCF.

    A record is kept iff its FILTER is exactly PASS, it has a single-base
    REF, exactly one single-base ALT allele, and both alleles are plain
    A/C/G/T (excludes symbolic alleles such as ``*``). Indels, MNPs, and
    multiallelic records are skipped. Genotype/phasing is intentionally
    not used to select sites - haplotype assignment comes from the BAM's
    HP tag, not from VCF phase.
    """
    with pysam.VariantFile(vcf_path) as vcf_file:
        for record in vcf_file:
            if list(record.filter.keys()) != ["PASS"]:
                continue

            alts = record.alts
            if alts is None or len(alts) != 1:
                continue

            ref, alt = record.ref, alts[0]
            if len(ref) != 1 or len(alt) != 1:
                continue
            if ref not in _VALID_BASES or alt not in _VALID_BASES:
                continue

            yield SnvSite(chrom=record.chrom, pos=record.pos, ref=ref, alt=alt)
