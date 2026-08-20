"""Reading and filtering of VCF records."""

from typing import Iterator, List, NamedTuple, Optional

import pysam

_VALID_BASES = frozenset("ACGT")


class SnvSite(NamedTuple):
    chrom: str
    pos: int
    ref: str
    alt: str
    genotype: str


class VcfSampleSelectionError(Exception):
    """Raised when the sample to use for genotyping can't be unambiguously resolved."""


def format_genotype(gt: tuple, phased: bool) -> str:
    """Format a pysam GT tuple as a genotype string, e.g. (0, 1) -> "0|1".

    Each allele index is rendered as-is ("0"/"1"/...), missing alleles
    (``None``) as ".". Alleles are joined with "|" if phased else "/". A
    single-allele (haploid) GT renders as just that one allele, no separator.
    """
    alleles = ["." if allele is None else str(allele) for allele in gt]
    separator = "|" if phased else "/"
    return separator.join(alleles)


def get_vcf_sample_names(vcf_path: str) -> List[str]:
    """Return the list of sample names in the VCF's header, in header order."""
    with pysam.VariantFile(vcf_path) as vcf_file:
        return list(vcf_file.header.samples)


def resolve_vcf_sample(vcf_path: str, sample: Optional[str]) -> str:
    """Resolve which VCF sample to use for genotyping.

    If ``sample`` is given, it must be present in the VCF's header. If
    ``sample`` is None, the VCF must have exactly one sample. Raises
    VcfSampleSelectionError (listing the available sample names) when the
    selection is invalid or ambiguous.
    """
    samples = get_vcf_sample_names(vcf_path)

    if sample is not None:
        if sample not in samples:
            raise VcfSampleSelectionError(
                f"Sample {sample!r} not found in VCF {vcf_path!r}. "
                f"Available samples: {samples}"
            )
        return sample

    if len(samples) == 1:
        return samples[0]

    raise VcfSampleSelectionError(
        f"VCF {vcf_path!r} has {len(samples)} samples ({samples}); "
        "use --sample to pick one."
    )


def read_biallelic_snvs(vcf_path: str, sample: str) -> Iterator[SnvSite]:
    """Yield biallelic SNV sites from a (optionally gzipped) VCF.

    A record is kept iff its FILTER is exactly PASS, it has a single-base
    REF, exactly one single-base ALT allele, and both alleles are plain
    A/C/G/T (excludes symbolic alleles such as ``*``). Indels, MNPs, and
    multiallelic records are skipped. Genotype/phasing is not used to
    select sites - haplotype assignment for the counting columns comes
    from the BAM's HP tag, not from VCF phase. The requested sample's GT
    (phased/unphased, "." for any missing allele) is carried through as
    SnvSite.genotype for the genotyping_from_input_vcf output column.
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

            sample_data = record.samples[sample]
            gt = sample_data.get("GT")
            genotype = format_genotype(gt, sample_data.phased) if gt else "."

            yield SnvSite(
                chrom=record.chrom, pos=record.pos, ref=ref, alt=alt, genotype=genotype
            )
