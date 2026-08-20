"""Command-line entrypoint for hap-counter."""

import csv
import warnings
from typing import List, Optional

import click
import pysam

from hap_counter.bam_io import get_bam_sample_names
from hap_counter.counting import (
    COUNT_FIELDS,
    DEFAULT_GENOTYPE_ALLELE_FRACTION,
    DEFAULT_MIN_HAPLOTAGGED_READS,
    call_bam_genotype,
    count_haplotype_support_batches,
)
from hap_counter.discrepancy import DiscrepancyAccumulator
from hap_counter.vcf_io import VcfSampleSelectionError, read_biallelic_snvs, resolve_vcf_sample


def _sample_mismatch_warning(bam_sample_names: List[str], vcf_sample: str) -> Optional[str]:
    """Return a warning message if the BAM's sample name(s) don't line up with
    the VCF sample being used, else None.

    A BAM with no sample name at all gets a distinct message from one whose
    sample name(s) actively disagree with the VCF sample.
    """
    if not bam_sample_names:
        return (
            f"BAM header has no sample name (@RG/SM); cannot compare against "
            f"VCF sample {vcf_sample!r}."
        )
    if vcf_sample not in bam_sample_names:
        return (
            f"BAM sample name(s) {bam_sample_names} do not match VCF sample {vcf_sample!r}."
        )
    return None


@click.command()
@click.option(
    "--bam",
    "-b",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Coordinate-sorted, indexed BAM with primary alignments optionally haplotagged (HP tag).",
)
@click.option(
    "--vcf",
    "-v",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Phased VCF (optionally bgzipped) with variant calls.",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False),
    help="Path to write the output TSV to.",
)
@click.option(
    "--sample",
    "-s",
    default=None,
    help="VCF sample to use for genotyping_from_input_vcf. Required if the VCF has "
    "more than one sample; defaults to the only sample otherwise.",
)
@click.option(
    "--min-genotyping-hap-reads",
    "-m",
    type=click.IntRange(min=0),
    default=DEFAULT_MIN_HAPLOTAGGED_READS,
    show_default=True,
    help="Minimum REF+ALT reads a haplotype needs at a site to attempt a "
    "genotyping_from_input_bam call.",
)
@click.option(
    "--genotype-allele-fraction",
    "-f",
    type=click.FloatRange(min=0, max=1, min_open=True),
    default=DEFAULT_GENOTYPE_ALLELE_FRACTION,
    show_default=True,
    help="Fraction of a haplotype's REF+ALT reads that must support one allele "
    "for genotyping_from_input_bam to call it.",
)
@click.option(
    "--discrepancy-plot",
    "-d",
    type=click.Path(dir_okay=False),
    default=None,
    help="If set, also compute per-(chromosome,haplotype) discrepant-vote "
    "statistics across confidently-genotyped SNVs, print the aggregate summary, "
    "and save a histogram of discrepancy fractions to this path.",
)
def main(
    bam: str,
    vcf: str,
    output: str,
    sample: Optional[str],
    min_genotyping_hap_reads: int,
    genotype_allele_fraction: float,
    discrepancy_plot: Optional[str],
) -> None:
    """Compute per-haplotype REF/ALT support for biallelic SNVs.

    For every biallelic SNV in VCF, counts how many primary alignments
    assigned to each haplotype (via the BAM's HP tag) support the REF vs.
    ALT allele, and writes the result as a TSV.
    """
    try:
        resolved_sample = resolve_vcf_sample(vcf, sample)
    except VcfSampleSelectionError as err:
        raise click.ClickException(str(err))

    with pysam.AlignmentFile(bam, "rb") as bam_header_check:
        bam_sample_names = get_bam_sample_names(bam_header_check.header.to_dict())
    mismatch_warning = _sample_mismatch_warning(bam_sample_names, resolved_sample)
    if mismatch_warning is not None:
        warnings.warn(mismatch_warning)

    fieldnames = [
        "chrom",
        "pos",
        "ref",
        "alt",
        *COUNT_FIELDS,
        "genotyping_from_input_vcf",
        "genotyping_from_input_bam",
    ]

    accumulator = DiscrepancyAccumulator() if discrepancy_plot else None

    with pysam.AlignmentFile(bam, "rb") as bam_file, open(output, "w", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t")
        writer.writerow(fieldnames)

        for site, counts in count_haplotype_support_batches(
            bam_file, read_biallelic_snvs(vcf, resolved_sample)
        ):
            bam_genotype = call_bam_genotype(
                counts, min_genotyping_hap_reads, genotype_allele_fraction
            )
            writer.writerow(
                [site.chrom, site.pos, site.ref, site.alt]
                + [counts[field] for field in COUNT_FIELDS]
                + [site.genotype, bam_genotype]
            )

            if accumulator is not None:
                h1_call, h2_call = bam_genotype.split("|")
                accumulator.add(site.chrom, 1, counts["h1_REF"], counts["h1_ALT"], h1_call)
                accumulator.add(site.chrom, 2, counts["h2_REF"], counts["h2_ALT"], h2_call)

    if accumulator is not None:
        from hap_counter.plotting import plot_discrepancy_histogram

        summary_text = accumulator.format_summary_text()
        click.echo(summary_text)
        plot_discrepancy_histogram(accumulator.fractions_by_group(), summary_text, discrepancy_plot)


if __name__ == "__main__":
    main()
