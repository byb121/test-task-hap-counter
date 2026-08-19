"""Command-line entrypoint for hap-counter."""

import csv

import click
import pysam

from hap_counter.counting import COUNT_FIELDS, count_haplotype_support_batches
from hap_counter.vcf_io import read_biallelic_snvs


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
def main(bam: str, vcf: str, output: str) -> None:
    """Compute per-haplotype REF/ALT support for biallelic SNVs.

    For every biallelic SNV in VCF, counts how many primary alignments
    assigned to each haplotype (via the BAM's HP tag) support the REF vs.
    ALT allele, and writes the result as a TSV.
    """
    fieldnames = ["chrom", "pos", "ref", "alt", *COUNT_FIELDS]

    with pysam.AlignmentFile(bam, "rb") as bam_file, open(output, "w", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t")
        writer.writerow(fieldnames)

        for site, counts in count_haplotype_support_batches(bam_file, read_biallelic_snvs(vcf)):
            writer.writerow(
                [site.chrom, site.pos, site.ref, site.alt]
                + [counts[field] for field in COUNT_FIELDS]
            )


if __name__ == "__main__":
    main()
