"""Standalone script to dump raw pileup reads at one position, for eyeballing
against hap-counter's output.

Usage:
    python scripts/inspect_pileup.py --bam test_data/*.bam --chrom chr16 --pos 28043697 --ref C --alt T
"""

import click
import pysam


@click.command()
@click.option("--bam", "-b", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--chrom", "-c", required=True)
@click.option("--pos", "-p", required=True, type=int, help="1-based position.")
@click.option("--ref", default=None, help="REF allele, to flag matching reads.")
@click.option("--alt", default=None, help="ALT allele, to flag matching reads.")
def main(bam: str, chrom: str, pos: int, ref: str, alt: str) -> None:
    ref = ref.upper() if ref else None
    alt = alt.upper() if alt else None

    with pysam.AlignmentFile(bam, "rb") as bam_file:
        columns = bam_file.pileup(
            chrom,
            pos - 1,
            pos,
            truncate=True,
            ignore_overlaps=False,
        )
        column = next(iter(columns), None)

    if column is None:
        click.echo(f"No pileup column at {chrom}:{pos} (zero coverage).")
        return

    rows = []
    tally = {}
    for pileup_read in column.pileups:
        alignment = pileup_read.alignment
        name = alignment.query_name
        flag = alignment.flag
        is_supp = alignment.is_supplementary
        is_sec = alignment.is_secondary
        hp = alignment.get_tag("HP") if alignment.has_tag("HP") else None

        if pileup_read.is_del:
            base, qual = "DEL", None
        elif pileup_read.is_refskip:
            base, qual = "SKIP", None
        else:
            base = alignment.query_sequence[pileup_read.query_position].upper()
            qual = alignment.query_qualities[pileup_read.query_position]

        allele = None
        if base in ("DEL", "SKIP"):
            allele = None
        elif ref and base == ref:
            allele = "REF"
        elif alt and base == alt:
            allele = "ALT"
        elif ref or alt:
            allele = "other"

        rows.append(
            {
                "read_name": name,
                "flag": flag,
                "supplementary": is_supp,
                "secondary": is_sec,
                "HP": hp,
                "base": base,
                "base_qual": qual,
                "mapq": alignment.mapping_quality,
                "allele": allele,
            }
        )

        key = (hp, allele if (ref or alt) else base)
        tally[key] = tally.get(key, 0) + 1

    header = ["read_name", "flag", "supp", "sec", "HP", "base", "base_qual", "mapq", "allele"]
    click.echo("\t".join(header))
    for row in rows:
        click.echo(
            "\t".join(
                str(row[field])
                for field in ["read_name", "flag", "supplementary", "secondary", "HP", "base", "base_qual", "mapq", "allele"]
            )
        )

    click.echo("")
    click.echo(f"Total reads in pileup column (incl. supplementary): {len(rows)}")
    click.echo("Tally by (HP, allele/base):")
    for key in sorted(tally, key=lambda k: (k[0] is None, k[0] or 0, str(k[1]))):
        click.echo(f"  HP={key[0]}, {'allele' if (ref or alt) else 'base'}={key[1]}: {tally[key]}")


if __name__ == "__main__":
    main()
