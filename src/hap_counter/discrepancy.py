"""Aggregate discrepant-vote statistics for confidently-genotyped haplotypes.

A "discrepant vote" is a read that supported the allele *not* called for its
haplotype at a site, among sites where genotyping_from_input_bam actually made
a confident call (i.e. "." no-call sites are excluded).
"""

from collections import defaultdict
from typing import Dict, List, NamedTuple, Tuple


class GroupSummary(NamedTuple):
    chrom: str
    hap: int
    n_sites: int
    discrepant_votes: int
    total_votes: int
    overall_rate: float


class DiscrepancyAccumulator:
    def __init__(self) -> None:
        self._fractions: Dict[Tuple[str, int], List[float]] = defaultdict(list)
        self._discrepant_votes: Dict[Tuple[str, int], int] = defaultdict(int)
        self._total_votes: Dict[Tuple[str, int], int] = defaultdict(int)

    def add(self, chrom: str, hap: int, ref_count: int, alt_count: int, call: str) -> None:
        """Record one haplotype's vote at one site, if it got a confident call."""
        if call == ".":
            return

        total = ref_count + alt_count
        discrepant = min(ref_count, alt_count)
        key = (chrom, hap)

        self._fractions[key].append(discrepant / total)
        self._discrepant_votes[key] += discrepant
        self._total_votes[key] += total

    def fractions_by_group(self) -> Dict[Tuple[str, int], List[float]]:
        return dict(self._fractions)

    def summary(self) -> List[GroupSummary]:
        summaries = []
        for key, fractions in self._fractions.items():
            chrom, hap = key
            discrepant = self._discrepant_votes[key]
            total = self._total_votes[key]
            summaries.append(
                GroupSummary(
                    chrom=chrom,
                    hap=hap,
                    n_sites=len(fractions),
                    discrepant_votes=discrepant,
                    total_votes=total,
                    overall_rate=discrepant / total,
                )
            )
        return summaries

    def format_summary_text(self) -> str:
        summaries = self.summary()
        if not summaries:
            return "No confidently-genotyped haplotype calls to compute discrepancy stats from."

        lines = ["Discrepant-vote summary (per chromosome/haplotype):"]
        for group in sorted(summaries, key=lambda g: (g.chrom, g.hap)):
            lines.append(
                f"  {group.chrom} h{group.hap}: {group.discrepant_votes} discrepant votes "
                f"/ {group.total_votes} total votes across {group.n_sites} sites "
                f"(rate {group.overall_rate:.3f})"
            )
        return "\n".join(lines)
