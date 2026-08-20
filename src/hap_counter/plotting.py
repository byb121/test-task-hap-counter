"""Plotting for discrepancy-fraction distributions."""

import matplotlib

matplotlib.use("Agg")

from typing import Dict, List, Tuple  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import AutoMinorLocator, MultipleLocator  # noqa: E402

_X_MAX = 0.5


def plot_discrepancy_histogram(
    fractions_by_group: Dict[Tuple[str, int], List[float]],
    summary_text: str,
    output_path: str,
) -> None:
    """Save a histogram of discrepancy fractions, one overlaid series per
    (chrom, hap) group, with summary_text rendered below the axes.
    """
    fig, ax = plt.subplots()

    items = sorted(fractions_by_group.items())
    if items:
        labels = [f"{chrom} h{hap}" for (chrom, hap), _ in items]
        data = [fractions for _, fractions in items]
        # Bin only up to the data's actual max, not the fixed axis range: if
        # bins spanned the full (0, _X_MAX) axis, a bin edge could land
        # partway through that unused space, and a value sitting exactly on
        # the true max would spill into the following (otherwise-empty) bin
        # (numpy bins an edge value into the bin to its right, except for the
        # last bin overall, which is right-closed). Binning only the
        # data's own span makes its true max the last, right-closed edge.
        data_max = max((value for group in data for value in group), default=0.0)
        bin_range = (0, max(data_max, 1e-6))
        # Passing multiple datasets to a single hist() call dodges the bars
        # side-by-side within each bin, instead of overlaying/alpha-blending
        # them (which made overlapping groups hard to tell apart).
        ax.hist(data, bins=20, range=bin_range, label=labels)

    ax.set_xlim(0, _X_MAX)
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(axis="x", which="major", length=7)
    ax.tick_params(axis="x", which="minor", length=4)

    ax.set_xlabel("Discrepancy fraction")
    ax.set_ylabel("Number of SNVs")
    ax.set_title("Per-haplotype discrepancy fraction")
    if fractions_by_group:
        ax.legend()

    fig.text(0.01, 0.01, summary_text, ha="left", va="bottom", fontsize=6, family="monospace")
    fig.tight_layout(rect=(0, 0.15, 1, 1))
    fig.savefig(output_path)
    plt.close(fig)
