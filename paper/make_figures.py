"""Figures for the manuscript that are not already in results/a1/ (run from the repository root).

usage: .venv/bin/python paper/make_figures.py
Inputs: results/a1/siraj_ref/siraj_by_locus.tsv (per-locus paired differences, pre-registered Siraj test, labels
        relative to the measured reference oligo; docs/internal/siraj_refnorm_zh.md)
Outputs: paper/figures/fig_siraj.pdf
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID, SURFACE, METHOD = "#0b0b0b", "#52514e", "#e4e3df", "#ffffff", "#2a78d6"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
d = pd.read_csv("results/a1/siraj_ref/siraj_by_locus.tsv", sep="\t", index_col=0)
panels = [("C1_D_k8", "C1 (matched kept), $k=8$\nv2 $-$ matched mean"),
          ("C2_D_k16", "C2 (own cell line hidden), $k=16$\nv2 $-$ related-cell mean"),
          ("C3_D_k8", "C3 (own and related hidden), $k=8$\nv2 $-$ all-track mean")]
fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.3), sharey=True)
bins = np.arange(-0.15, 0.45001, 0.015)   # covers every locus (range -0.128 to +0.439)
for ax, (col, title) in zip(axes, panels):
    x = d[col].dropna()
    ax.hist(x, bins=bins, color=METHOD, edgecolor=SURFACE, linewidth=0.6)
    ax.axvline(0, color=INK2, linewidth=0.8)
    if col == "C1_D_k8":
        ax.axvline(-0.05, color=INK2, linewidth=0.8, linestyle=":")
        ax.text(-0.052, ax.get_ylim()[1] * 0.92, "harm\nthreshold", ha="right", va="top", color=INK2, fontsize=6.5)
    ax.set_title(title, fontsize=8, color=INK)
    ax.set_xlabel("difference in Spearman correlation", fontsize=7)
    ax.text(0.97, 0.95, f"mean {x.mean():+.3f}\n{(x > 0).sum()}/{len(x)} > 0", transform=ax.transAxes,
            ha="right", va="top", fontsize=7, color=INK)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
axes[0].set_ylabel("loci")
fig.tight_layout()
fig.savefig("paper/figures/fig_siraj.pdf")
print("wrote paper/figures/fig_siraj.pdf")
