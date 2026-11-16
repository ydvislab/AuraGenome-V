"""Figure 1a, the problem: at one element, the 1,018 AlphaGenome tracks disagree about the measured effects.
IRF4 enhancer in SK-MEL-28 (Kircher test set; the same element as Figure 4), all substitutions, ALL_FOLDS.
Histogram of each track's Spearman correlation with the measured effects; ticks mark the metadata-matched tracks.
Run from the repository root (package installed; needs the local data under pilot/):
    python paper/make_fig1.py      -> paper/figures/fig1_problem.pdf  (numbers printed for the caption)
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "tests")
from paper_data import load, rank_corr  # noqa: E402

E = "IRF4"
X, y, matched, related = load(E, "data/variants_grch38.tsv", "predictions/alphagenome_all")
yy = y.to_numpy()
rho = np.array([rank_corr(X.iloc[:, j].to_numpy(), yy) for j in range(X.shape[1])])
m = np.asarray(matched)
rho_all = rank_corr(X.mean(axis=1).to_numpy(), yy)
rho_mat = rank_corr(X.loc[:, m].mean(axis=1).to_numpy(), yy)
print(f"{E}: {len(yy)} variants, {len(rho)} tracks, {m.sum()} matched; per-track rho {rho.min():+.2f} to "
      f"{rho.max():+.2f}; negative {np.mean(rho < 0):.0%} ({(rho < 0).sum()}); matched mean {rho_mat:+.2f}; "
      f"all-track mean {rho_all:+.2f}; matched tracks {rho[m].min():+.2f} to {rho[m].max():+.2f}")

plt.rcParams.update({"font.size": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "font.family": "sans-serif"})
fig, ax = plt.subplots(figsize=(2.55, 1.6))
bins = np.arange(-0.5, 0.75, 0.025)
ax.hist(rho[rho < 0], bins=bins, color="#c0504d", alpha=0.85, edgecolor="white", linewidth=0.3)
ax.hist(rho[rho >= 0], bins=bins, color="#9a9a9a", edgecolor="white", linewidth=0.3)
top = ax.get_ylim()[1]
for r in rho[m]:
    ax.plot([r, r], [-0.10 * top, -0.02 * top], color="#1f5fa8", lw=0.9, solid_capstyle="butt")
ax.set_ylim(-0.12 * top, top * 1.12)
for x, lab, c, dx, ha in [(rho_all, f"all-track\nmean {rho_all:+.2f}", "#404040", -0.025, "right"),
                          (rho_mat, f"matched\n{rho_mat:+.2f}", "#1f5fa8", 0.015, "left")]:
    ax.axvline(x, color=c, lw=0.8, ls="--", ymax=0.9)
    ax.text(x + dx, top * 1.08, lab, color=c, ha=ha, va="top", fontsize=6.3, linespacing=1.0)
ax.text(-0.47, top * 0.62, f"wrong sign\n{np.mean(rho < 0):.0%} of tracks", color="#a33b38", fontsize=6.3,
        va="top", linespacing=1.0)
ax.axvline(0, color="black", lw=0.5)
ax.set_xlabel("Spearman correlation of a track with the assay", labelpad=1.5)
ax.set_ylabel("Tracks", labelpad=1.5)
ax.set_xlim(-0.5, 0.92)
ax.set_xticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.tick_params(length=2, pad=1.5)
fig.tight_layout(pad=0.2)
fig.text(0.005, 0.985, "(a)", fontsize=7.5, fontweight="bold", va="top", family="serif")
fig.savefig("paper/figures/fig1_problem.pdf")
