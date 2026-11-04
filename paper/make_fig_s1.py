"""Supplementary Fig. S1: Kircher test set, learning curve averaged over the 11 test regions (experiments in one
region averaged first), as in Table 2 of the paper. Input: pilot/results/proxy_select_ce_test/summary.tsv (common
evaluation set; the pre-registered version). Run from the repository root:
    python paper/make_fig_s1.py
Writes paper/figures/fig_s1_kircher_curve.pdf and results/a1/fig_s1_kircher_curve_by_region.tsv."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

s = pd.read_csv("pilot/results/proxy_select_ce_test/summary.tsv", sep="\t")
scen = s.groupby("element").scenario.agg(lambda x: "hide_matched" if (x == "hide_matched").any() else "all")
s = s[s.scenario == s.element.map(scen)]
s = s.assign(region=s.element.map(lambda e: "TERT" if e.startswith("TERT") else e.split(".")[0].split("-flip")[0]))
assert s.region.nunique() == 11
curve = (s[s.k > 0].groupby(["strategy", "k", "region"]).rho_bma.mean().groupby(["strategy", "k"]).mean()
         .unstack("k"))
ref = s[s.k == 0].groupby(["strategy", "region"]).rho_pick.mean().groupby("strategy").mean()
curve.assign(**{f"ref_{k}": v for k, v in ref.items()}).round(4).to_csv(
    "results/a1/fig_s1_kircher_curve_by_region.tsv", sep="\t")

STYLE = {"disc": ("pilot panel (TrackPilot)", "#d95f02", 2.6), "disagree": ("largest unweighted disagreement", "#7570b3", 1.8),
         "high_effect": ("largest predicted effect", "#1b9e77", 1.8), "random": ("random", "#666666", 1.8)}
fig, ax = plt.subplots(figsize=(6.4, 3.8))
ax.grid(axis="y", color="#e5e5e5", lw=0.8)
ax.set_axisbelow(True)
for name, label in [("oracle", "best single track\n(uses all labels)"), ("all_mean", "all-track mean\n(no labels)")]:
    ax.axhline(ref[name], color="#555555", lw=1, ls="--")
    ax.text(16.4, ref[name], label, va="center", fontsize=8, color="#555555")
for st, (label, color, lw) in STYLE.items():
    y = curve.loc[st].dropna()
    ax.plot(y.index, y.values, color=color, lw=lw, marker="o", ms=5, label=label, zorder=3 if st == "disc" else 2)
ax.set_xticks([2, 4, 8, 12, 16])
ax.set_xlim(1.5, 16.3)
ax.set_xlabel("measured variants k")
ax.set_ylabel("Spearman correlation\n(held-out variants)")
ax.legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig("paper/figures/fig_s1_kircher_curve.pdf")
print(curve.round(3).to_string())
print(ref.round(3).to_string())
