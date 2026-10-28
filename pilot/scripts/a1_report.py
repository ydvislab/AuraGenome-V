"""Collect the A1 (proxy-track selection) results into small tracked tables and figures under results/a1/,
so collaborators can read them without re-running the API scoring.

Inputs (untracked, produced by proxy_select.py): <pilot_dir>/results/proxy_select[_noise|_ensemble|_full|_test]/
Outputs (tracked): results/a1/*.tsv, results/a1/*.png, results/a1/test_decision.txt (when the test set exists)

usage: python a1_report.py <pilot_dir> <out_dir>
"""
import os, sys, glob, subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

pilot_dir, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
res = lambda tag="": os.path.join(pilot_dir, "results", "proxy_select" + (f"_{tag}" if tag else ""))

# reference palette (dataviz skill, light mode): fixed categorical order, method first
COLOR = {"disc": "#2a78d6", "random": "#eb6834", "high_effect": "#1baf7a", "disagree": "#eda100"}
LABEL = {"disc": "disc (sequential discrimination)", "random": "random", "high_effect": "high predicted effect",
         "disagree": "static disagreement"}
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "axes.spines.top": False, "axes.spines.right": False})


def primary(s):
    """hide_matched where the element has matched tracks, else all"""
    scen = s.groupby("element").scenario.agg(lambda x: "hide_matched" if (x == "hide_matched").any() else "all")
    return s[s.scenario == s.element.map(scen)]


def learning_curve(s, title, path):
    p = primary(s)
    curve = p[p.k > 0].pivot_table(index="strategy", columns="k", values="rho_bma", aggfunc="mean")
    ref = p[p.k == 0].groupby("strategy").rho_pick.mean()
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for name, y in [("oracle", ref["oracle"]), ("all_mean", ref["all_mean"])]:
        ax.axhline(y, color=INK2, lw=1, ls="--")
        ax.text(16.4, y, {"oracle": "best single track\n(uses all labels)", "all_mean": "all-track mean\n(no labels)"}[name],
                va="center", fontsize=8, color=INK2)
    for st in ["random", "high_effect", "disagree", "disc"]:
        y = curve.loc[st]
        ax.plot(y.index, y.values, color=COLOR[st], lw=2.6 if st == "disc" else 2, marker="o", ms=5,
                label=LABEL[st], zorder=3 if st == "disc" else 2)
    ax.set_xticks(curve.columns)
    ax.set_xlim(1.5, 16.3)
    ax.set_xlabel("measured variants k")
    ax.set_ylabel("Spearman with measured effect\n(held-out variants)")
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return curve.round(4), ref.round(4)


def per_unit(table, title, path):
    """dot plot: one row per element/region, strategies at k = 8 plus label-free references"""
    table = table.sort_values("oracle")
    fig, ax = plt.subplots(figsize=(6.4, 0.42 * len(table) + 1.4))
    yy = np.arange(len(table))
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    span = table[["all_mean", "oracle", "disc", "random", "high_effect"]]
    ax.hlines(yy, span.min(axis=1), span.max(axis=1), color=GRID, lw=1.5, zorder=1)
    ax.scatter(table.oracle, yy, marker="|", s=160, color=INK2, lw=2, label="best single track (all labels)", zorder=2)
    ax.scatter(table.all_mean, yy, marker="o", s=60, facecolor="none", edgecolor=INK2, lw=1.5,
               label="all-track mean (no labels)", zorder=2)
    for st in ["random", "high_effect", "disc"]:
        ax.scatter(table[st], yy, s=55 if st == "disc" else 40, color=COLOR[st], edgecolor=SURFACE, lw=1.5,
                   label=LABEL[st] + ", k = 8", zorder=3)
    ax.set_yticks(yy, table.index)
    ax.set_xlabel("Spearman with measured effect (held-out variants)")
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def k8_table(s):
    p = primary(s)
    row = lambda st, k, col: p[(p.strategy == st) & (p.k == k)].set_index("element")[col]
    return pd.DataFrame({"disc": row("disc", 8, "rho_bma"), "random": row("random", 8, "rho_bma"),
                         "high_effect": row("high_effect", 8, "rho_bma"), "disagree": row("disagree", 8, "rho_bma"),
                         "all_mean": row("all_mean", 0, "rho_pick"), "oracle": row("oracle", 0, "rho_pick"),
                         "scenario": p.groupby("element").scenario.first()})


# development set
dev = pd.read_csv(os.path.join(res(), "summary.tsv"), sep="\t")
dev.to_csv(os.path.join(out, "dev_summary.tsv"), sep="\t", index=False, float_format="%.4f")
curve, ref = learning_curve(dev, "Development set (5 elements), own cell type's tracks hidden",
                            os.path.join(out, "fig1_dev_learning_curve.png"))
curve.assign(**{f"ref_{k}": v for k, v in ref.items()}).to_csv(os.path.join(out, "dev_learning_curve.tsv"), sep="\t")
k8 = k8_table(dev)
k8.round(4).to_csv(os.path.join(out, "dev_k8_by_element.tsv"), sep="\t")
per_unit(k8, "Development set, k = 8 (posterior-weighted prediction)", os.path.join(out, "fig2_dev_by_element_k8.png"))

# ablations (mean over development elements, primary scenario, posterior-weighted prediction)
rows = []
for tag, name in [("", "v1 (frozen)"), ("noise", "+ measurement noise"), ("ensemble", "+ ensemble hypotheses"),
                  ("full", "+ both")]:
    f = os.path.join(res(tag), "summary.tsv")
    if not os.path.exists(f):
        continue
    p = primary(pd.read_csv(f, sep="\t"))
    for st in ["disc", "random", "high_effect"]:
        for k in [2, 8, 16]:
            q = p[(p.strategy == st) & (p.k == k)]
            rows.append(dict(version=name, strategy=st, k=k, rho_bma=q.rho_bma.mean(), rho_pick=q.rho_pick.mean()))
pd.DataFrame(rows).round(4).to_csv(os.path.join(out, "dev_ablation.tsv"), sep="\t", index=False)

# per-track Spearman: top 10 and the metadata-matched tracks per development element
tops = []
for f in sorted(os.listdir(os.path.join(res(), "tracks"))):
    t = pd.read_csv(os.path.join(res(), "tracks", f), sep="\t").assign(element=f[:-4])
    t["rank"] = np.arange(1, len(t) + 1)
    tops.append(t[(t["rank"] <= 10) | t.matched])
pd.concat(tops)[["element", "rank", "track", "output_type", "biosample_name", "biosample_type", "matched", "spearman",
                 "dir_acc_sig", "cluster_size"]].round(4).to_csv(os.path.join(out, "dev_top_tracks.tsv"), sep="\t", index=False)

# held-out test set (pre-registered), when it exists
if os.path.exists(os.path.join(res("test"), "summary.tsv")):
    test = pd.read_csv(os.path.join(res("test"), "summary.tsv"), sep="\t")
    test.to_csv(os.path.join(out, "test_summary.tsv"), sep="\t", index=False, float_format="%.4f")
    learning_curve(test, "Held-out test set (17 experiments), own cell type's tracks hidden",
                   os.path.join(out, "fig3_test_learning_curve.png"))
    k8 = k8_table(test)
    k8.round(4).to_csv(os.path.join(out, "test_k8_by_experiment.tsv"), sep="\t")
    per_unit(k8, "Held-out test set, k = 8 (posterior-weighted prediction)", os.path.join(out, "fig4_test_by_experiment_k8.png"))
    decision = subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "prereg_test.py"), pilot_dir],
                              capture_output=True, text=True, check=True).stdout
    open(os.path.join(out, "test_decision.txt"), "w").write(decision)
# few-label regression baselines (added after the pre-registration; first 50 repeats), region level
reg_dirs = sorted(glob.glob(os.path.join(pilot_dir, "results", "proxy_select_reg*")))
if reg_dirs:
    r = pd.concat([pd.read_csv(os.path.join(d, "runs.tsv.gz"), sep="\t",
                               usecols=["element", "scenario", "strategy", "k", "repeat", "rho_bma"]) for d in reg_dirs])
    r = r[(r.repeat >= 0) & (r.repeat < 50)]
    r = r[r.scenario == r.element.map(r.groupby("element").scenario.agg(
        lambda x: "hide_matched" if (x == "hide_matched").any() else "all"))]
    var = pd.read_csv(os.path.join(pilot_dir, "data", "variants_grch38.tsv"), sep="\t", low_memory=False)
    reg = var.groupby("element").agg(c=("chrom", "first"), a=("pos", "min"), b=("pos", "max"))
    r["region"] = r.element.map(reg.c + ":" + reg.a.astype(str) + "-" + reg.b.astype(str))
    dev_elements = set(dev.element)
    r["set"] = np.where(r.element.isin(dev_elements), "dev", "test")
    strategies = ["random", "disc", "random+ridge", "random+lasso", "disc+ridge", "disc+lasso"]
    e = (r[r.strategy.isin(strategies)].groupby(["set", "region", "element", "strategy", "k"]).rho_bma.mean()
         .groupby(["set", "region", "strategy", "k"]).mean().reset_index())
    e.round(4).to_csv(os.path.join(out, "regression_baselines_by_region.tsv"), sep="\t", index=False)
    lines = []
    for st in ["dev", "test"]:
        t = e[(e.set == st) & (e.k == 8)].pivot(index="region", columns="strategy", values="rho_bma")
        for b in ["random+lasso", "random+ridge", "disc+lasso", "disc+ridge"]:
            d = t.disc - t[b]
            p = wilcoxon(t.disc, t[b], alternative="greater").pvalue
            lines.append(f"{st}\tk=8\tdisc (posterior-weighted) - {b}\t{d.mean():+.3f}\t{int((d > 0).sum())}/{len(d)}\t{p:.4f}")
    open(os.path.join(out, "regression_baselines_k8.tsv"), "w").write(
        "set\tbudget\tcomparison\tmean_diff\tregions_positive\twilcoxon_p_one_sided\n" + "\n".join(lines) + "\n")
print("wrote", sorted(os.listdir(out)))
