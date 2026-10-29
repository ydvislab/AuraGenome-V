"""C1 robustness (docs/internal/c1_robustness_zh.md): A1 with the held-out-fold checkpoint versus ALL_FOLDS.

Inputs (untracked): <pilot_dir>/results/proxy_select_c1/{v1,v2}_{all,fold}/<element>/summary.tsv from run_c1.sh.
Outputs: <out_dir>/c1_by_region.tsv (every quantity per region and checkpoint, k = 8 and 16),
         <out_dir>/c1_summary.tsv (mean under each checkpoint, change with 95% bootstrap interval),
         <out_dir>/c1_decision.txt (the declared reading on the test set)
Unit: the genomic region (experiments sharing a region are averaged first).

usage: python c1_report.py <pilot_dir> <out_dir>
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

pilot_dir, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
DEV = ["PKLR-48h", "F9", "HNF4A", "IRF6", "BCL11A"]
B, HARM = 10000, 0.05
root = os.path.join(pilot_dir, "results", "proxy_select_c1")
var = pd.read_csv(os.path.join(pilot_dir, "data", "variants_grch38.tsv"), sep="\t", low_memory=False)
reg = var.groupby("element").agg(c=("chrom", "first"), a=("pos", "min"), b=("pos", "max"))
region_of = reg.c + ":" + reg.a.astype(str) + "-" + reg.b.astype(str)


def read(run):
    s = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(root, run, "*", "summary.tsv")))])
    return s


def primary(s):
    """hide_matched where the element has matched tracks, else all (as in the pre-registered test)"""
    scen = s.groupby("element").scenario.agg(lambda x: "hide_matched" if (x == "hide_matched").any() else "all")
    return s[s.scenario == s.element.map(scen)]


def get(s, st, k, col):
    return s[(s.strategy == st) & (s.k == k)].set_index("element")[col]


rows = []
for side in ["all", "fold"]:
    v1, v2 = read(f"v1_{side}"), read(f"v2_{side}")
    p1, a2, h2 = primary(v1), v2[v2.scenario == "all"], v2[v2.scenario == "hide_matched"]
    a1 = v1[v1.scenario == "all"]
    for k in [8, 16]:
        q = {
            "ref matched": get(a1, "matched", 0, "rho_pick"),
            "ref all_mean": get(a1, "all_mean", 0, "rho_pick"),
            "ref oracle": get(a1, "oracle", 0, "rho_pick"),
            "v1 disc": get(p1, "disc", k, "rho_bma"),
            "v1 disc - random": get(p1, "disc", k, "rho_bma") - get(p1, "random", k, "rho_bma"),
            "v1 disc - all_mean": get(p1, "disc", k, "rho_bma") - get(p1, "all_mean", 0, "rho_pick"),
            "v1 disc - high_effect": get(p1, "disc", k, "rho_bma") - get(p1, "high_effect", k, "rho_bma"),
            "v2 all: disc - matched": get(a2, "disc", k, "rho_bma") - get(a2, "matched", 0, "rho_pick"),
            "v2 hide_matched: disc - all_mean": get(h2, "disc", k, "rho_bma") - get(h2, "all_mean", 0, "rho_pick"),
        }
        for name, x in q.items():
            x = x.dropna()
            x = x.groupby(x.index.map(region_of)).mean()
            for r, val in x.items():
                rows.append(dict(side=side, k=k, quantity=name, region=r, value=val))
long = pd.DataFrame(rows)
dev_regions = set(region_of[DEV])
long["set"] = np.where(long.region.isin(dev_regions), "dev", "test")
wide = long.pivot_table(index=["set", "k", "quantity", "region"], columns="side", values="value").reset_index()
wide["change"] = wide["fold"] - wide["all"]
wide.round(4).to_csv(os.path.join(out, "c1_by_region.tsv"), sep="\t", index=False)


def ci(d):
    rng = np.random.default_rng(0)
    m = rng.choice(d, (B, len(d)), replace=True).mean(1)
    return np.quantile(m, [0.025, 0.975])


summ = []
for (st, k, qn), g in pd.concat([wide, wide.assign(set="all")]).groupby(["set", "k", "quantity"]):
    g = g.dropna(subset=["all", "fold"])
    lo, hi = ci(g.change.to_numpy())
    summ.append(dict(set=st, k=k, quantity=qn, regions=len(g), all_folds=g["all"].mean(), held_out=g["fold"].mean(),
                     change=g.change.mean(), ci_low=lo, ci_high=hi, change_positive=int((g.change > 0).sum()),
                     harmed_all=int((g["all"] < -HARM).sum()) if "disc - matched" in qn else np.nan,
                     harmed_fold=int((g["fold"] < -HARM).sum()) if "disc - matched" in qn else np.nan))
summ = pd.DataFrame(summ)
summ.round(4).to_csv(os.path.join(out, "c1_summary.tsv"), sep="\t", index=False)

# declared reading: test set, k = 8, contrasts keep their sign; H1 and H2 one-sided Wilcoxon (descriptive)
lines, robust = [], True
t = summ[(summ.set == "test") & (summ.k == 8)].set_index("quantity")
for qn in ["v1 disc - random", "v1 disc - all_mean", "v1 disc - high_effect", "v2 all: disc - matched",
           "v2 hide_matched: disc - all_mean"]:
    a, f = t.loc[qn, "all_folds"], t.loc[qn, "held_out"]
    same = np.sign(a) == np.sign(f)
    robust &= bool(same)
    lines.append(f"{qn}: ALL_FOLDS {a:+.4f} -> held-out fold {f:+.4f} (change {f - a:+.4f}, 95% CI "
                 f"[{t.loc[qn, 'ci_low']:+.4f}, {t.loc[qn, 'ci_high']:+.4f}]) sign {'kept' if same else 'FLIPPED'}")
w = wide[(wide.set == "test") & (wide.k == 8)]
for qn, tag in [("v1 disc - random", "H1"), ("v1 disc - all_mean", "H2")]:
    for side in ["all", "fold"]:
        d = w[w.quantity == qn][side].dropna()
        p = wilcoxon(d, alternative="greater").pvalue
        if side == "fold":
            robust &= bool(p < 0.05)
        lines.append(f"{tag} ({qn}) under {'ALL_FOLDS' if side == 'all' else 'held-out fold'}: mean {d.mean():+.4f}, "
                     f"{int((d > 0).sum())}/{len(d)} positive, one-sided Wilcoxon p = {p:.4g}")
lines.append("READING: " + ("ROBUST - every test-set contrast keeps its sign and H1, H2 stay at p < 0.05 with the "
                            "held-out checkpoint." if robust else
                            "DEPENDS ON TRAINING EXPOSURE - a contrast flips sign or H1/H2 lose p < 0.05."))
text = "\n".join(lines) + "\n"
open(os.path.join(out, "c1_decision.txt"), "w").write(text)
print(text)
with pd.option_context("display.width", 220, "display.max_rows", 200):
    print(summ[summ.k == 8].round(3).to_string(index=False))
