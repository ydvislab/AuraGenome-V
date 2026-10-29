"""Collect the method-v2 runs (docs/internal/method_v2_zh.md) into small tracked tables under results/a1/.

Inputs (untracked): <pilot_dir>/results/proxy_select_v2/<config>/<element>/summary.tsv, from
  proxy_select.py --common-eval --tag v2/<config>/<element> [v2 flags]
and the v1 common-evaluation runs <pilot_dir>/results/proxy_select_ce_<element>/summary.tsv for reference.
Outputs: <out_dir>/v2_<set>_long.tsv (every config x region x strategy at k = 8 and 16),
         <out_dir>/v2_<set>_configs.tsv (deployment and control comparisons per config, k = 8 and 16),
         <out_dir>/v2_<set>_pi_choice.tsv (the declared pi rule on all elements and leaving one out; dev only)
Analysis unit: the genomic region (experiments sharing a region are averaged first).

usage: python v2_report.py <pilot_dir> <out_dir> <set: dev|test>
"""
import os, sys, glob
import numpy as np
import pandas as pd

pilot_dir, out, which = sys.argv[1], sys.argv[2], sys.argv[3]
DEV = ["PKLR-48h", "F9", "HNF4A", "IRF6", "BCL11A"]
root = os.path.join(pilot_dir, "results", "proxy_select_v2")
var = pd.read_csv(os.path.join(pilot_dir, "data", "variants_grch38.tsv"), sep="\t", low_memory=False)
reg = var.groupby("element").agg(c=("chrom", "first"), a=("pos", "min"), b=("pos", "max"))
region_of = reg.c + ":" + reg.a.astype(str) + "-" + reg.b.astype(str)
TIE = 0.005     # declared: pi values within this of the best count as tied, and the smallest tied pi wins
HARM = 0.05     # a region is harmed when v2 - matched < -HARM

rows = []
for d in sorted(glob.glob(os.path.join(root, "*", "*", "summary.tsv"))):
    config, element = d.split(os.sep)[-3:-1]
    if (element in DEV) != (which == "dev"):
        continue
    rows.append(pd.read_csv(d, sep="\t").assign(config=config))
for element in sorted({r.element.iloc[0] for r in rows}):
    rows.append(pd.read_csv(os.path.join(pilot_dir, "results", f"proxy_select_ce_{element}", "summary.tsv"),
                            sep="\t").assign(config="v1"))
s = pd.concat(rows)
s["region"] = s.element.map(region_of)
s = s[(s.k.isin([0, 8, 16])) | (s.strategy == "disc_stop")]
cols = ["rho_bma", "rho_pick", "ensemble_mass", "labels"]
long = (s.groupby(["config", "region", "scenario", "strategy", "k"])[cols].mean().reset_index())
long.round(4).to_csv(os.path.join(out, f"v2_{which}_long.tsv"), sep="\t", index=False)

a = long[long.scenario == "all"]
v = lambda cfg, st, k, col="rho_bma": a[(a.config == cfg) & (a.strategy == st) & (a.k == k)].set_index("region")[col]
matched = a[(a.strategy == "matched") & (a.config == "v1")].set_index("region").rho_pick
lines = []
for cfg in sorted(a.config.unique()):
    for k in [8, 16]:
        disc = v(cfg, "disc", k)
        if disc.empty:
            continue
        d = (disc - matched).dropna()
        row = dict(config=cfg, k=k, regions=len(d), disc=disc.mean(), matched=matched[d.index].mean(),
                   disc_minus_matched=d.mean(), positive=(d > 0).sum(), harmed=(d < -HARM).sum(), worst=d.min(),
                   disc_map_minus_matched=(v(cfg, "disc", k, "rho_pick") - matched).mean(),
                   random_minus_matched=(v(cfg, "random", k) - matched).mean(),
                   anchor_mass=v(cfg, "disc", k, "ensemble_mass").mean())
        for ctl in ["disc+meta_ridge", "random+meta_ridge", "disc+matched_cal", "random+matched_cal"]:
            c = v(cfg, ctl, k)
            if not c.empty:
                row[f"disc_minus_{ctl}"] = (disc - c).mean()
                row[f"{ctl}_minus_matched"] = (c - matched).mean()
        lines.append(row)
configs = pd.DataFrame(lines)
# hide_matched: v2 there is v1 plus the family prior (no anchor)
h = long[(long.scenario == "hide_matched") & (long.strategy.isin(["disc", "random"])) & (long.k.isin([8, 16]))]
hide = h.pivot_table(index=["config", "k"], columns="strategy", values="rho_bma", aggfunc="mean").reset_index()
hide.columns = ["config", "k", "hide_disc", "hide_random"]
configs = configs.merge(hide, on=["config", "k"], how="outer")
configs.round(4).to_csv(os.path.join(out, f"v2_{which}_configs.tsv"), sep="\t", index=False)

if which == "dev":
    # declared rule: among the family-prior configs with an anchor, maximise mean(disc - matched) at k = 8;
    # ties (within TIE) go to the smallest pi. Then the same rule leaving out each region.
    fam = sorted(c for c in a.config.unique() if c.endswith("_fam") and c != "p0_fam")
    pi = lambda c: float(c[1:].split("_")[0])
    D = pd.DataFrame({c: v(c, "disc", 8) - matched for c in fam}).dropna()

    def choose(rows):
        m = D.loc[rows].mean()
        return min((c for c in fam if m[c] >= m.max() - TIE), key=pi)
    choice = [dict(left_out="(none)", chosen=choose(D.index), left_out_D=np.nan)]
    for r in D.index:
        c = choose(D.index.drop(r))
        choice.append(dict(left_out=r, chosen=c, left_out_D=D.loc[r, c]))
    pd.DataFrame(choice).round(4).to_csv(os.path.join(out, f"v2_{which}_pi_choice.tsv"), sep="\t", index=False)
    print(pd.DataFrame(choice).round(3).to_string(index=False))

with pd.option_context("display.width", 250, "display.max_columns", 30):
    print(configs.round(3).to_string(index=False))
