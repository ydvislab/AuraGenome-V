"""Directional accuracy on significant held-out variants (docs/internal/direction_accuracy_zh.md; post hoc,
descriptive). Reads the existing runs; no model is rerun.

Labelled methods are scored by the posterior-average prediction (acc_bma), label-free references by the sign of the
track prediction (acc_pick), on the same held-out variants. Repeats without a significant held-out variant are
skipped; experiments are averaged within a Kircher region and the two cell types within a Siraj locus.

usage: python direction_report.py <pilot_dir> <out_dir>      e.g. pilot results/a1/direction
Outputs: <out_dir>/direction_contrasts.tsv, direction_levels.tsv
Diagnostic added after the first results (set kircher_v2_no_offset): the same v2 average without the fitted offsets
(acc_bma0), from a rerun of the frozen v2 on the Kircher test set tagged diag_dir/<element>.
"""
import os, sys, glob
import numpy as np
import pandas as pd

pilot_dir, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
res = os.path.join(pilot_dir, "results")
TEST = ["FOXE1", "GP1BA", "HBB", "HBG1", "IRF4", "LDLR", "LDLR.2", "MSMB", "MYCrs11986220", "MYCrs6983267", "SORT1",
        "SORT1-flip", "SORT1.2", "TERT-GAa", "TERT-GBM", "TERT-GSc", "TERT-HEK"]
B = 10000


def load(paths, unit_of, col="acc_bma"):
    r = pd.concat([pd.read_csv(f, sep="\t", usecols=["element", "scenario", "strategy", "k", "repeat", "acc_pick",
                                                     col]) for f in paths])
    r["acc"] = r[col].where(r.k > 0, r.acc_pick)
    r["unit"] = r.element.map(unit_of)
    # repeats first (NaN when no significant held-out variant), then experiments within a unit
    return r.groupby(["unit", "element", "scenario", "strategy", "k"]).acc.mean().groupby(
        ["unit", "scenario", "strategy", "k"]).mean()


def boot_ci(d):
    m = np.random.default_rng(0).choice(d, (B, len(d)), replace=True).mean(1)
    return np.quantile(m, [0.025, 0.975])


var = pd.read_csv(os.path.join(pilot_dir, "data", "variants_grch38.tsv"), sep="\t", low_memory=False)
reg = var.groupby("element").agg(c=("chrom", "first"), a=("pos", "min"), b=("pos", "max"))
region_of = (reg.c + ":" + reg.a.astype(str) + "-" + reg.b.astype(str)).to_dict()
sets = {
    "siraj": (load(glob.glob(os.path.join(res, "proxy_select_siraj_ref", "*", "runs.tsv.gz")),
                   lambda e: e.rsplit(".", 1)[0]),
              [("P1", "all", ("disc", 8), ("matched", 0)), ("P1 (S4)", "all", ("disc", 16), ("matched", 0)),
               ("P3", "hide_matched", ("disc", 8), ("all_mean", 0)), ("P4", "hide_exact", ("disc", 16), ("related", 0)),
               ("S1", "hide_matched", ("disc", 8), ("random", 8)),
               ("S3", "all", ("random+matched_cal", 8), ("matched", 0)),
               ("S3b", "all", ("disc+matched_cal", 8), ("matched", 0))]),
    "kircher_v1": (load([os.path.join(res, f"proxy_select_ce_{e}", "runs.tsv.gz") for e in TEST], region_of),
                   [("H1", "hide_matched", ("disc", 8), ("random", 8)),
                    ("H2", "hide_matched", ("disc", 8), ("all_mean", 0))]),
    "kircher_v2": (load([os.path.join(res, "proxy_select_v2", "p0.95_fam", e, "runs.tsv.gz") for e in TEST], region_of),
                   [("v2 C1 k=8", "all", ("disc", 8), ("matched", 0)), ("v2 C1 k=16", "all", ("disc", 16), ("matched", 0))]),
}
diag = [os.path.join(res, "proxy_select_diag_dir", e, "runs.tsv.gz") for e in TEST]
if all(os.path.exists(f) for f in diag):
    sets["kircher_v2_no_offset"] = (load(diag, region_of, "acc_bma0"),
                                    [("v2 C1 k=8", "all", ("disc", 8), ("matched", 0)),
                                     ("v2 C1 k=16", "all", ("disc", 16), ("matched", 0)),
                                     ("v2 C3 k=8", "hide_matched", ("disc", 8), ("all_mean", 0)),
                                     ("v2 C3 k=8 vs random", "hide_matched", ("disc", 8), ("random", 8))])
    sets["kircher_v2_with_offset_c3"] = (load(diag, region_of),
                                         [("v2 C3 k=8", "hide_matched", ("disc", 8), ("all_mean", 0)),
                                          ("v2 C3 k=8 vs random", "hide_matched", ("disc", 8), ("random", 8))])
rows, levels = [], []
for name, (acc, contrasts) in sets.items():
    for tag, scen, (s1, k1), (s2, k2) in contrasts:
        a1, a2 = acc.xs((scen, s1, k1), level=[1, 2, 3]), acc.xs((scen, s2, k2), level=[1, 2, 3])
        D = (a1 - a2).dropna()
        lo, hi = boot_ci(D.to_numpy())
        rows.append(dict(set=name, contrast=tag, view=scen, first=f"{s1}@{k1}", second=f"{s2}@{k2}", units=len(D),
                         acc_first=a1[D.index].mean(), acc_second=a2[D.index].mean(), mean=D.mean(), ci_low=lo,
                         ci_high=hi, positive=int((D > 0).sum()), negative=int((D < 0).sum()),
                         worse_by_005=int((D < -0.05).sum())))
    lv = acc.groupby(level=[1, 2, 3]).agg(["mean", "count"]).reset_index()
    levels.append(lv[lv.k.isin([0, 8, 16])].assign(set=name))
rows = pd.DataFrame(rows)
rows.round(4).to_csv(os.path.join(out, "direction_contrasts.tsv"), sep="\t", index=False)
pd.concat(levels).round(4).to_csv(os.path.join(out, "direction_levels.tsv"), sep="\t", index=False)
with pd.option_context("display.width", 220):
    print(rows.round(3).to_string(index=False))
