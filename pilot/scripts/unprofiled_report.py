"""Kircher experiments in human cell lines that the models do not profile (post hoc subgroup; the flag
biosample_is_proxy was set in build_variant_table.py before any TrackPilot run): HEL 92.1.7 (K562 tracks as the proxy),
SK-MEL-28 (melanocyte), three glioblastoma lines (SK-N-SH and astrocyte). Frozen v2 runs, k = 8, 200 repeats.

usage (from the repository root): python pilot/scripts/unprofiled_report.py
Output: results/a1/unprofiled/subgroup.tsv (per test region and model) and summary.txt
"""
import os, glob
import numpy as np
import pandas as pd

R = "pilot/results"
OUT = "results/a1/unprofiled"
DEV = ["PKLR-48h", "F9", "HNF4A", "IRF6", "BCL11A"]
v = pd.read_csv("pilot/data/variants_grch38.tsv", sep="\t", low_memory=False)
proxy = v.groupby("element").agg(proxy=("biosample_is_proxy", "first"), cell=("mpra_cell_line", "first"))


def region(e):
    return "TERT (glioblastoma)" if e in ("TERT-GAa", "TERT-GBM", "TERT-GSc") else e


def val(g, sc, st, k, col):
    q = g[(g.scenario == sc) & (g.strategy == st) & (g.k == k)]
    return q[col].mean()


rows = []
for model, files in [("AlphaGenome", glob.glob(f"{R}/proxy_select_v2/p0.95_fam/*/summary.tsv")),
                     ("Enformer", [f"{R}/proxy_select_enformer_kircher/summary.tsv"])]:
    s = pd.concat([pd.read_csv(f, sep="\t") for f in files])
    for e, g in s.groupby("element"):
        if e in DEV or e not in proxy.index or not proxy.proxy[e]:
            continue
        rows.append(dict(model=model, region=region(e), cell=proxy.cell[e],
                         proxy_mean=val(g, "all", "matched", 0, "rho_pick"), tp_c1=val(g, "all", "disc", 8, "rho_bma"),
                         all_mean=val(g, "hide_matched", "all_mean", 0, "rho_pick"),
                         tp_c3=val(g, "hide_matched", "disc", 8, "rho_bma"),
                         random_c3=val(g, "hide_matched", "random", 8, "rho_bma")))
d = pd.DataFrame(rows).groupby(["model", "region"]).agg(
    {"cell": "first", "proxy_mean": "mean", "tp_c1": "mean", "all_mean": "mean", "tp_c3": "mean", "random_c3": "mean"})
os.makedirs(OUT, exist_ok=True)
d.round(4).to_csv(f"{OUT}/subgroup.tsv", sep="\t")
lines = []
for model, t in d.groupby(level=0):
    for a, b, what in [("tp_c1", "proxy_mean", "with proxy tracks: TrackPilot - proxy mean"),
                       ("tp_c3", "all_mean", "proxy hidden: TrackPilot - all-track mean"),
                       ("tp_c3", "random_c3", "proxy hidden: pilot - random panel")]:
        x = t[a] - t[b]
        lines.append(f"{model}, {len(x)} regions, {what} = {x.mean():+.3f} ({int((x > 0).sum())}/{len(x)} > 0)")
open(f"{OUT}/summary.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
print(d.round(3).to_string())
