"""EIG convergence check (docs/internal/eig_convergence_zh.md; post hoc): the information-gain panel with the
published approximation (E0: 100 hypotheses, 96 quadrature points) against a wider truncation (E1: 300) and one that
keeps 99% of the posterior mass at every step with 192 points (E2), and each against disc.

Inputs (untracked): <root>/<unit>/runs.tsv.gz from proxy_select.py ... --eig-check (run_eig_check.sh). Views and
units as in cmp_report.py: Kircher test regions (deploy, no_anchor), Siraj loci (C1, C2, C3).

usage: python eig_check_report.py kircher <pilot_dir> <root> <out_dir>
       python eig_check_report.py siraj   <pilot_dir> <root> <out_dir>
Outputs: <out_dir>/eigconv_<set>_contrasts.tsv, eigconv_<set>_levels.tsv, eigconv_<set>_cover.tsv
"""
import os, sys, glob
import numpy as np
import pandas as pd

data, pilot_dir, root, out = sys.argv[1:5]
os.makedirs(out, exist_ok=True)
DEV = ["PKLR-48h", "F9", "HNF4A", "IRF6", "BCL11A"]
B, DELTA = 10000, 0.02
E = {"E0": "eig@100/96", "E1": "eig@300/96", "E2": "eig@0.99/192"}

runs = pd.concat([pd.read_csv(f, sep="\t", low_memory=False)
                  for f in sorted(glob.glob(os.path.join(root, "*", "runs.tsv.gz")))])
runs = runs[runs.repeat < runs[runs.strategy == E["E0"]].repeat.max() + 1]   # the repeats carrying the check
if data == "kircher":
    var = pd.read_csv(os.path.join(pilot_dir, "data", "variants_grch38.tsv"), sep="\t", low_memory=False)
    reg = var.groupby("element").agg(c=("chrom", "first"), a=("pos", "min"), b=("pos", "max"))
    runs = runs[~runs.element.isin(DEV)]
    runs["unit"] = runs.element.map(reg.c + ":" + reg.a.astype(str) + "-" + reg.b.astype(str))
    # as in cmp_report.py: no_anchor = hide_matched where the element has matched tracks, else all; deploy = all
    has_hidden = runs.groupby("element").scenario.agg(lambda x: (x == "hide_matched").any())
    noanchor = np.where(runs.element.map(has_hidden), "hide_matched", "all")
    runs["view"] = np.where(runs.scenario == noanchor, "no_anchor", "")
    deploy = runs[(runs.scenario == "all") & runs.element.map(has_hidden)].assign(view="deploy")
    runs = pd.concat([runs[runs.view != ""], deploy])
else:
    runs["unit"] = runs.element.str.rsplit(".", n=1).str[0]
    runs["view"] = runs.scenario.map({"all": "C1", "hide_exact": "C2", "hide_matched": "C3"})


def per_unit(strategy, k, col="rho_bma"):
    q = runs[(runs.strategy == strategy) & (runs.k == k)]
    return q.groupby(["view", "element", "unit"])[col].mean().groupby(["view", "unit"]).mean()


def boot(d):
    m = d[np.random.default_rng(0).integers(0, len(d), size=(B, len(d)))].mean(1)
    return np.quantile(m, [0.025, 0.975])


def classify(lo, hi):   # the reading rule of cmp_report.py
    if lo > -DELTA and hi < DELTA:
        return "equivalent within 0.02"
    return "first better" if lo > 0 else "second better" if hi < 0 else "undetermined"


rows, levels, cover = [], [], []
for k in [8, 16]:
    lv = {"disc": per_unit("disc", k)} | {e: per_unit(s, k) for e, s in E.items()}
    for s, t in lv.items():
        for view, g in t.groupby(level="view"):
            levels.append(dict(view=view, k=k, strategy=s, units=len(g), rho_bma=g.mean()))
    for first, second in [("E2", "E0"), ("E2", "disc"), ("E0", "disc"), ("E1", "E0")]:
        for view, d in (lv[first] - lv[second]).dropna().groupby(level="view"):
            d = d.to_numpy()
            lo, hi = boot(d)
            rows.append(dict(view=view, k=k, contrast=f"{first} - {second}", units=len(d), mean=d.mean(), ci_low=lo,
                             ci_high=hi, positive=int((d > 0).sum()), negative=int((d < 0).sum()),
                             reading=classify(lo, hi)))
    for e, s in E.items():
        for col in ["eig_cover_first", "eig_cover_mean", "eig_cover_min"]:
            for view, g in per_unit(s, k, col).groupby(level="view"):
                cover.append(dict(view=view, k=k, setting=e, quantity=col, units=len(g), median=g.median(),
                                  minimum=g.min(), below_half=int((g < 0.5).sum())))
rows, levels, cover = pd.DataFrame(rows), pd.DataFrame(levels), pd.DataFrame(cover)
tag = "test" if data == "kircher" else "siraj"
rows.to_csv(os.path.join(out, f"eigconv_{tag}_contrasts.tsv"), sep="\t", index=False, float_format="%.4f")
levels.to_csv(os.path.join(out, f"eigconv_{tag}_levels.tsv"), sep="\t", index=False, float_format="%.4f")
cover.to_csv(os.path.join(out, f"eigconv_{tag}_cover.tsv"), sep="\t", index=False, float_format="%.4f")
with pd.option_context("display.width", 200, "display.max_rows", 200):
    print(levels.round(4).to_string(index=False))
    print(rows.round(4).to_string(index=False))
    print(cover[cover.k == 8].round(3).to_string(index=False))
