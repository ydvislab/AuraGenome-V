"""Check the trackpilot package against every frozen v2 run of the paper (proxy_select.py output).

For every element and scenario, on the given repeats, the package recomputes the pool, the sequential pilot panel
and the posterior-average prediction, and the label-free anchor / all-track mean, and compares the held-out Spearman
correlations with the rows of the frozen runs. Prints the largest absolute difference per dataset.

usage (from the repository root, with the package installed: pip install -e .):
  python pilot/scripts/trackpilot_reproduce.py siraj   [--repeats 0,1]
  python pilot/scripts/trackpilot_reproduce.py kircher [--repeats 0,1]
"""
import os, sys, glob, argparse, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "tests"))
from paper_data import ROOT, SIRAJ, load, rank_corr, scenarios  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("dataset", choices=["siraj", "kircher"])
ap.add_argument("--repeats", default="0")
a = ap.parse_args()
reps = [int(r) for r in a.repeats.split(",")]

if a.dataset == "siraj":
    jobs = [(f"{os.path.basename(os.path.dirname(f))}.{c}", f)
            for f in sorted(glob.glob(os.path.join(ROOT, SIRAJ["runs"].format(locus="*")))) for c in ["K562", "HEPG2"]]
    vpath, pdir = SIRAJ["variants"], SIRAJ["pred_dir"]
else:
    jobs = [(os.path.basename(os.path.dirname(f)), f)
            for f in sorted(glob.glob(os.path.join(ROOT, "results/proxy_select_v2/p0.95_fam/*/runs.tsv.gz")))]
    vpath, pdir = "data/variants_grch38.tsv", "predictions/alphagenome_all"

t0, diffs, n_cmp = time.time(), [], 0
for element, runs_path in jobs:
    runs = pd.read_csv(runs_path, sep="\t")
    runs = runs[runs.element == element]
    X, y, matched, related = load(element, vpath, pdir)
    n = len(y)
    for scen, tp in scenarios(X, matched, related).items():
        if scen not in set(runs.scenario):
            continue
        for rep in reps:
            rng = np.random.default_rng([0, rep])
            pool = rng.choice(n, int(0.5 * n), replace=False)
            out = np.ones(n, bool)
            out[pool] = False
            order = tp.suggest(16, y.index[pool], oracle=y)
            r = runs[(runs.scenario == scen) & (runs.repeat == rep)]
            disc = r[r.strategy == "disc"].set_index("k").rho_bma
            for k in disc.index:
                diffs.append(abs(rank_corr(tp.fit(y.loc[order[:k]]).predict()[out], y[out]) - disc.loc[k]))
            ref = r[r.k == 0].set_index("strategy").rho_pick
            name = {"all": "matched", "hide_exact": "related"}[scen] if tp.n_anchor_tracks else "all_mean"
            prior = tp.prior_prediction()
            diffs.append(abs(rank_corr(prior[out], y[out]) - ref.loc[name]))
            n_cmp += len(disc) + 1
    print(f"{element}: max |diff| so far {max(diffs):.2e} ({time.time() - t0:.0f}s)", flush=True)
print(f"{a.dataset}: {len(jobs)} elements, {n_cmp} numbers compared, max |difference| = {max(diffs):.3e}")
