"""Cross-locus prior (docs/internal/cross_locus_prior_zh.md): carry the track posterior of other loci in the same cell
type into a new locus, as a mixture with the frozen prior. Everything else is the frozen TrackPilot (package).

For each cell type and scenario, every locus first gets the single-locus result (sequential pilot panel of 8 on its
pool, repeats 0..R-1). For a target locus i and repeat rep, m other loci are drawn (default_rng([1, rep, i]).permutation,
first m), their posteriors after 8 labels are averaged (q), and the target's prior becomes (1 - lam) prior + lam q.
The target then chooses its own panel and fits as usual. The target's own labels are only its k panel variants.

usage (from the repository root, package installed):
  python pilot/scripts/cross_locus.py run --stage dev            # lam grid, m = 10, on the 40 development loci
  python pilot/scripts/cross_locus.py choose                     # prints the lam chosen on the development loci
  python pilot/scripts/cross_locus.py run --stage test --lam L   # m in {1,3,10,30} on the 85 test loci
  python pilot/scripts/cross_locus.py evaluate                   # X1, X2 and secondary contrasts
Output: pilot/results/cross_locus/<stage>/<cell>_<scenario>.tsv.gz, results/a1/cross_locus/*.
"""
import os, sys, glob, argparse, time
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, beta

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "tests"))
from paper_data import ROOT, SIRAJ, load, rank_corr  # noqa: E402
from trackpilot import TrackPilot  # noqa: E402

LOCI = sorted(os.path.basename(os.path.dirname(f)) for f in glob.glob(os.path.join(ROOT, SIRAJ["runs"].format(locus="*"))))
DEV = set(np.random.default_rng(2611).choice(LOCI, 40, replace=False))
CELLS, SCENS = ["K562", "HEPG2"], ["all", "hide_exact", "hide_matched"]
OUT = os.path.join(ROOT, "results", "cross_locus")
REPORT = "results/a1/cross_locus"


def build(locus, cell, scen):
    X, y, matched, related = load(f"{locus}.{cell}", SIRAJ["variants"], SIRAJ["pred_dir"])
    if scen == "all":
        tp = TrackPilot(X, matched)
    elif scen == "hide_exact":
        tp = TrackPilot(X, related, matched)
    else:
        tp = TrackPilot(X, None, matched | related)
    return tp, y


def split(n, rep):
    pool = np.random.default_rng([0, rep]).choice(n, int(0.5 * n), replace=False)
    out = np.ones(n, bool)
    out[pool] = False
    return pool, out


def panel_result(tp, y, pool, out, ks):
    order = tp.suggest(max(ks), y.index[pool], oracle=y)
    rho, w = {}, None
    for k in ks:
        fit = tp.fit(y.loc[order[:k]])
        rho[k] = rank_corr(fit.predict()[out], y[out])
        if k == 8:
            w = pd.Series(fit.weights, index=tp.hypotheses)
    return rho, w


def job(args):
    cell, scen, stage, lams, ms, reps = args
    t0 = time.time()
    data = [build(l, cell, scen) for l in LOCI]
    rows, post = [], {}
    for i, (tp, y) in enumerate(data):            # single-locus results (and posteriors) for every locus
        for rep in range(reps):
            pool, out = split(len(y), rep)
            rho, post[i, rep] = panel_result(tp, y, pool, out, [4, 8])
            for k, r in rho.items():
                rows.append((LOCI[i], cell, scen, rep, "single", 0.0, 0, k, r))
            rows.append((LOCI[i], cell, scen, rep, "default", 0.0, 0, 0, rank_corr(tp.prior_prediction()[out], y[out])))
    targets = [i for i, l in enumerate(LOCI) if (l in DEV) == (stage == "dev")]
    for i in targets:
        tp, y = data[i]
        base = np.exp(tp.log_prior)
        others = [j for j in range(len(LOCI)) if j != i]
        for rep in range(reps):
            pool, out = split(len(y), rep)
            perm = np.random.default_rng([1, rep, i]).permutation(others)
            for m in ms:
                q = sum(post[j, rep].reindex(tp.hypotheses, fill_value=0.0).to_numpy() for j in perm[:m]) / m
                for lam in lams:
                    p = (1 - lam) * base + lam * q
                    tp.log_prior = np.log(p / p.sum())
                    rho, _ = panel_result(tp, y, pool, out, [4, 8])
                    for k, r in rho.items():
                        rows.append((LOCI[i], cell, scen, rep, "cross", lam, m, k, r))
        tp.log_prior = np.log(base)
    df = pd.DataFrame(rows, columns=["locus", "cell", "scenario", "rep", "method", "lam", "m", "k", "rho"])
    os.makedirs(os.path.join(OUT, stage), exist_ok=True)
    df.to_csv(os.path.join(OUT, stage, f"{cell}_{scen}.tsv.gz"), sep="\t", index=False, float_format="%.6g")
    return f"{cell} {scen} {stage}: {len(df)} rows, {time.time() - t0:.0f}s"


def read(stage):
    return pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(os.path.join(OUT, stage, "*.tsv.gz"))])


def by_locus(d, method, k, scen, lam=None, m=None):
    q = d[(d.method == method) & (d.k == k) & (d.scenario == scen)]
    if lam is not None:
        q = q[(q.lam == lam) & (q.m == m)]
    return q.groupby(["locus", "cell"]).rho.mean().groupby("locus").mean()   # repeats, then the two cells


def choose():
    d = read("dev")
    dev = d[d.locus.isin(DEV)]
    single = by_locus(dev, "single", 8, "hide_matched").mean()
    res = {lam: by_locus(dev, "cross", 8, "hide_matched", lam, 10).mean() for lam in sorted(dev[dev.method == "cross"].lam.unique())}
    best = max(res.values())
    lam = min(l for l, v in res.items() if v >= best - 0.002)
    for l, v in res.items():
        print(f"lam {l}: C3 k=8 mean rho {v:.4f} (single-locus {single:.4f})")
    print(f"chosen lam = {lam}")
    return lam


def evaluate():
    d = read("test")
    d = d[~d.locus.isin(DEV)]
    lam = d[d.method == "cross"].lam.iloc[0]
    os.makedirs(REPORT, exist_ok=True)
    rng = np.random.default_rng(0)
    lines = [f"cross-locus prior, lam = {lam}, m = 10, {d.locus.nunique()} test loci, repeats 0-{d.rep.max()}"]
    D1 = by_locus(d, "cross", 8, "hide_matched", lam, 10) - by_locus(d, "single", 8, "hide_matched")
    p1 = wilcoxon(D1, alternative="greater").pvalue
    x1 = D1.mean() >= 0.02 and p1 < 0.05
    lines.append(f"X1 C3 k=8: cross - single = {D1.mean():+.4f} (>= 0.02?), {int((D1 > 0).sum())}/{len(D1)} > 0, "
                 f"one-sided Wilcoxon p = {p1:.3g} -> {'HOLDS' if x1 else 'FAILS'}")
    D2 = by_locus(d, "cross", 8, "all", lam, 10) - by_locus(d, "default", 0, "all")
    lb = np.quantile(rng.choice(D2.to_numpy(), (10000, len(D2)), replace=True).mean(1), 0.05)
    x, n = int((D2 < -0.05).sum()), len(D2)
    ub = beta.ppf(0.95, x + 1, n - x) if x < n else 1.0
    x2 = lb > -0.02 and ub <= 0.15
    lines.append(f"X2 C1 k=8: cross - matched = {D2.mean():+.4f}, bootstrap lower bound {lb:+.4f} (> -0.02?), "
                 f"harmed {x}/{n}, upper bound {ub:.3f} (<= 0.15?) -> {'HOLDS' if x2 else 'FAILS'}")
    lines.append(f"DECISION: {'SUPPORTED' if x1 and x2 else 'NOT SUPPORTED'}")
    # secondary, descriptive
    sec = []
    for scen, ref in [("hide_matched", "all-track mean"), ("hide_exact", "related mean"), ("all", "matched mean")]:
        dflt = by_locus(d, "default", 0, scen)
        for k in [4, 8]:
            s_, c_ = by_locus(d, "single", k, scen), by_locus(d, "cross", k, scen, lam, 10)
            sec.append(dict(scenario=scen, k=k, default=dflt.mean(), single=s_.mean(), cross=c_.mean(),
                            cross_minus_single=(c_ - s_).mean(), positive=int((c_ - s_ > 0).sum()), loci=len(s_)))
        for m in sorted(d[d.method == "cross"].m.unique()):
            c_ = by_locus(d, "cross", 8, scen, lam, m)
            sec.append(dict(scenario=scen, k=8, m=int(m), cross=c_.mean(),
                            cross_minus_single=(c_ - by_locus(d, "single", 8, scen)).mean()))
    sec = pd.DataFrame(sec)
    text = "\n".join(lines) + "\n"
    open(f"{REPORT}/decision.txt", "w").write(text)
    sec.round(4).to_csv(f"{REPORT}/secondary.tsv", sep="\t", index=False)
    pd.DataFrame({"X1_C3_D_k8": D1, "X2_C1_D_k8": D2}).round(4).to_csv(f"{REPORT}/by_locus.tsv", sep="\t")
    print(text)
    print(sec.round(3).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "choose", "evaluate"])
    ap.add_argument("--stage", choices=["dev", "test"], default="dev")
    ap.add_argument("--lam", type=float, default=None)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    if a.cmd == "run":
        lams = [0.25, 0.5, 0.75, 0.9] if a.stage == "dev" else [a.lam]
        ms = [10] if a.stage == "dev" else [1, 3, 10, 30]
        assert a.stage == "dev" or a.lam is not None, "--lam is required for the test stage"
        with Pool(a.jobs) as pool:
            for msg in pool.imap_unordered(job, [(c, s, a.stage, lams, ms, a.reps) for c in CELLS for s in SCENS]):
                print(msg, flush=True)
    elif a.cmd == "choose":
        choose()
    else:
        evaluate()
