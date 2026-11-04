"""Applicability check (docs/internal/applicability_check_zh.md): after the pilot panel is measured, a permutation test
on the mixture's log evidence decides whether any track carries the assay's signal; if not (P > p*), the label-free
default is returned instead of the TrackPilot prediction. Everything else is the frozen TrackPilot (package).

usage (from the repository root, package installed):
  python pilot/scripts/applicability_check.py run --stage dev      # development units of the three datasets
  python pilot/scripts/applicability_check.py choose               # p* chosen on the development units
  python pilot/scripts/applicability_check.py run --stage test     # test units (only after p* is committed)
  python pilot/scripts/applicability_check.py evaluate --pstar P   # A1, A2 and secondary results
Output: pilot/results/applicability/<stage>.tsv.gz, results/a1/applicability/*.
"""
import os, sys, glob, argparse
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import wilcoxon

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "tests"))
from paper_data import ROOT, SIRAJ, load, rank_corr, scenarios  # noqa: E402

K, B, REPS = 8, 200, range(20)
GRID = [0.05, 0.1, 0.2, 0.3, 0.5]
SCENS = ["all", "hide_matched"]
OUT = os.path.join(ROOT, "results", "applicability")
REPORT = "results/a1/applicability"
KIRCHER_DEV = ["PKLR-48h", "F9", "HNF4A", "IRF6", "BCL11A"]


def kircher_region(e):
    return e.split(".")[0].split("-flip")[0] if not e.startswith("TERT") else "TERT"


def units():
    """(dataset, element, unit, stage, variants path, prediction dir) for every element."""
    out = []
    for f in sorted(glob.glob(os.path.join(ROOT, "results/proxy_select_v2/p0.95_fam/*"))):
        e = os.path.basename(f)
        out.append(("kircher", e, kircher_region(e), "dev" if e in KIRCHER_DEV else "test",
                    "data/variants_grch38.tsv", "predictions/alphagenome_all"))
    loci = sorted(os.path.basename(os.path.dirname(f))
                  for f in glob.glob(os.path.join(ROOT, SIRAJ["runs"].format(locus="*"))))
    dev = set(np.random.default_rng(2611).choice(loci, 40, replace=False))
    for locus in loci:
        for cell in ["K562", "HEPG2"]:
            out.append(("siraj", f"{locus}.{cell}", locus, "dev" if locus in dev else "test", SIRAJ["variants"],
                        SIRAJ["pred_dir"]))
    v = pd.read_csv(os.path.join(ROOT, "data/kosicki/proxy_variants_alphagenome.tsv"), sep="\t",
                    usecols=["element", "pred_element"]).drop_duplicates()
    regions = sorted(v.pred_element.unique())
    dev = set(np.random.default_rng(2612).choice(regions, 50, replace=False))
    for e, r in v.itertuples(index=False):
        out.append(("neuron", e, r, "dev" if r in dev else "test", "data/kosicki/proxy_variants_alphagenome.tsv",
                    "predictions/kosicki_alphagenome"))
    return out


def check(tp, rows, ys, rep):
    """Permutation P value of the mixture's log evidence for the measured panel."""
    obs = logsumexp(tp._posterior(rows, ys)[0])
    rng = np.random.default_rng([2, rep])
    perm = np.array([logsumexp(tp._posterior(rows, rng.permutation(ys))[0]) for _ in range(B)])
    return (1 + int((perm >= obs).sum())) / (B + 1)


def job(u):
    dataset, element, unit, stage, vpath, pdir = u
    X, y, matched, related = load(element, vpath, pdir)
    n, rows = len(y), []
    for scen, tp in scenarios(X, matched, related).items():
        if scen not in SCENS:
            continue
        prior = tp.prior_prediction()
        for rep in REPS:
            pool = np.random.default_rng([0, rep]).choice(n, int(0.5 * n), replace=False)
            out = np.ones(n, bool)
            out[pool] = False
            order = tp.suggest(K, y.index[pool], oracle=y)
            lab = y.loc[order]
            p = check(tp, tp._rows(lab.index), lab.to_numpy(float), rep)
            rows.append(dict(dataset=dataset, element=element, unit=unit, scenario=scen, rep=rep, p=p,
                             rho_tp=rank_corr(tp.fit(lab).predict()[out], y[out]),
                             rho_default=rank_corr(prior[out], y[out])))
    return rows


def gated(d, pstar):
    return np.where(d.p > pstar, d.rho_default, d.rho_tp)


def unit_means(d, pstar):
    d = d.assign(gated=gated(d, pstar))
    e = d.groupby(["dataset", "scenario", "element", "unit"])[["rho_tp", "rho_default", "gated"]].mean()
    return e.groupby(["dataset", "scenario", "unit"]).mean().reset_index()


def lower_bound(x):
    x = np.asarray(x)
    return np.quantile(np.random.default_rng(0).choice(x, (10000, len(x)), replace=True).mean(1), 0.05)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "choose", "evaluate"])
    ap.add_argument("--stage", choices=["dev", "test"])
    ap.add_argument("--pstar", type=float)
    ap.add_argument("--jobs", type=int, default=10)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if a.cmd == "run":
        todo = [u for u in units() if u[3] == a.stage]
        print(f"{a.stage}: {len(todo)} elements", flush=True)
        with Pool(a.jobs) as pool:
            res = pool.map(job, todo, chunksize=1)
        pd.DataFrame([r for rr in res for r in rr]).to_csv(f"{OUT}/{a.stage}.tsv.gz", sep="\t", index=False)
        print(f"wrote {OUT}/{a.stage}.tsv.gz")
    elif a.cmd == "choose":
        d = pd.read_csv(f"{OUT}/dev.tsv.gz", sep="\t")
        rows = []
        for ps in GRID:
            u = unit_means(d, ps)
            g = u.assign(diff=u.gated - u.rho_default).groupby(["dataset", "scenario"])["diff"].mean()
            rows.append(dict(pstar=ps, objective=g.groupby("dataset").mean().mean(), **{f"{i}_{j}": v for (i, j), v in g.items()}))
        t = pd.DataFrame(rows)
        u = unit_means(d, 1.0)
        base = u.assign(diff=u.rho_tp - u.rho_default).groupby(["dataset", "scenario"])["diff"].mean()
        print("no check (TrackPilot - default):", base.round(4).to_dict())
        print(t.round(4).to_string(index=False))
        best = t.objective.max()
        print("p* =", t[t.objective >= best - 0.002].pstar.max())
        print("share judged 'no signal' at each p*:")
        print(pd.DataFrame({ps: d.groupby(["dataset", "scenario"]).p.apply(lambda x: (x > ps).mean()) for ps in GRID}).round(2))
    else:
        d = pd.read_csv(f"{OUT}/test.tsv.gz", sep="\t")
        u = unit_means(d, a.pstar)
        both = u.groupby(["dataset", "unit"])[["rho_tp", "rho_default", "gated"]].mean()
        lines = [f"applicability check, p* = {a.pstar}, test units, k = {K}, {len(REPS)} repeats"]
        n = both.loc["neuron"]
        d1, d2 = n.gated - n.rho_tp, n.gated - n.rho_default
        p1 = wilcoxon(d1, alternative="greater").pvalue
        a1 = d1.mean() > 0 and p1 < 0.05 and lower_bound(d2) > -0.02
        lines.append(f"A1 neurons ({len(n)} regions, C1 and C3 averaged): gated - TrackPilot = {d1.mean():+.4f}, "
                     f"{int((d1 > 0).sum())}/{len(d1)} > 0, P = {p1:.3g}; gated - default = {d2.mean():+.4f}, lower "
                     f"bound {lower_bound(d2):+.4f} (> -0.02?) -> {'HOLDS' if a1 else 'FAILS'}")
        s = u[u.dataset == "siraj"].set_index(["scenario", "unit"])
        c3 = s.loc["hide_matched"].gated - s.loc["hide_matched"].rho_default
        c1 = s.loc["all"].gated - s.loc["all"].rho_default
        sb = both.loc["siraj"]
        keep = sb.gated - sb.rho_tp
        p3 = wilcoxon(c3, alternative="greater").pvalue
        a2 = c3.mean() >= 0.03 and p3 < 0.05 and lower_bound(c1) > -0.02 and lower_bound(keep) > -0.01
        lines.append(f"A2 Siraj ({len(sb)} loci): C3 gated - all-track mean = {c3.mean():+.4f} (>= 0.03?), P = {p3:.3g}; "
                     f"C1 gated - matched mean = {c1.mean():+.4f}, lower bound {lower_bound(c1):+.4f} (> -0.02?); "
                     f"gated - TrackPilot = {keep.mean():+.4f}, lower bound {lower_bound(keep):+.4f} (> -0.01?) "
                     f"-> {'HOLDS' if a2 else 'FAILS'}")
        lines.append(f"DECISION: {'SUPPORTED' if a1 and a2 else 'NOT SUPPORTED'} (A1 and A2)")
        sec = u.groupby(["dataset", "scenario"])[["rho_tp", "rho_default", "gated"]].mean()
        sec["share_no_signal"] = d.groupby(["dataset", "scenario"]).p.apply(lambda x: (x > a.pstar).mean())
        os.makedirs(REPORT, exist_ok=True)
        text = "\n".join(lines) + "\n"
        open(f"{REPORT}/decision.txt", "w").write(text)
        sec.round(4).to_csv(f"{REPORT}/secondary.tsv", sep="\t")
        print(text)
        print(sec.round(3))
