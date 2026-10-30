"""Post hoc acquisition comparators (docs/internal/acquisition_comparators_zh.md): one-round panels, expected
information gain and the posterior-average stopping rule, against the frozen v2 sequential panel (disc).

Inputs (untracked): <root>/<unit>/runs.tsv.gz from proxy_select.py ... --comparators (run_comparators.sh).
Only repeats < --comparator-repeats carry the comparators, so every contrast is paired on those repeats.
Units: Kircher, the genomic region (experiments sharing a region averaged first); Siraj, the locus (its two cell
types averaged first).

usage: python cmp_report.py kircher <pilot_dir> <root> <out_dir> [dev|test]
       python cmp_report.py siraj   <pilot_dir> <root> <out_dir>
Outputs: <out_dir>/cmp_<set>_contrasts.tsv, cmp_<set>_levels.tsv, cmp_<set>_stop.tsv, and for the
development set the frozen stopping threshold in cmp_dev_stop_choice.txt
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

data, pilot_dir, root, out = sys.argv[1:5]
which = sys.argv[5] if len(sys.argv) > 5 else "siraj"
os.makedirs(out, exist_ok=True)
DEV = ["PKLR-48h", "F9", "HNF4A", "IRF6", "BCL11A"]
B, DELTA, TOL = 10000, 0.02, 0.01

runs = pd.concat([pd.read_csv(f, sep="\t", low_memory=False)
                  for f in sorted(glob.glob(os.path.join(root, "*", "runs.tsv.gz")))])
runs = runs[runs.repeat < runs[runs.strategy == "eig"].repeat.max() + 1]   # the repeats carrying comparators
if data == "kircher":
    var = pd.read_csv(os.path.join(pilot_dir, "data", "variants_grch38.tsv"), sep="\t", low_memory=False)
    reg = var.groupby("element").agg(c=("chrom", "first"), a=("pos", "min"), b=("pos", "max"))
    unit_of = reg.c + ":" + reg.a.astype(str) + "-" + reg.b.astype(str)
    keep = runs.element.isin(DEV) if which == "dev" else ~runs.element.isin(DEV)
    runs = runs[keep]
    runs["unit"] = runs.element.map(unit_of)
    # no anchor: hide_matched where the element has matched tracks, else all (as in the pre-registered test)
    has_hidden = runs.groupby("element").scenario.agg(lambda x: (x == "hide_matched").any())
    noanchor = np.where(runs.element.map(has_hidden), "hide_matched", "all")
    runs["view"] = np.where(runs.scenario == noanchor, "no_anchor", "")
    deploy = runs[(runs.scenario == "all") & runs.element.map(has_hidden)].assign(view="deploy")
    runs = pd.concat([runs[runs.view != ""], deploy])
    tag = which
else:
    runs["unit"] = runs.element.str.rsplit(".", n=1).str[0]
    runs["view"] = runs.scenario.map({"all": "C1", "hide_exact": "C2", "hide_matched": "C3"})
    tag = "siraj"


def per_unit(strategy, k):
    q = runs[(runs.strategy == strategy) & (runs.k == k)]
    return q.groupby(["view", "element", "unit"]).agg(rho=("rho_bma", "mean"), labels=("labels", "mean")) \
            .groupby(["view", "unit"]).mean()


def boot(d):
    rng = np.random.default_rng(0)
    m = d[rng.integers(0, len(d), size=(B, len(d)))].mean(1)
    return np.quantile(m, 0.025), np.quantile(m, 0.975)


def classify(lo, hi):
    if lo > -DELTA and hi < DELTA:
        return "equivalent within 0.02"
    if lo > 0:
        return "first better"
    if hi < 0:
        return "second better"
    return "undetermined"


rows, levels = [], []
for k in [8, 16]:
    lv = {s: per_unit(s, k) for s in ["disc", "batch_kb", "disc_static", "eig", "random", "high_effect", "disagree"]}
    for s, t in lv.items():
        for view, g in t.groupby(level="view"):
            levels.append(dict(view=view, k=k, strategy=s, units=len(g), rho_bma=g.rho.mean()))
    for first, second in [("disc", "batch_kb"), ("disc", "disc_static"), ("eig", "disc"), ("batch_kb", "disc_static"),
                          ("disc", "high_effect"), ("disc", "disagree")]:
        d_all = (lv[first].rho - lv[second].rho).dropna()
        for view, d in d_all.groupby(level="view"):
            d = d.to_numpy()
            lo, hi = boot(d)
            p = wilcoxon(d).pvalue if np.any(d != 0) and len(d) > 5 else np.nan
            rows.append(dict(view=view, k=k, contrast=f"{first} - {second}", units=len(d), mean=d.mean(), ci_low=lo,
                             ci_high=hi, positive=int((d > 0).sum()), wilcoxon_p=p, reading=classify(lo, hi)))
con = pd.DataFrame(rows)
con.to_csv(os.path.join(out, f"cmp_{tag}_contrasts.tsv"), sep="\t", index=False, float_format="%.4f")
pd.DataFrame(levels).to_csv(os.path.join(out, f"cmp_{tag}_levels.tsv"), sep="\t", index=False, float_format="%.4f")

# stopping rule: labels used and accuracy against disc with 16 and 8 labels, per threshold
d16, d8 = per_unit("disc", 16), per_unit("disc", 8)
stop = []
for s in sorted(x for x in runs.strategy.unique() if x.startswith("disc_stop")):
    t = per_unit(s, 16)
    for view, g in t.groupby(level="view"):
        lo, hi = boot((g.rho - d16.rho.loc[g.index]).to_numpy())
        stop.append(dict(view=view, rule=s, units=len(g), labels=g.labels.mean(), rho_bma=g.rho.mean(),
                         minus_disc16=(g.rho - d16.rho.loc[g.index]).mean(), ci_low=lo, ci_high=hi,
                         minus_disc8=(g.rho - d8.rho.loc[g.index]).mean()))
stop = pd.DataFrame(stop)
stop.to_csv(os.path.join(out, f"cmp_{tag}_stop.tsv"), sep="\t", index=False, float_format="%.4f")

with pd.option_context("display.width", 200, "display.max_rows", 200):
    print(con.round(4).to_string(index=False))
    print(stop.round(4).to_string(index=False))

if tag == "dev":
    # declared rule: the threshold with the fewest mean labels whose posterior-average accuracy is within 0.01 of
    # disc with 16 labels in every view (means over development regions); none -> the stopping rule fails
    g = stop[stop.rule.str.startswith("disc_stop_bma@")]
    ok = g.groupby("rule").filter(lambda x: (x.minus_disc16 >= -TOL).all())
    if len(ok):
        choice = ok.groupby("rule").labels.mean().idxmin()
        text = f"FROZEN STOPPING RULE: {choice} (mean labels {ok[ok.rule == choice].labels.mean():.2f})"
    else:
        text = "STOPPING RULE FAILS on the development set: no threshold within 0.01 of disc at k = 16 in every view"
    open(os.path.join(out, "cmp_dev_stop_choice.txt"), "w").write(text + "\n")
    print(text)
