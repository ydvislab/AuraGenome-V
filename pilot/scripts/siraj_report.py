"""Pre-registered analysis of method v2 on Siraj SatMut (docs/internal/preregistration_v2_siraj_zh.md, section 5-6).

Inputs (untracked): <root>/<locus>/summary.tsv from
  proxy_select.py --variants data/siraj/proxy_variants.tsv --pred-dir predictions/siraj_all --common-eval
                  --anchor-prior 0.95 --family-prior --controls --elements <locus>.K562,<locus>.HEPG2 --tag siraj/<locus>
and, for the family-prior ablation (S5), <root>_uni/<locus>/summary.tsv (same without --family-prior, scenario all).
Scenarios: all = C1 (matched cell line kept), hide_exact = C2 (related-cell anchor), hide_matched = C3 (no anchor).
Unit: the locus; its two cell types are averaged first.

usage: python siraj_report.py <root> <out_dir>        e.g. pilot/results/proxy_select_siraj results/a1
Outputs: <out_dir>/siraj_decision.txt, siraj_by_locus.tsv, siraj_secondary.tsv
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, beta

root, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
DELTA, TAU, Q, B = 0.02, 0.05, 0.15, 10000


def read(r):
    s = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(r, "*", "summary.tsv")))])
    s["locus"] = s.element.str.rsplit(".", n=1).str[0]
    s["cell"] = s.element.str.rsplit(".", n=1).str[1]
    return s


s = read(root)


def val(scen, st, k, col="rho_bma", data=None, by="locus"):
    d = s if data is None else data
    q = d[(d.scenario == scen) & (d.strategy == st) & (d.k == k)]
    return q.groupby(by)[col].mean()


def boot_lower(d, level=0.95):
    rng = np.random.default_rng(0)
    m = rng.choice(d, (B, len(d)), replace=True).mean(1)
    return np.quantile(m, 1 - level)


def boot_ci(d):
    rng = np.random.default_rng(0)
    m = rng.choice(d, (B, len(d)), replace=True).mean(1)
    return np.quantile(m, [0.025, 0.975])


lines, rows = [], {}
# P1, P2 (and S4 at k = 16): C1, v2 (disc, posterior-weighted) - unchanged matched mean
for k in [8, 16]:
    D = (val("all", "disc", k) - val("all", "matched", 0, "rho_pick")).dropna()
    rows[f"C1_D_k{k}"] = D
    lb = boot_lower(D.to_numpy())
    x, n = int((D < -TAU).sum()), len(D)
    ub = beta.ppf(0.95, x + 1, n - x) if x < n else 1.0
    tag = "" if k == 8 else " (S4, secondary)"
    lines.append(f"P1{tag} k={k}: mean D = {D.mean():+.4f}, one-sided 95% bootstrap lower bound {lb:+.4f} "
                 f"(> -{DELTA}?) -> {'HOLDS' if lb > -DELTA else 'FAILS'}   [{n} loci]")
    lines.append(f"P2{tag} k={k}: harmed (D < -{TAU}) {x}/{n} = {x / n:.3f}, one-sided 95% Clopper-Pearson upper "
                 f"bound {ub:.3f} (<= {Q}?) -> {'HOLDS' if ub <= Q else 'FAILS'}")
p1 = boot_lower(rows["C1_D_k8"].to_numpy()) > -DELTA
x, n = int((rows["C1_D_k8"] < -TAU).sum()), len(rows["C1_D_k8"])
p2 = (beta.ppf(0.95, x + 1, n - x) if x < n else 1.0) <= Q

# P3: C3, k = 8, v2 - label-free all-track mean
D3 = (val("hide_matched", "disc", 8) - val("hide_matched", "all_mean", 0, "rho_pick")).dropna()
rows["C3_D_k8"] = D3
p = wilcoxon(D3, alternative="greater").pvalue
p3 = D3.mean() >= 0.03 and p < 0.05
lines.append(f"P3 C3 k=8: v2 - all_mean = {D3.mean():+.4f} (>= 0.03?), one-sided Wilcoxon p = {p:.4g} -> "
             f"{'HOLDS' if p3 else 'FAILS'}   [{len(D3)} loci]")

# P4: C2, k = 16, v2 - unchanged related-cell mean
D4 = (val("hide_exact", "disc", 16) - val("hide_exact", "related", 0, "rho_pick")).dropna()
rows["C2_D_k16"] = D4
p = wilcoxon(D4, alternative="greater").pvalue
p4 = D4.mean() > 0 and p < 0.05
lines.append(f"P4 C2 k=16: v2 - related mean = {D4.mean():+.4f} (> 0?), one-sided Wilcoxon p = {p:.4g} -> "
             f"{'HOLDS' if p4 else 'FAILS'}   [{len(D4)} loci]")

if p1 and p2 and (p3 or p4):
    verdict = "SAFETY NET SUPPORTED: P1 and P2 hold, and P3 or P4 holds."
elif p1 and p2:
    verdict = "SAFE BUT NO RESCUE: P1 and P2 hold; P3 and P4 fail. Claim only safe metadata-anchored inference (+ S1)."
else:
    verdict = "SAFETY NET DROPPED: P1 or P2 fails (reviewer stop rule)."
lines.append("DECISION: " + verdict)

# secondary (effect sizes and 95% bootstrap intervals only)
sec = []
def add(name, d, test="none"):
    d = d.dropna()
    lo, hi = boot_ci(d.to_numpy())
    p = wilcoxon(d, alternative="two-sided").pvalue if test == "two-sided" else np.nan
    sec.append(dict(comparison=name, loci=len(d), mean=d.mean(), ci_low=lo, ci_high=hi, positive=(d > 0).sum(),
                    wilcoxon_p_two_sided=p))
for scen, c in [("all", "C1"), ("hide_matched", "C3")]:
    add(f"S1 {c} k=8: v2 disc panel - v2 random panel", val(scen, "disc", 8) - val(scen, "random", 8))
add("S2 C1 k=8: v2 - meta_ridge on the same disc panel", val("all", "disc", 8) - val("all", "disc+meta_ridge", 8),
    "two-sided")
add("S3 C1 k=8: matched_cal (random panel) - matched", val("all", "random+matched_cal", 8)
    - val("all", "matched", 0, "rho_pick"))
uni = os.path.normpath(root) + "_uni"
if glob.glob(os.path.join(uni, "*", "summary.tsv")):
    su = read(uni)
    add("S5 C1 k=8: family prior - uniform prior", val("all", "disc", 8) - val("all", "disc", 8, data=su), "two-sided")
for cell in sorted(s.cell.unique()):
    sc = s[s.cell == cell]
    add(f"S7 {cell} C1 k=8: v2 - matched",
        val("all", "disc", 8, data=sc) - val("all", "matched", 0, "rho_pick", data=sc))
    add(f"S7 {cell} C3 k=8: v2 - all_mean",
        val("hide_matched", "disc", 8, data=sc) - val("hide_matched", "all_mean", 0, "rho_pick", data=sc))
sec = pd.DataFrame(sec)

text = "\n".join(lines) + "\n"
open(os.path.join(out, "siraj_decision.txt"), "w").write(text)
pd.DataFrame(rows).round(4).to_csv(os.path.join(out, "siraj_by_locus.tsv"), sep="\t")
sec.round(4).to_csv(os.path.join(out, "siraj_secondary.tsv"), sep="\t", index=False)
print(text)
with pd.option_context("display.width", 200):
    print(sec.round(3).to_string(index=False))
