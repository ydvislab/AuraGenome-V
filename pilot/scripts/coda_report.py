"""CODA / MODEL SELECTOR style baseline (docs/internal/coda_baseline_zh.md; post hoc).

Inputs (untracked): <root>/<element or locus>/summary.tsv from proxy_select.py with the frozen v2 flags plus --coda
(--repeats 50 --comparator-repeats 50). v2 = disc panel, posterior average (rho_bma); coda@kappa = consensus prior,
eig panel, the single most probable track (rho_pick), and as an ablation its posterior average (rho_bma).

usage: python coda_report.py dev   <root> <out_dir>               chooses kappa (declared rule) -> coda_dev_kappa.txt
       python coda_report.py test  <root> <out_dir> <kappa>       Kircher test regions
       python coda_report.py siraj <root> <out_dir> <kappa>       Siraj loci (two cell types averaged)
Outputs: <out_dir>/coda_<mode>_levels.tsv, coda_<mode>_contrasts.tsv (test, siraj)
"""
import os, sys, glob
import numpy as np
import pandas as pd

mode, root, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)
MARGIN, TIE, B = 0.02, 0.005, 10000
s = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(root, "*", "summary.tsv")))])
if mode == "siraj":
    s["unit"] = s.element.str.rsplit(".", n=1).str[0]
else:
    DEV = ["PKLR-48h", "F9", "HNF4A", "IRF6", "BCL11A"]
    s = s[s.element.isin(DEV) == (mode == "dev")]
    var = pd.read_csv("pilot/data/variants_grch38.tsv", sep="\t", low_memory=False)
    reg = var.groupby("element").agg(c=("chrom", "first"), a=("pos", "min"), b=("pos", "max"))
    s["unit"] = s.element.map(reg.c + ":" + reg.a.astype(str) + "-" + reg.b.astype(str))


def val(view, strategy, k, col):
    q = s[(s.scenario == view) & (s.strategy == strategy) & (s.k == k)]
    return q.groupby("unit")[col].mean()


coda = sorted({x for x in s.strategy if x.startswith("coda@")}, key=lambda x: float(x[5:]))
views = [v for v in ["all", "hide_exact", "hide_matched"] if v in set(s.scenario)]
lv = [dict(view=v, k=k, method=m, output=c, units=len(x), rho=x.mean())
      for v in views for k in [8, 16] for m, c in [("disc", "rho_bma")] + [(m, c) for m in coda for c in
                                                                         ["rho_pick", "rho_bma"]]
      for x in [val(v, m, k, c).dropna()]]
lv = pd.DataFrame(lv)
lv.round(4).to_csv(os.path.join(out, f"coda_{mode}_levels.tsv"), sep="\t", index=False)

if mode == "dev":
    # declared rule: mean rho_pick at k = 8 over C1 and C3; ties within TIE go to the smallest kappa
    m = lv[(lv.k == 8) & (lv.output == "rho_pick") & (lv.method != "disc")].groupby("method").rho.mean()
    tied = [x for x in coda if m[x] >= m.max() - TIE]
    text = "\n".join(f"{x}: mean rho_pick k=8 over C1,C3 = {m[x]:.4f}" for x in coda) + \
        f"\nchosen: {tied[0]} (best {m.idxmax()}, ties within {TIE} -> smallest kappa)\n"
    open(os.path.join(out, "coda_dev_kappa.txt"), "w").write(text)
    print(lv.round(3).to_string(index=False))
    print(text)
    sys.exit(0)

kap = f"coda@{float(sys.argv[4]):g}"
rows = []
pairs = [("v2 - coda (single track)", ("disc", "rho_bma"), (kap, "rho_pick")),
         ("v2 - model_selector (single track)", ("disc", "rho_bma"), ("coda@0", "rho_pick")),
         ("coda: average - single track", (kap, "rho_bma"), (kap, "rho_pick")),
         ("v2 - coda (average)", ("disc", "rho_bma"), (kap, "rho_bma"))]
for v in views:
    for k in [8, 16]:
        for name, (m1, c1), (m2, c2) in pairs:
            d = (val(v, m1, k, c1) - val(v, m2, k, c2)).dropna()
            bm = np.random.default_rng(0).choice(d.to_numpy(), (B, len(d)), replace=True).mean(1)
            lo, hi = np.quantile(bm, [0.025, 0.975])
            # the reading rule of cmp_report.py
            reading = (f"equivalent within {MARGIN}" if lo > -MARGIN and hi < MARGIN else "first better" if lo > 0
                       else "second better" if hi < 0 else "undetermined")
            rows.append(dict(view=v, k=k, contrast=name, units=len(d), mean=d.mean(), ci_low=lo, ci_high=hi,
                             positive=int((d > 0).sum()), reading=reading))
rows = pd.DataFrame(rows)
rows.round(4).to_csv(os.path.join(out, f"coda_{mode}_contrasts.tsv"), sep="\t", index=False)
with pd.option_context("display.width", 200, "display.max_rows", 200):
    print(lv.round(3).to_string(index=False))
    print(rows.round(4).to_string(index=False))
