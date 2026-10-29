"""C1 on Siraj: change of each pre-registered contrast D_i between ALL_FOLDS and the held-out-fold checkpoint
(docs/internal/c1_robustness_zh.md, "扩展到 Siraj"). Unit: the locus; paired over loci present in both tables.

usage: python c1_siraj_compare.py <all_folds_by_locus.tsv> <fold_by_locus.tsv> <out.tsv>
       e.g. results/a1/siraj_by_locus.tsv results/a1/c1_siraj/siraj_by_locus.tsv results/a1/c1_siraj/c1_siraj_delta.tsv
Output columns: contrast, n, mean D under ALL_FOLDS, mean D held out, mean change (held out - ALL_FOLDS) with a
95% locus-bootstrap interval.
"""
import sys
import numpy as np
import pandas as pd

a = pd.read_csv(sys.argv[1], sep="\t", index_col="locus")
f = pd.read_csv(sys.argv[2], sep="\t", index_col="locus")
loci = a.index.intersection(f.index)
rng = np.random.default_rng(0)
idx = rng.integers(0, len(loci), size=(10000, len(loci)))
rows = []
for c in [c for c in a.columns if c in f.columns]:
    d = (f.loc[loci, c] - a.loc[loci, c]).to_numpy()
    bm = d[idx].mean(axis=1)
    rows.append(dict(contrast=c, n=len(loci), all_folds=a.loc[loci, c].mean(), held_out=f.loc[loci, c].mean(),
                     change=d.mean(), ci_low=np.quantile(bm, 0.025), ci_high=np.quantile(bm, 0.975)))
t = pd.DataFrame(rows)
t.to_csv(sys.argv[3], sep="\t", index=False, float_format="%.4f")
print(t.to_string(index=False, float_format=lambda x: f"{x:+.4f}"))
