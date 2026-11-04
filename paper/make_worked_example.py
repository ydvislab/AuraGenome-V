"""Figure 4, worked example: IRF4 enhancer in SK-MEL-28 (Kircher test set), unprofiled scenario, repeat 0, AlphaGenome.
(a) measured effects along the element with the pilot panel; (b) posterior over track families after 0, 4 and 8
measurements; (c) held-out variants, measured against predicted. Run from the repository root (package installed):
    python paper/make_worked_example.py      -> paper/figures/fig_worked_example.pdf
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "tests")
from paper_data import ROOT, load, rank_corr, variants  # noqa: E402
from trackpilot import TrackPilot, read_predictions  # noqa: E402

E, REP, KS = "IRF4", 0, [0, 4, 8]
X, y, matched, related = load(E, "data/variants_grch38.tsv", "predictions/alphagenome_all")
_, meta = read_predictions(os.path.join(ROOT, "predictions/alphagenome_all", f"{E}__ALL_FOLDS.tsv.gz"))
tp = TrackPilot(X, None, matched | related)
n = len(y)
pool = np.random.default_rng([0, REP]).choice(n, int(0.5 * n), replace=False)
out = np.ones(n, bool)
out[pool] = False
panel = tp.suggest(8, y.index[pool], oracle=y)
fit = tp.fit(y.loc[panel])
pred, default = fit.predict(), tp.prior_prediction()
rho_tp, rho_def = rank_corr(pred[out], y[out]), rank_corr(default[out], y[out])
p_check = tp.signal_check(y.loc[panel])

# family weights after k measurements; a family is named by its most common biosample and output type
fam = pd.Series(tp.families[:len(tp.tracks)], index=tp.tracks)
m = meta.reindex(tp.tracks)
label = m.biosample_name.fillna("?").groupby(fam.to_numpy()).agg(
    lambda s: s.value_counts().index[0] + (f" (+{len(s) - 1})" if len(s) > 1 else ""))
W = {k: pd.Series(np.exp(tp.log_prior) if k == 0 else tp.fit(y.loc[panel[:k]]).weights,
                  index=tp.hypotheses).groupby(tp.families).sum() for k in KS}
top = W[8].sort_values(ascending=False).index[:4]

pos = variants("data/variants_grch38.tsv")
pos = pos[pos.element == E].set_index("key").pos.reindex(y.index)
pos = pos - pos.min() + 1
plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.0), gridspec_kw=dict(width_ratios=[1.35, 1.0, 0.9]), layout="constrained")
a = ax[0]
a.scatter(pos, y, s=2, c="0.75", lw=0)
a.scatter(pos[panel], y[panel], s=22, c="#c2410c", lw=0.5, edgecolor="white", zorder=3)
for i, v in enumerate(panel):
    a.annotate(str(i + 1), (pos[v], y[v]), xytext=(2, 2), textcoords="offset points", fontsize=5.5, color="#c2410c")
a.axhline(0, c="0.4", lw=0.5)
a.set_xlabel("position in the element (bp)")
a.set_ylabel("measured effect (log$_2$)")
a.set_title("a  pilot panel (orange, in order)", loc="left", fontsize=7)
a = ax[1]
w = 0.26
for j, k in enumerate(KS):
    a.barh(np.arange(len(top)) + (j - 1) * w, W[k].reindex(top).to_numpy(), height=w,
           color=["0.8", "#fdba74", "#c2410c"][j], label=f"k = {k}")
a.set_yticks(np.arange(len(top)))
a.set_yticklabels([label[t][:28] for t in top], fontsize=6)
a.invert_yaxis()
a.set_xlabel("posterior weight")
a.set_xlim(0, 1.45 * W[8].max())
a.legend(frameon=False, fontsize=6, loc="center right")
a.set_title("b  posterior over track families", loc="left", fontsize=7)
a = ax[2]
a.scatter(pred[out], y[out], s=2, c="#c2410c", lw=0, alpha=0.6)
a.set_xlabel("prediction, k = 8")
a.set_ylabel("measured effect")
a.text(0.03, 0.97, f"held-out variants\nTrackPilot $\\rho$ = {rho_tp:.2f}\nall-track mean $\\rho$ = {rho_def:.2f}\n"
       f"check P = {p_check:.3f}", transform=a.transAxes, va="top", fontsize=6)
a.set_title("c  prediction", loc="left", fontsize=7)
fig.savefig("paper/figures/fig_worked_example.pdf")
print(f"rho TrackPilot {rho_tp:.3f}, all-track mean {rho_def:.3f}, check P {p_check:.3f}; top families:")
for t in top:
    print(" ", label[t], [round(W[k][t], 3) for k in KS])
