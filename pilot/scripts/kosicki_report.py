"""Pre-registered analysis of the neuron MPRA (docs/internal/neuron_mpra_zh.md, section 5).

Inputs: pilot/results/proxy_select_kosicki_<model>/<region>/summary.tsv (one proxy_select run per region, all tiles of
the region) and pilot/data/kosicki/proxy_variants_<model>.tsv (tile -> region). Unit: the region; tiles averaged first.

usage (from the repository root): python pilot/scripts/kosicki_report.py <model>      model: alphagenome | enformer
Output: results/a1/kosicki_<model>/{decision.txt, by_region.tsv, levels.tsv}
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, beta

model = sys.argv[1]
root = f"pilot/results/proxy_select_kosicki_{model}"
out = f"results/a1/kosicki_{model}"
os.makedirs(out, exist_ok=True)
s = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{root}/*/summary.tsv"))])
var = pd.read_csv(f"pilot/data/kosicki/proxy_variants_{model}.tsv", sep="\t", usecols=["element", "pred_element"])
s["region"] = s.element.map(var.drop_duplicates("element").set_index("element").pred_element)
DELTA, TAU, Q, B = 0.02, 0.05, 0.15, 10000


def val(scen, strategy, k, col="rho_bma"):
    q = s[(s.scenario == scen) & (s.strategy == strategy) & (s.k == k)]
    return q.groupby("region")[col].mean()


lines, rows = [f"neuron MPRA, {model}: {s.region.nunique()} regions, {s.element.nunique()} tiles"], {}
D1 = (val("hide_matched", "disc", 8) - val("hide_matched", "all_mean", 0, "rho_pick")).dropna()
p = wilcoxon(D1, alternative="greater").pvalue
n1 = D1.mean() >= 0.03 and p < 0.05
lines.append(f"N1 C3 k=8: TrackPilot - all-track mean = {D1.mean():+.4f} (>= 0.03?), {int((D1 > 0).sum())}/{len(D1)} > 0, "
             f"one-sided Wilcoxon p = {p:.3g} -> {'HOLDS' if n1 else 'FAILS'}")
D2 = (val("hide_matched", "disc", 8) - val("hide_matched", "random", 8)).dropna()
p = wilcoxon(D2, alternative="greater").pvalue
lines.append(f"N2 C3 k=8: pilot - random panel = {D2.mean():+.4f}, {int((D2 > 0).sum())}/{len(D2)} > 0, "
             f"one-sided Wilcoxon p = {p:.3g} -> {'HOLDS' if p < 0.05 else 'FAILS'}")
D3 = (val("all", "disc", 8) - val("all", "matched", 0, "rho_pick")).dropna()
lb = np.quantile(np.random.default_rng(0).choice(D3.to_numpy(), (B, len(D3)), replace=True).mean(1), 0.05)
x, n = int((D3 < -TAU).sum()), len(D3)
ub = beta.ppf(0.95, x + 1, n - x) if x < n else 1.0
n3 = lb > -DELTA and ub <= Q
lines.append(f"N3 C1 k=8: TrackPilot - matched mean = {D3.mean():+.4f}, bootstrap lower bound {lb:+.4f} (> -{DELTA}?), "
             f"harmed {x}/{n}, upper bound {ub:.3f} (<= {Q}?) -> {'HOLDS' if n3 else 'FAILS'}")
D4 = (val("hide_exact", "disc", 8) - val("hide_exact", "related", 0, "rho_pick")).dropna()
p = wilcoxon(D4, alternative="greater").pvalue
lines.append(f"N4 C2 k=8: TrackPilot - related mean = {D4.mean():+.4f}, {int((D4 > 0).sum())}/{len(D4)} > 0, "
             f"one-sided Wilcoxon p = {p:.3g} -> {'HOLDS' if p < 0.05 else 'FAILS'}")
lines.append(f"DECISION: {'HOLDS in a poorly profiled cell type' if n1 and n3 else 'NOT SUPPORTED'} (N1 and N3)")

lev = []
for scen, ref in [("all", "matched"), ("hide_exact", "related"), ("hide_matched", "all_mean")]:
    base, orc = val(scen, ref, 0, "rho_pick"), val(scen, "oracle", 0, "rho_pick")
    for k in sorted(s[s.strategy == "disc"].k.unique()):
        lev.append(dict(scenario=scen, k=k, pilot=val(scen, "disc", k).mean(), random=val(scen, "random", k).mean(),
                        default=base.mean(), oracle=orc.mean()))
text = "\n".join(lines) + "\n"
open(f"{out}/decision.txt", "w").write(text)
pd.DataFrame({"N1_C3": D1, "N2_C3": D2, "N3_C1": D3, "N4_C2": D4}).round(4).to_csv(f"{out}/by_region.tsv", sep="\t")
pd.DataFrame(lev).round(4).to_csv(f"{out}/levels.tsv", sep="\t", index=False)
print(text)
print(pd.DataFrame(lev).round(3).to_string(index=False))
