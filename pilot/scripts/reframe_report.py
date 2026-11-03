"""Descriptive views of the pre-registered Siraj runs (post hoc, no new runs): label efficiency and where the gain is.

Reads pilot/results/proxy_select_siraj_ref/<locus>/summary.tsv (labels relative to the reference oligo) and writes
  results/a1/reframe/label_efficiency.tsv  mean rho_bma of pilot and random panels per k, per scenario, with the
                                           label-free default and the best single track chosen with all labels
                                           (oracle); 'gap_closed' = (pilot - default) / (oracle - default)
  results/a1/reframe/gain_by_gap.tsv       gain of the pilot panel over the default at k = 8, by quartile of the
                                           locus's gap between the default and the oracle (the oracle uses all
                                           labels, so this is descriptive only)
usage (from the repository root): python pilot/scripts/reframe_report.py
"""
import os, glob
import pandas as pd

s = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob("pilot/results/proxy_select_siraj_ref/*/summary.tsv"))])
s["locus"] = s.element.str.rsplit(".", n=1).str[0]
out = "results/a1/reframe"
os.makedirs(out, exist_ok=True)


def v(scen, strategy, k, col="rho_bma"):
    q = s[(s.scenario == scen) & (s.strategy == strategy) & (s.k == k)]
    return q.groupby("locus")[col].mean()


eff, gap_rows = [], []
for scen, ref in [("all", "matched"), ("hide_exact", "related"), ("hide_matched", "all_mean")]:
    base, orc = v(scen, ref, 0, "rho_pick"), v(scen, "oracle", 0, "rho_pick")
    for k in [2, 4, 8, 12, 16]:
        d, r = v(scen, "disc", k), v(scen, "random", k)
        eff.append(dict(scenario=scen, k=k, pilot=d.mean(), random=r.mean(), default=base.mean(), oracle=orc.mean(),
                        gap_closed=(d - base).mean() / (orc - base).mean()))
    gap, gain = orc - base, v(scen, "disc", 8) - base
    q = pd.qcut(gap, 4, labels=[1, 2, 3, 4])
    for i in [1, 2, 3, 4]:
        m = q == i
        gap_rows.append(dict(scenario=scen, quartile=i, loci=int(m.sum()), gap=gap[m].mean(), gain_k8=gain[m].mean(),
                             share_closed=gain[m].mean() / gap[m].mean(), worse_by_005=int((gain[m] < -0.05).sum())))
pd.DataFrame(eff).round(4).to_csv(f"{out}/label_efficiency.tsv", sep="\t", index=False)
pd.DataFrame(gap_rows).round(4).to_csv(f"{out}/gain_by_gap.tsv", sep="\t", index=False)
print(pd.DataFrame(eff).round(3).to_string(index=False))
print(pd.DataFrame(gap_rows).round(3).to_string(index=False))
