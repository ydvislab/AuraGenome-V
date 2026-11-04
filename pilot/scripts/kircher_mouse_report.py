"""Pre-registered analysis of the Kircher mouse-cell experiments (docs/internal/kircher_mouse_zh.md, section 5).

Inputs: pilot/results/proxy_select_kircher_mouse_<model>/<experiment>/summary.tsv. Unit: the region (the two ZRS
experiments averaged), 5 regions. C1 = 'all' (anchor: the proxy tracks), C3 = 'hide_matched' (proxy hidden, no anchor).

usage (from the repository root): python pilot/scripts/kircher_mouse_report.py <alphagenome|enformer>
Output: results/a1/kircher_mouse_<model>/{decision.txt, by_region.tsv}
"""
import os, sys, glob
import pandas as pd

model = sys.argv[1]
out = f"results/a1/kircher_mouse_{model}"
os.makedirs(out, exist_ok=True)
s = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"pilot/results/proxy_select_kircher_mouse_{model}/*/summary.tsv")])
s["region"] = s.element.str.replace(r"^ZRS.*", "ZRS", regex=True)


def val(scen, strategy, k, col="rho_bma"):
    q = s[(s.scenario == scen) & (s.strategy == strategy) & (s.k == k)]
    return q.groupby("region")[col].mean()


d = pd.DataFrame({"tp_c3": val("hide_matched", "disc", 8), "all_mean": val("hide_matched", "all_mean", 0, "rho_pick"),
                  "random_c3": val("hide_matched", "random", 8), "tp_c1": val("all", "disc", 8),
                  "proxy_mean": val("all", "matched", 0, "rho_pick"), "oracle": val("all", "oracle", 0, "rho_pick")})
m1, m2, m3 = d.tp_c3 - d.all_mean, d.tp_c3 - d.random_c3, d.tp_c1 - d.proxy_mean
h1 = m1.mean() >= 0.03 and (m1 > 0).sum() >= 4
h2 = (m2 > 0).sum() >= 4
h3 = m3.mean() >= -0.02 and (m3 >= -0.05).all()
lines = [f"Kircher mouse-cell experiments, {model}: {len(d)} regions",
         f"M1 C3 k=8: TrackPilot - all-track mean = {m1.mean():+.4f} (>= 0.03?), {int((m1 > 0).sum())}/{len(m1)} > 0 "
         f"(>= 4?) -> {'HOLDS' if h1 else 'FAILS'}",
         f"M2 C3 k=8: pilot - random panel = {m2.mean():+.4f}, {int((m2 > 0).sum())}/{len(m2)} > 0 (>= 4?) -> "
         f"{'HOLDS' if h2 else 'FAILS'}",
         f"M3 C1 k=8: TrackPilot - proxy mean = {m3.mean():+.4f} (>= -0.02?), worst region {m3.min():+.4f} (>= -0.05?) "
         f"-> {'HOLDS' if h3 else 'FAILS'}",
         f"DECISION: {'SUPPORTED' if h1 and h3 else 'NOT SUPPORTED'} (M1 and M3)"]
text = "\n".join(lines) + "\n"
open(f"{out}/decision.txt", "w").write(text)
d.round(4).to_csv(f"{out}/by_region.tsv", sep="\t")
print(text)
print(d.round(3).to_string())
