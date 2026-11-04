"""Descriptive report of the fourth test (docs/internal/kircher_mouse_zh.md, section 5): how often the frozen
applicability check (p* = 0.5) calls the Kircher mouse-cell experiments 'no signal', and what the checked prediction
gives. Same runs as applicability_check.py (20 repeats, k = 8, 200 permutations); no test is made.

usage (from the repository root): python pilot/scripts/kircher_mouse_check.py <alphagenome|enformer>
Output: results/a1/kircher_mouse_<model>/check.txt
"""
import os, sys
from multiprocessing import Pool
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from applicability_check import job, unit_means  # noqa: E402

PSTAR = 0.5
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def main(model):
    os.chdir(os.path.join(ROOT, "pilot"))  # data and prediction paths are relative to pilot/
    pdir, ver = (("predictions/alphagenome_all", "ALL_FOLDS") if model == "alphagenome"
                 else ("predictions/enformer_kircher_mouse", "ENFORMER"))
    todo = [("kircher_mouse", e, "ZRS" if e.startswith("ZRS") else e, "test",
             f"data/kircher_mouse/proxy_variants_{model}.tsv", pdir, ver)
            for e in ["RET", "UC88", "TCF7L2", "ZFAND3", "ZRSh-13", "ZRSh-13h2"]]
    with Pool(6) as pool:
        d = pd.DataFrame([r for rr in pool.map(job, todo, chunksize=1) for r in rr])
    u = unit_means(d, PSTAR)
    lines = [f"Kircher mouse-cell experiments, {model}, applicability check p* = {PSTAR} (descriptive)"]
    for scen, g in u.groupby("scenario"):
        share = (d[d.scenario == scen].p > PSTAR).mean()
        lines.append(f"{scen}: share no signal {share:.2f}; mean rho TrackPilot {g.rho_tp.mean():.3f}, checked "
                     f"{g.gated.mean():.3f}, default {g.rho_default.mean():.3f} ({len(g)} regions)")
    os.makedirs(f"{ROOT}/results/a1/kircher_mouse_{model}", exist_ok=True)
    open(f"{ROOT}/results/a1/kircher_mouse_{model}/check.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(u.round(3).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1])
