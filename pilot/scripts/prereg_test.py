"""Pre-registered A1 test (docs/internal/preregistration_A1_zh.md), written before the test data existed.

Reads the proxy_select summary of the held-out experiments, averages experiments that share a genomic
region (they share one set of predictions), and evaluates at k = 8 on the posterior-weighted
prediction (rho_bma), scenario hide_matched where it exists (else all):
  H1  disc > random               one-sided Wilcoxon signed-rank over regions, p < 0.05
  H2  disc - all_mean >= 0.03      on average, and positive in >= 7 of 11 regions
  H3  disc > high_effect           one-sided Wilcoxon over regions, p < 0.05 (secondary)
Decision: H1 and H2 -> GO; H1 only -> reframe; not H1 -> drop the A1 claim.

usage:
  python proxy_select.py <pilot_dir> --elements <test elements> --tag test
  python prereg_test.py <pilot_dir>
  python prereg_test.py <pilot_dir> ce_test     # same test on results/proxy_select_ce_test (common evaluation
                                                # set, added 2026-10-29; not part of the pre-registration)
"""
import os, sys
import pandas as pd
from scipy.stats import wilcoxon

pilot_dir = sys.argv[1]
tag = sys.argv[2] if len(sys.argv) > 2 else "test"
K = 8
s = pd.read_csv(os.path.join(pilot_dir, "results", f"proxy_select_{tag}", "summary.tsv"), sep="\t")
var = pd.read_csv(os.path.join(pilot_dir, "data", "variants_grch38.tsv"), sep="\t", low_memory=False)
reg = var.groupby("element").agg(chrom=("chrom", "first"), start=("pos", "min"), end=("pos", "max"))
region = (reg.chrom + ":" + reg.start.astype(str) + "-" + reg.end.astype(str)).rename("region")

# primary scenario per experiment
scen = s.groupby("element").scenario.agg(lambda x: "hide_matched" if (x == "hide_matched").any() else "all")
s = s[s.scenario == s.element.map(scen)]
row = lambda strategy, k: s[(s.strategy == strategy) & (s.k == k)].set_index("element")
exp = pd.DataFrame({"disc": row("disc", K).rho_bma, "random": row("random", K).rho_bma,
                    "high_effect": row("high_effect", K).rho_bma, "all_mean": row("all_mean", 0).rho_pick,
                    "oracle": row("oracle", 0).rho_pick, "scenario": scen})
exp["region"] = exp.index.map(region)
by_region = exp.drop(columns="scenario").groupby("region").agg(
    lambda x: x.mean() if x.dtype.kind == "f" else ",".join(x))
by_region["elements"] = exp.reset_index().groupby("region").element.agg(",".join)

with pd.option_context("display.width", 200):
    print("experiments (k = 8, rho_bma; all_mean and oracle are label-free references)")
    print(exp.round(3).to_string())
    print(f"\nregions: {len(by_region)}")
    print(by_region[["elements", "disc", "random", "high_effect", "all_mean", "oracle"]].round(3).to_string())

d = by_region
h1 = wilcoxon(d.disc, d.random, alternative="greater")
gap = d.disc - d.all_mean
need = 7 if len(d) == 11 else -(-len(d) * 7 // 11)   # 7 of 11, scaled if a region is missing
h2 = gap.mean() >= 0.03 and (gap > 0).sum() >= need
h3 = wilcoxon(d.disc, d.high_effect, alternative="greater")
print(f"\nH1 disc > random:       mean diff {(d.disc - d.random).mean():+.3f}, "
      f"{(d.disc > d.random).sum()}/{len(d)} regions, Wilcoxon p = {h1.pvalue:.4f} -> {h1.pvalue < 0.05}")
print(f"H2 disc - all_mean:     mean {gap.mean():+.3f}, positive in {(gap > 0).sum()}/{len(d)} "
      f"(need mean >= 0.03 and >= {need}) -> {h2}")
print(f"H3 disc > high_effect:  mean diff {(d.disc - d.high_effect).mean():+.3f}, Wilcoxon p = {h3.pvalue:.4f} "
      f"-> {h3.pvalue < 0.05}")
decision = "GO" if h1.pvalue < 0.05 and h2 else ("REFRAME (H1 only)" if h1.pvalue < 0.05 else "DROP A1 claim")
print(f"DECISION: {decision}")
