"""TrackPilot reproduces the pre-registered Siraj runs of the paper number for number.

For a few loci, both cell types, all three scenarios and the first repeats, the package recomputes the variant pool,
the sequential pilot panel, the one-round panel and the posterior-average prediction, and compares the held-out
Spearman correlation with the rows written by pilot/scripts/proxy_select.py (frozen method v2). The full
reproduction of P1-P4 is pilot/scripts/trackpilot_reproduce.py.
"""
import os

import numpy as np
import pandas as pd
import pytest

from paper_data import ROOT, SIRAJ, available, load, rank_corr, scenarios

LOCI = ["10_121010256_A_G_wC", "19_1493207_C_T_wC", "7_27717716_G_C_wC"]
REPEATS = 3
BUDGETS = [2, 4, 8, 12, 16]

pytestmark = pytest.mark.skipif(not available(SIRAJ["variants"], SIRAJ["pred_dir"]),
                                reason="paper data not present (pilot/data, pilot/predictions are local only)")


def held_out_rho(tp, y, panel, outside):
    pred = tp.fit(y.loc[panel]).predict()
    return rank_corr(pred[outside], y[outside])


@pytest.mark.parametrize("locus", LOCI)
def test_siraj_runs(locus):
    runs = pd.read_csv(os.path.join(ROOT, SIRAJ["runs"].format(locus=locus)), sep="\t")
    cmp_path = os.path.join(ROOT, SIRAJ["cmp"].format(locus=locus))
    cmp = pd.read_csv(cmp_path, sep="\t") if os.path.exists(cmp_path) else None
    n_checked = 0
    for cell in ["K562", "HEPG2"]:
        element = f"{locus}.{cell}"
        X, y, matched, related = load(element, SIRAJ["variants"], SIRAJ["pred_dir"])
        n = len(y)
        for scen, tp in scenarios(X, matched, related).items():
            for rep in range(REPEATS):
                rng = np.random.default_rng([0, rep])
                pool = rng.choice(n, int(0.5 * n), replace=False)
                outside = np.ones(n, bool)
                outside[pool] = False
                ids = y.index[pool]
                panels = {"disc": tp.suggest(16, ids, oracle=y), "batch_kb": tp.suggest(16, ids)}
                for strategy, order in panels.items():
                    ref = runs if strategy == "disc" else cmp
                    if ref is None:
                        continue
                    sel = ref[(ref.element == element) & (ref.scenario == scen) & (ref.strategy == strategy)
                              & (ref.repeat == rep)].set_index("k").rho_bma
                    for k in BUDGETS:
                        assert held_out_rho(tp, y, order[:k], outside) == pytest.approx(sel.loc[k], abs=1e-9), \
                            (element, scen, strategy, rep, k)
                        n_checked += 1
                # label-free references: the anchor (matched or related mean) and the all-track mean
                sel = runs[(runs.element == element) & (runs.scenario == scen) & (runs.repeat == rep)
                           & (runs.k == 0)].set_index("strategy").rho_pick
                prior = tp.prior_prediction()
                name = {"all": "matched", "hide_exact": "related"}.get(scen, "all_mean")
                assert rank_corr(prior[outside], y[outside]) == pytest.approx(sel.loc[name], abs=1e-9)
    assert n_checked > 0
