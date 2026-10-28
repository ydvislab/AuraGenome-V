"""Pilot step 3: reliability signals per variant and the go/no-go test (research plan R1).

Reads the filtered model scores written by alphagenome_score.py (and borzoi_score.py when
present), summarises each variant's predicted effect on the MPRA cell type, derives the
reliability signals, and asks the pilot question:

    Among variants with a measured effect (truth_class = +1/-1), do the reliability signals
    separate predictions with the right direction from predictions with the wrong direction?

Signals (higher = more reliable)
  ag_fold_agree     fraction of AlphaGenome FOLD_0..3 whose effect sign matches their mean sign
  ag_fold_snr       |mean| / sd of the four fold effects
  ag_track_agree    fraction of matched tracks (same cell type, different experiments) agreeing in sign
  ag_magnitude      |ALL_FOLDS effect|
  bz_rep_agree      fraction of Borzoi rep0..3 agreeing in sign            (if Borzoi present)
  cross_model_agree 1 if AlphaGenome and Borzoi signs agree, else 0          (if Borzoi present)
  combined          mean of the per-signal ranks (no fitting; a learned combination comes later)
Report: AUROC of each signal for "direction correct", per-element Spearman between predicted and
measured effect, and the accuracy among the top / bottom reliability tertiles.

--simulate writes synthetic model files with a known reliability structure and runs the same
code, to test the pipeline before any API or GPU time is spent. Simulated numbers mean nothing
biologically.

usage:
  python compute_signals.py <pilot_dir> --simulate
  python compute_signals.py <pilot_dir> [--promoter-output CAGE] [--enhancer-output DNASE]
"""
import os, sys, glob, json, argparse
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("pilot_dir")
ap.add_argument("--promoter-output", default="CAGE")
ap.add_argument("--enhancer-output", default="DNASE")
ap.add_argument("--simulate", action="store_true")
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()

var = pd.read_csv(os.path.join(a.pilot_dir, "data", "variants_grch38.tsv"), sep="\t", low_memory=False)
var = var[var.human_cells & ~var.biosample_is_proxy].copy()
var["vkey"] = var.chrom + ":" + var.pos.astype(str) + ":" + var.ref + ">" + var.alt
pred_root = os.path.join(a.pilot_dir, "predictions_simulated" if a.simulate else "predictions")
VERSIONS = ["ALL_FOLDS", "FOLD_0", "FOLD_1", "FOLD_2", "FOLD_3"]

def simulate():
    """Synthetic AlphaGenome/Borzoi files: each variant has a latent reliability r; the model
    effect is truth*|signal|*r + noise, and fold/replicate noise grows as r falls."""
    rng = np.random.default_rng(a.seed)
    os.makedirs(os.path.join(pred_root, "alphagenome"), exist_ok=True)
    os.makedirs(os.path.join(pred_root, "borzoi"), exist_ok=True)
    for element, ev in var.groupby("element"):
        n = len(ev)
        r = rng.uniform(0, 1, n)
        truth = np.where(ev.truth_class.values == 0, rng.normal(0, 0.2, n), ev.effect.values)
        base = truth * r + rng.normal(0, 0.3, n) * (1 - r)
        out_type = a.promoter_output if ev.element_type.iloc[0] == "promoter" else a.enhancer_output
        kw = ev.model_biosample_keywords.iloc[0].split(";")[0]
        for v in VERSIONS:
            rows = []
            for t in range(3):  # three matched tracks
                eff = base + rng.normal(0, 0.05 + 0.6 * (1 - r), n) + rng.normal(0, 0.1, n)
                rows.append(pd.DataFrame(dict(variant_id=ev.vkey.values, output_type=out_type,
                                              track_name=f"{kw} track{t}", biosample_name=kw, raw_score=eff)))
            pd.concat(rows).to_csv(os.path.join(pred_root, "alphagenome", f"{element}__{v}.tsv.gz"), sep="\t", index=False)
        for rep in range(4):
            eff = base * 0.8 + rng.normal(0, 0.05 + 0.7 * (1 - r), n)
            pd.DataFrame(dict(variant_id=ev.variant_id.values, track_index=0, identifier="x",
                              description=f"CAGE:{kw}", raw_score=eff)).to_csv(
                os.path.join(pred_root, "borzoi", f"{element}__rep{rep}.tsv.gz"), sep="\t", index=False)

if a.simulate:
    simulate()

def load_ag(element, version, out_type):
    f = os.path.join(pred_root, "alphagenome", f"{element}__{version}.tsv.gz")
    if not os.path.exists(f):
        return None
    df = pd.read_csv(f, sep="\t")
    df = df[df.output_type == out_type]
    return df

rows = []
for element, ev in var.groupby("element"):
    out_type = a.promoter_output if ev.element_type.iloc[0] == "promoter" else a.enhancer_output
    ag = {v: load_ag(element, v, out_type) for v in VERSIONS}
    if ag["ALL_FOLDS"] is None:
        continue
    per_v = {v: df.groupby("variant_id").raw_score.mean() for v, df in ag.items() if df is not None}
    track_agree = ag["ALL_FOLDS"].groupby("variant_id").raw_score.apply(lambda s: max((s > 0).mean(), (s < 0).mean()))
    bz = [pd.read_csv(f, sep="\t").groupby("variant_id").raw_score.mean()
          for f in sorted(glob.glob(os.path.join(pred_root, "borzoi", f"{element}__rep*.tsv.gz")))]
    for v in ev.itertuples():
        if v.vkey not in per_v["ALL_FOLDS"].index:
            continue
        main = per_v["ALL_FOLDS"][v.vkey]
        folds = np.array([per_v[f][v.vkey] for f in VERSIONS[1:] if f in per_v and v.vkey in per_v[f].index])
        row = dict(variant_id=v.variant_id, element=element, element_type=v.element_type, output_type=out_type,
                   measured=v.effect, truth_class=v.truth_class, ag_effect=main, ag_magnitude=abs(main),
                   ag_track_agree=track_agree.get(v.vkey, np.nan))
        if len(folds) >= 2:
            sgn = np.sign(folds.mean())
            row["ag_fold_agree"] = float((np.sign(folds) == sgn).mean())
            row["ag_fold_snr"] = float(abs(folds.mean()) / (folds.std(ddof=1) + 1e-9))
        reps = np.array([b[v.variant_id] for b in bz if v.variant_id in b.index])
        if len(reps) >= 2:
            row["bz_effect"] = float(reps.mean())
            row["bz_rep_agree"] = float((np.sign(reps) == np.sign(reps.mean())).mean())
            row["cross_model_agree"] = float(np.sign(reps.mean()) == np.sign(main))
        rows.append(row)
sig = pd.DataFrame(rows)
if sig.empty:
    sys.exit(f"no predictions under {pred_root}; run alphagenome_score.py first (or --simulate)")

SIGNALS = [c for c in ["ag_fold_agree", "ag_fold_snr", "ag_track_agree", "ag_magnitude", "bz_rep_agree", "cross_model_agree"] if c in sig]
sig["combined"] = sig[SIGNALS].rank(pct=True).mean(axis=1)
SIGNALS.append("combined")

def auroc(score, label):
    """Mann-Whitney AUROC; label 1 = positive class."""
    s, y = np.asarray(score, float), np.asarray(label, int)
    ok = ~np.isnan(s); s, y = s[ok], y[ok]
    n1, n0 = y.sum(), (1 - y).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = pd.Series(s).rank().values
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))

eff = sig[sig.truth_class != 0].copy()
eff["correct"] = (np.sign(eff.ag_effect) == eff.truth_class).astype(int)
report = dict(
    mode="SIMULATED (pipeline test only)" if a.simulate else "real predictions",
    n_variants_scored=int(len(sig)), n_with_measured_effect=int(len(eff)),
    direction_accuracy=float(eff.correct.mean()),
    auroc_direction_correct={s: auroc(eff[s], eff.correct) for s in SIGNALS},
    accuracy_by_combined_tertile={k: float(g.correct.mean()) for k, g in eff.groupby(pd.qcut(eff.combined, 3, labels=["low", "mid", "high"]), observed=True)},
    spearman_by_element={e: float(g.ag_effect.corr(g.measured, method="spearman")) for e, g in sig.groupby("element")},
)
# element-level bootstrap: variants within an element are not independent, so resample elements
rng = np.random.default_rng(a.seed + 1)
els = eff.element.unique()
boot = []
for _ in range(1000):
    pick = rng.choice(els, len(els), replace=True)
    b = pd.concat([eff[eff.element == e] for e in pick])
    boot.append(auroc(b.combined, b.correct))
lo, hi = np.nanpercentile(boot, [2.5, 97.5])
tert = report["accuracy_by_combined_tertile"]
report["combined_auroc_95ci_by_element"] = [float(lo), float(hi)]
report["go_no_go"] = dict(
    rule="GO if combined AUROC 95% CI (element bootstrap) excludes 0.5 AND lower bound >= 0.6 AND "
         "high-tertile accuracy exceeds low-tertile accuracy by >= 0.15",
    auroc_ci_ok=bool(lo >= 0.6), tertile_gap=float(tert.get("high", np.nan) - tert.get("low", np.nan)),
    decision="GO" if (lo >= 0.6 and tert.get("high", 0) - tert.get("low", 0) >= 0.15) else "NO-GO / revise")
out = os.path.join(a.pilot_dir, "results_simulated" if a.simulate else "results")
os.makedirs(out, exist_ok=True)
sig.to_csv(os.path.join(out, "signals.tsv"), sep="\t", index=False)
json.dump(report, open(os.path.join(out, "r1_report.json"), "w"), indent=1)
print(json.dumps(report, indent=1))
