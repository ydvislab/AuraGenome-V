"""Command-line interface.

  trackpilot suggest  PREDICTIONS -k 8 --anchor K562 [--pool pool.txt] [-o panel.tsv]
  trackpilot predict  PREDICTIONS LABELS --anchor K562 [-o predictions.tsv] [--posterior tracks.tsv] [--report r.html]
  trackpilot score-alphagenome --region chr11:5227000-5228000 -o predictions.tsv.gz    (needs ALPHAGENOME_API_KEY)
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

from . import __version__
from .core import ANCHOR, SIGNAL_P, TrackPilot
from .io import OUTPUTS, match_tracks, read_labels, read_predictions, split_list


def _common(sp):
    sp.add_argument("predictions", help="predicted effects: long TSV (variant_id, track_name, raw_score, ...) or wide")
    sp.add_argument("--anchor", default="", help="keywords naming the cell type in track metadata, ';'-separated "
                    "(e.g. 'K562' or 'liver;hepatocyte'); their tracks' mean becomes the anchor hypothesis")
    sp.add_argument("--anchor-biosamples", default="", help="exact biosample names added to the anchor")
    sp.add_argument("--exclude", default="", help="keywords of tracks to leave out (simulate an unprofiled cell type)")
    sp.add_argument("--outputs", default=",".join(OUTPUTS), help="output types to use (long input only)")
    sp.add_argument("--anchor-prior", type=float, default=0.95, help="prior mass of the anchor (paper: 0.95)")


def _pilot(a):
    t0 = time.time()
    X, meta = read_predictions(a.predictions, split_list(a.outputs))
    anchor = match_tracks(meta, split_list(a.anchor), split_list(a.anchor_biosamples))
    exclude = match_tracks(meta, split_list(a.exclude)) if a.exclude else None
    tp = TrackPilot(X, anchor if anchor.any() else None, exclude, anchor_prior=a.anchor_prior)
    T = len(tp.tracks)
    msg = f"{len(tp.variants)} variants, {T} tracks in {len(np.unique(tp.families[:T]))} families"
    msg += (f"; anchor = mean of {tp.n_anchor_tracks} matched tracks" if tp.n_anchor_tracks else
            f"; no anchor ({'no track matches --anchor' if a.anchor or a.anchor_biosamples else 'no --anchor given'}):"
            " equal prior weight on each track family")
    print(msg, file=sys.stderr)
    return tp, meta, t0


def suggest(a):
    tp, _, t0 = _pilot(a)
    pool = None
    if a.pool:
        pool = [v for v in pd.read_csv(a.pool, sep="\t", header=None)[0].astype(str) if v != "variant_id"]
    measured = read_labels(a.measured) if a.measured else None
    panel = tp.suggest(a.k, pool, measured)
    out = pd.DataFrame(dict(order=np.arange(1, len(panel) + 1), variant_id=panel))
    out.to_csv(a.out or sys.stdout, sep="\t", index=False)
    print(f"pilot panel of {len(panel)} variants ({time.time() - t0:.1f}s)", file=sys.stderr)


def predict(a):
    tp, meta, t0 = _pilot(a)
    y = read_labels(a.labels, a.column)
    y = y[y.index.isin(tp.variants)]
    fit = tp.fit(y)
    p = tp.signal_check(y)
    fallback = p > a.check_p
    checked = (tp.prior_prediction() if fallback else fit.predict()).rename("prediction_checked")
    out = pd.concat([fit.predict(), checked, fit.predict_no_offset(), fit.predict_map(), tp.prior_prediction(),
                     y.rename("measured")], axis=1)
    # the fallback is the label-free default in standardised track units: a ranking score, not a calibrated effect
    out["checked_source"] = "default_ranking_score" if fallback else "posterior_average"
    out.index.name = "variant_id"
    out.to_csv(a.out or sys.stdout, sep="\t", float_format="%.6g")
    post = posterior_table(tp, fit, meta)
    if a.posterior:
        post.to_csv(a.posterior, sep="\t", index=False, float_format="%.6g")
    if a.report:
        from .report import write_report
        write_report(a.report, tp, fit, post, out)
    verdict = ("the labels fit no track better than shuffled labels: prediction_checked is the label-free default, a "
               "ranking score in standardised track units, not a calibrated effect" if fallback
               else "the labels fit the tracks better than shuffled labels: prediction_checked is the posterior average")
    print(summary(tp, fit, post) + f"\napplicability check P = {p:.2f} ({verdict})\n({time.time() - t0:.1f}s)",
          file=sys.stderr)


def posterior_table(tp, fit, meta):
    m = meta.reindex(tp.hypotheses)
    if ANCHOR in tp.hypotheses:
        m.loc[ANCHOR, "biosample_name"] = f"mean of {tp.n_anchor_tracks} matched tracks"
    post = pd.DataFrame(dict(hypothesis=tp.hypotheses, weight=fit.weights, family=tp.families,
                             intercept=fit.intercept, slope=fit.slope))
    post["family_weight"] = post.groupby("family").weight.transform("sum")
    for c in ["output_type", "biosample_name"]:
        if c in m:
            post[c] = m[c].to_numpy()
    return post.sort_values("weight", ascending=False)


def summary(tp, fit, post):
    fam = post.groupby("family").weight.sum().sort_values(ascending=False)
    top = post[post.family == fam.index[0]].iloc[0]
    lines = [f"leading track family: {top.get('biosample_name', top.hypothesis)} ({fam.iloc[0]:.2f} of the posterior;"
             f" top hypothesis {top.hypothesis})",
             f"track families holding 90% of the posterior: {int(np.searchsorted(fam.cumsum().to_numpy(), 0.9)) + 1}"]
    if ANCHOR in tp.hypotheses:
        lines.append(f"anchor weight: {float(post.weight[post.hypothesis == ANCHOR].iloc[0]):.2f}")
    return "\n".join(lines)


def score_alphagenome(a):
    from .alphagenome import score_region
    df = score_region(a.region, split_list(a.outputs), a.chunk_bp)
    df.to_csv(a.out, sep="\t", index=False)
    print(f"wrote {len(df)} rows to {a.out}", file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="trackpilot", description="Choose which genomic foundation-model tracks to "
                                 "trust for a cell type from a few measured variants.")
    ap.add_argument("--version", action="version", version=f"trackpilot {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("suggest", help="choose the pilot panel: the variants to measure")
    _common(s)
    s.add_argument("-k", type=int, default=8, help="panel size (paper: 8)")
    s.add_argument("--pool", default="", help="file of variant ids that can be measured (default: all)")
    s.add_argument("--measured", default="", help="labels already measured (variant_id, effect): suggest k more")
    s.add_argument("-o", "--out", default="")
    s.set_defaults(func=suggest)
    p = sub.add_parser("predict", help="posterior over tracks and predictions for every variant")
    _common(p)
    p.add_argument("labels", help="measured effects: TSV/CSV with variant_id and effect")
    p.add_argument("--column", default=None, help="label column (default: effect, else the second column)")
    p.add_argument("-o", "--out", default="", help="predictions per variant (default: stdout)")
    p.add_argument("--posterior", default="", help="write the posterior over tracks here")
    p.add_argument("--report", default="", help="write an HTML report here")
    p.add_argument("--check-p", type=float, default=SIGNAL_P, help="applicability check: fall back to the default "
                   f"when the permutation P exceeds this (paper: {SIGNAL_P}; 0.2 is stricter, for assays that may be "
                   "weak)")
    p.set_defaults(func=predict)
    g = sub.add_parser("score-alphagenome", help="score every substitution of a region with AlphaGenome (ISM)")
    g.add_argument("--region", required=True, help="chrom:start-end, 1-based inclusive, GRCh38")
    g.add_argument("--outputs", default=",".join(OUTPUTS))
    g.add_argument("--chunk-bp", type=int, default=100)
    g.add_argument("-o", "--out", required=True)
    g.set_defaults(func=score_alphagenome)
    a = ap.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
