"""Load one element of the paper's data the way pilot/scripts/proxy_select.py does, as TrackPilot inputs.

The data are local only (pilot/data is gitignored; the Siraj data are CC-BY-NC and AlphaGenome outputs fall under its
terms of use), so the tests that use them skip when the files are absent.
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from trackpilot import TrackPilot, match_tracks, read_predictions

ROOT = os.path.join(os.path.dirname(__file__), "..", "pilot")
SIRAJ = dict(variants="data/siraj_ref/proxy_variants.tsv", pred_dir="predictions/siraj_all",
             runs="results/proxy_select_siraj_ref/{locus}/runs.tsv.gz",
             cmp="results/proxy_select_siraj_ref_cmp/{locus}/runs.tsv.gz")
_var = {}


def variants(path):
    if path not in _var:
        v = pd.read_csv(os.path.join(ROOT, path), sep="\t", low_memory=False)
        v["key"] = v.chrom + ":" + v.pos.astype(str) + ":" + v.ref + ">" + v.alt
        if "pred_element" not in v:
            v["pred_element"] = v.element
        _var[path] = v
    return _var[path]


def load(element, variants_path, pred_dir, version="ALL_FOLDS"):
    """Predictions (rows in label order), labels, and the matched / related track masks."""
    var = variants(variants_path)
    lab = var[var.element == element]
    X, meta = read_predictions(os.path.join(ROOT, pred_dir, f"{lab.pred_element.iloc[0]}__{version}.tsv.gz"))
    lab = lab[lab.key.isin(X.index)]
    X = X.loc[lab.key]
    matched = match_tracks(meta, [k for k in lab.model_biosample_keywords.iloc[0].split(";") if k])
    rel = lab.related_biosamples.iloc[0] if "related_biosamples" in lab else ""
    related = match_tracks(meta, biosamples=[k for k in str(rel).split(";") if k]) & ~matched
    y = pd.Series(lab.effect.to_numpy(float), index=X.index)
    return X, y, matched, related


def scenarios(X, matched, related):
    """The paper's three scenarios: profiled (all), related-only (hide_exact), unprofiled (hide_matched)."""
    out = {"all": TrackPilot(X, matched if matched.any() else None)}
    if matched.any():
        out["hide_matched"] = TrackPilot(X, None, matched | related)
    if matched.any() and related.any():
        out["hide_exact"] = TrackPilot(X, related, matched)
    return out


def rank_corr(x, y):
    r1, r2 = rankdata(x), rankdata(y)
    r1, r2 = r1 - r1.mean(), r2 - r2.mean()
    return float((r1 * r2).sum() / np.sqrt((r1 ** 2).sum() * (r2 ** 2).sum()))


def available(*paths):
    return all(os.path.exists(os.path.join(ROOT, p)) for p in paths)
