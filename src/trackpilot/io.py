"""Reading predictions and labels, and choosing anchor tracks from metadata."""
from __future__ import annotations

import numpy as np
import pandas as pd

OUTPUTS = ("CAGE", "DNASE", "ATAC")


def read_predictions(path_or_df, outputs=OUTPUTS):
    """Predicted effects as (variants x tracks DataFrame, track metadata DataFrame).

    Input: a long table with one row per variant and track, columns variant_id, track_name and raw_score, and
    optionally output_type, track_strand and biosample_name (the tidy format of AlphaGenome's tidy_scores and of
    trackpilot score-alphagenome). Tracks are keyed output_type|track_name|track_strand and sorted by key.
    A wide table (first column variant_id, one column per track) is read as is.
    """
    if isinstance(path_or_df, pd.DataFrame):
        p = path_or_df
    else:
        sep = "," if str(path_or_df).endswith((".csv", ".csv.gz")) else "\t"
        p = pd.read_csv(path_or_df, sep=sep, low_memory=False)
    if "raw_score" not in p:                                  # wide
        X = p.set_index(p.columns[0])
        meta = pd.DataFrame(dict(track_name=X.columns, biosample_name=X.columns, output_type=""), index=X.columns)
        return X.astype(float), meta
    if "output_type" in p and outputs:
        p = p[p.output_type.isin(outputs)]
    for c, v in [("output_type", ""), ("track_strand", "."), ("biosample_name", "")]:
        if c not in p:
            p = p.assign(**{c: v})
    p = p.assign(track=p.output_type.astype(str) + "|" + p.track_name.astype(str) + "|"
                 + p.track_strand.fillna(".").astype(str))
    X = p.pivot_table(index="variant_id", columns="track", values="raw_score", aggfunc="mean")
    cols = [c for c in ["output_type", "track_name", "track_strand", "biosample_name", "biosample_type"] if c in p]
    meta = p.drop_duplicates("track").set_index("track").loc[X.columns, cols]
    return X, meta


def read_labels(path, column=None):
    """Measured effects as a Series indexed by variant_id (column `effect` unless named otherwise)."""
    lab = pd.read_csv(path, sep=None, engine="python")
    col = column or ("effect" if "effect" in lab else lab.columns[1])
    return lab.set_index("variant_id")[col].astype(float)


def match_tracks(meta, keywords=(), biosamples=()):
    """Tracks whose biosample name or track name contains any keyword (case-insensitive), or whose biosample name
    equals one of `biosamples` exactly."""
    text = meta.biosample_name.fillna("").astype(str) + " " + meta.track_name.astype(str)
    hit = np.zeros(len(meta), bool)
    for k in keywords:
        if k:
            hit |= text.str.contains(k, case=False, regex=False).to_numpy()
    if len(biosamples):
        hit |= meta.biosample_name.isin([b for b in biosamples if b]).to_numpy()
    return hit


def split_list(s):
    return [x.strip() for x in (s or "").replace(",", ";").split(";") if x.strip()]
