"""Unit tests on the synthetic example (examples/make_example.py) and on properties of the method."""
import os

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from trackpilot import TrackPilot, match_tracks, read_labels, read_predictions
from trackpilot.cli import main

EX = os.path.join(os.path.dirname(__file__), "..", "examples")
PRED = os.path.join(EX, "example_predictions.tsv.gz")
MEAS = os.path.join(EX, "example_measurements.tsv")


@pytest.fixture(scope="module")
def example():
    X, meta = read_predictions(PRED)
    return X, meta, read_labels(MEAS).loc[X.index]


def held_out(tp, y, panel):
    out = ~y.index.isin(panel)
    return spearmanr(tp.fit(y.loc[panel]).predict()[out], y[out])[0]


def test_read(example):
    X, meta, y = example
    assert X.shape == (600, 120) and len(meta) == 120 and y.notna().all()
    assert match_tracks(meta, ["lung"]).sum() == 10
    assert match_tracks(meta, ["fibroblast"]).sum() == 20
    assert match_tracks(meta, biosamples=["K562"]).sum() == 10


def test_unprofiled_recovers_the_right_family(example):
    X, meta, y = example
    tp = TrackPilot(X)
    assert len(np.unique(tp.families)) == 12
    panel = tp.suggest(8)
    fit = tp.fit(y.loc[panel])
    lung = match_tracks(meta.loc[tp.tracks], ["lung"])
    assert fit.weights[lung].sum() > 0.5
    assert held_out(tp, y, panel) > spearmanr(tp.prior_prediction(), y)[0] + 0.2


def test_first_variant_maximises_prior_disagreement(example):
    # (whether pilot panels beat random ones is an empirical question, answered on the paper's data in
    # test_reproduce.py; on this clean synthetic example random panels identify the family about as well)
    X, _, _ = example
    tp = TrackPilot(X)
    w = np.exp(tp.log_prior)
    var = ((tp.Z - (tp.Z @ w)[:, None]) ** 2) @ w
    assert tp.suggest(1) == [tp.variants[int(np.argmax(var))]]


def test_one_round_panel_reads_no_labels(example):
    X, _, y = example
    tp = TrackPilot(X)
    assert tp.suggest(8) == tp.suggest(8)
    y2 = y.sample(frac=1, random_state=0)
    assert tp.suggest(8, oracle=y) == tp.suggest(8, oracle=y2.loc[y.index])


def test_anchor(example):
    X, meta, y = example
    tp = TrackPilot(X, match_tracks(meta, ["lung"]))
    assert tp.hypotheses[-1] == "anchor" and tp.n_anchor_tracks == 10
    assert np.exp(tp.log_prior[-1]) == pytest.approx(0.95)
    assert np.exp(tp.log_prior).sum() == pytest.approx(1.0)
    # a wrong anchor loses its weight once the measurements contradict it
    wrong = TrackPilot(X, match_tracks(meta, ["keratinocyte"]))
    fit = wrong.fit(y.loc[wrong.suggest(8)])
    assert fit.weights[-1] < 0.5


def test_family_prior_ignores_duplicated_tracks(example):
    X, _, y = example
    dup = pd.concat([X, X.iloc[:, :10].add_suffix(" copy")], axis=1)   # ten more tracks in family 0
    w0 = np.exp(TrackPilot(X).log_prior)
    w1 = np.exp(TrackPilot(dup).log_prior)
    fam = TrackPilot(dup).families
    assert w1[fam == fam[0]].sum() == pytest.approx(w0[:10].sum())


def test_scale_invariance(example):
    X, _, y = example
    tp = TrackPilot(X)
    tp2 = TrackPilot(X * np.linspace(0.1, 10, X.shape[1]))
    panel = tp.suggest(8)
    assert panel == tp2.suggest(8)
    assert np.allclose(tp.fit(y.loc[panel]).predict(), tp2.fit(y.loc[panel]).predict())


def test_cli(tmp_path):
    panel = tmp_path / "panel.tsv"
    main(["suggest", PRED, "-k", "8", "--anchor", "lung", "-o", str(panel)])
    ids = pd.read_csv(panel, sep="\t").variant_id
    assert len(ids) == 8
    meas = read_labels(MEAS)
    lab = tmp_path / "measured.tsv"
    meas.loc[ids].reset_index().to_csv(lab, sep="\t", index=False)
    out, post, rep = tmp_path / "pred.tsv", tmp_path / "post.tsv", tmp_path / "report.html"
    main(["predict", PRED, str(lab), "--anchor", "lung", "-o", str(out), "--posterior", str(post), "--report", str(rep)])
    p = pd.read_csv(out, sep="\t")
    assert len(p) == 600 and p.measured.notna().sum() == 8
    assert pd.read_csv(post, sep="\t").weight.sum() == pytest.approx(1.0, abs=1e-4)
    assert "TrackPilot report" in rep.read_text()
    # suggest more after a first round
    more = tmp_path / "more.tsv"
    main(["suggest", PRED, "-k", "4", "--measured", str(lab), "-o", str(more)])
    assert not set(pd.read_csv(more, sep="\t").variant_id) & set(ids)
