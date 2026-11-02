"""TrackPilot inference: tracks as calibrated hypotheses, a metadata prior, and the pilot panel.

Every candidate track t is a hypothesis "measured effect y = a + b * z_t + noise", where z_t is the track's predicted
variant effect divided by its standard deviation over all variants of the region (not centred, so the sign keeps the
alt-vs-ref direction). With Gaussian priors on the intercept a and slope b, both integrate out and the measured panel
y_S under track t is distributed as

    y_S ~ N(s * m_t z_tS, s^2 [(1 - beta) I + alpha 11' + v_t z_tS z_tS'])        s^2 = mean of y_S^2

with slope prior N(m_t s, v_t s^2): N(0, beta s^2) for single tracks and N(sqrt(beta) s, (0.5 sqrt(beta) s)^2) for
the metadata anchor (the mean of the tracks whose metadata name the cell type). The posterior over hypotheses is that
marginal likelihood times the prior; predictions are the posterior average of every hypothesis' calibrated line.

The pilot panel is chosen by Box-Hill model discrimination: each next variant maximises the posterior-weighted
variance of the hypotheses' predictive means. The one-round panel replaces each chosen variant's label by the current
posterior-average prediction, so the whole panel is fixed before any measurement.

This is the frozen method of the paper (alpha = 0.1, beta = 0.5, anchor prior 0.95, family prior at r = 0.9),
ported from pilot/scripts/proxy_select.py without changing any computation; tests/test_reproduce.py checks it
against the pre-registered runs.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage

ANCHOR = "anchor"


def softmax(ll):
    w = np.exp(ll - ll.max())
    return w / w.sum()


def track_families(Z, r=0.9):
    """Cluster labels (0-based) of the columns of Z: average linkage on correlation distance, cut at 1 - r."""
    cl = fcluster(linkage(Z.T, method="average", metric="correlation"), t=1 - r, criterion="distance")
    return np.unique(cl, return_inverse=True)[1]


@dataclass
class Fit:
    """Posterior after a set of measurements."""
    hypotheses: pd.Index
    weights: np.ndarray          # posterior probability of each hypothesis
    intercept: np.ndarray        # posterior-mean intercept of each hypothesis' calibration line
    slope: np.ndarray            # posterior-mean slope (per unit of the scaled prediction)
    log_posterior: np.ndarray    # unnormalised
    s2: float                    # signal variance used to set the prior scale
    _Z: np.ndarray
    _variants: pd.Index

    def predict(self):
        """Posterior-average prediction for every variant."""
        return pd.Series((self.intercept + self._Z * self.slope) @ self.weights, index=self._variants,
                         name="prediction")

    def predict_map(self):
        """Prediction of the single most probable hypothesis, calibrated."""
        t = int(np.argmax(self.log_posterior))
        return pd.Series(self.intercept[t] + self.slope[t] * self._Z[:, t], index=self._variants, name="prediction_map")

    @property
    def map_hypothesis(self):
        return self.hypotheses[int(np.argmax(self.log_posterior))]


class TrackPilot:
    """Posterior over the tracks of one region.

    predictions  DataFrame, variants x tracks, raw predicted effects (e.g. from trackpilot.io.read_predictions)
    anchor       boolean mask over the tracks whose metadata name the cell type (or None); their mean becomes the
                 anchor hypothesis
    exclude      boolean mask of tracks to leave out (to simulate an unprofiled cell type, as in the paper); track
                 families are still computed over all tracks
    """

    def __init__(self, predictions: pd.DataFrame, anchor=None, exclude=None, *, anchor_prior=0.95, family_r=0.9,
                 family_prior=True, alpha=0.1, beta=0.5, anchor_slope_sd=0.5):
        X = predictions
        sd = X.std(ddof=0)
        ok = ((sd > 0) & X.notna().all()).to_numpy()
        self.dropped = list(X.columns[~ok])
        Z_all = (X.loc[:, ok] / X.loc[:, ok].std(ddof=0)).to_numpy(float)
        keep = np.ones(Z_all.shape[1], bool) if exclude is None else ~np.asarray(exclude, bool)[ok]
        anc = np.zeros(Z_all.shape[1], bool) if anchor is None else np.asarray(anchor, bool)[ok]
        anc = anc[keep]
        self.alpha, self.beta = alpha, beta
        self.variants = X.index
        self.tracks = X.columns[ok][keep]
        Z = Z_all[:, keep]
        clusters = np.unique(track_families(Z_all, family_r)[keep], return_inverse=True)[1]
        T = Z.shape[1]
        cols, names = [Z], list(self.tracks)
        mass = 0.0
        if anchor_prior > 0 and anc.any():
            m = Z[:, anc].mean(1)
            cols.append((m / m.std())[:, None])
            names.append(ANCHOR)
            mass = anchor_prior
        self.Z = np.column_stack(cols)
        self.hypotheses = pd.Index(names)
        self.families = np.concatenate([clusters, clusters.max() + 1 + np.arange(len(names) - T)])
        self.n_anchor_tracks = int(anc.sum()) if mass else 0
        lp = (-np.log(clusters.max() + 1) - np.log(np.bincount(clusters)[clusters]) if family_prior
              else np.full(T, -np.log(T)))
        n_extra = len(names) - T
        self.log_prior = np.concatenate([lp + np.log1p(-mass), np.full(n_extra, np.log(max(mass, 1e-300)))])
        self.slope_mean = np.concatenate([np.zeros(T), np.full(n_extra, np.sqrt(beta))])
        self.slope_var = np.concatenate([np.full(T, beta), np.full(n_extra, (anchor_slope_sd * np.sqrt(beta)) ** 2)])

    # ---- inference -------------------------------------------------------------------------------------------
    def _posterior(self, rows, ys):
        """Log posterior, posterior-mean intercepts and slopes of every hypothesis, and s^2, for labels ys at rows."""
        k = len(rows)
        zt = self.Z[rows].T
        mb, vb = self.slope_mean, self.slope_var
        s2 = max(float(np.mean(ys ** 2)), 1e-6)
        s = np.sqrt(s2)
        S = s2 * ((1 - self.beta) * np.eye(k) + self.alpha + vb[:, None, None] * zt[:, :, None] * zt[:, None, :])
        L = np.linalg.cholesky(S)
        r = ys - s * mb[:, None] * zt
        c = np.linalg.solve(S, r[..., None])[..., 0]
        ll = -0.5 * ((r * c).sum(1) + 2 * np.log(np.diagonal(L, axis1=1, axis2=2)).sum(1)) + self.log_prior
        return ll, s2 * self.alpha * c.sum(1), s * mb + s2 * vb * (zt * c).sum(1), s2

    def _rows(self, ids):
        pos = self.variants.get_indexer(pd.Index(ids))
        if (pos < 0).any():
            raise KeyError(f"{int((pos < 0).sum())} variants are not in the predictions")
        return pos

    def fit(self, labels: pd.Series) -> Fit:
        """Posterior given measured effects (a Series indexed by variant id)."""
        labels = labels.dropna()
        rows = self._rows(labels.index)
        ll, A, B, s2 = self._posterior(rows, labels.to_numpy(float))
        return Fit(self.hypotheses, softmax(ll), A, B, ll, s2, self.Z, self.variants)

    def prior_prediction(self):
        """Label-free default: the anchor (mean of the matched tracks) if there is one, else the all-track mean."""
        T = len(self.tracks)
        z = self.Z[:, T] if len(self.hypotheses) > T else self.Z[:, :T].mean(1)
        return pd.Series(z, index=self.variants, name="prior_prediction")

    # ---- choosing variants -----------------------------------------------------------------------------------
    def suggest(self, k, pool=None, measured: pd.Series | None = None, one_round=True, oracle: pd.Series | None = None):
        """Pilot panel of k variants from `pool` (default: every variant), after the already `measured` ones.

        one_round=True   the whole panel is fixed now: each chosen variant's label is replaced by the current
                         posterior-average prediction (what a laboratory ordering one batch does)
        oracle           a Series of true labels: choose sequentially, measuring each variant before choosing the
                         next (for benchmarks; one_round is then ignored)
        """
        pool_ids = self.variants if pool is None else pd.Index(pool)
        P = self._rows(pool_ids)
        measured = pd.Series(dtype=float) if measured is None else measured.dropna()
        order = list(self._rows(measured.index))
        yy = np.zeros(len(self.variants))
        yy[order] = measured.to_numpy(float)
        left = ~np.isin(P, order)
        chosen = []
        if order:
            ll, A, B, s2 = self._posterior(np.array(order), yy[order])
            w, mu, noise = softmax(ll), A + self.Z[P] * B, s2 * (1 - self.beta)
        else:
            w, mu, noise = softmax(self.log_prior), self.Z[P], 0.0
        for _ in range(min(k, int(left.sum()))):
            mbar = mu @ w
            score = ((mu - mbar[:, None]) ** 2) @ w / max(noise, 1e-12)
            score[~left] = -np.inf
            j = int(np.argmax(score))
            left[j] = False
            order.append(P[j])
            chosen.append(P[j])
            if oracle is not None:
                yy[P[j]] = oracle.loc[self.variants[P[j]]]
            else:
                yy[P[j]] = mbar[j] if one_round else np.nan
                if not one_round:
                    break           # sequential use without an oracle: return one variant, measure it, call again
            ll, A, B, s2 = self._posterior(np.array(order), yy[order])
            w, mu, noise = softmax(ll), A + self.Z[P] * B, s2 * (1 - self.beta)
        return list(self.variants[chosen])
