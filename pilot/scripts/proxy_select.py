"""Pilot step 3 (A1): label-efficient choice of a proxy track among all AlphaGenome output tracks.

Question: an element is measured in a cell type the model may not profile. Among the ~1,000
CAGE/DNase/ATAC tracks AlphaGenome predicts, which one should stand in for that cell type, and how
few measured variants are needed to find it?

Model. Each candidate track t is a hypothesis "measured effect y = a + b * z_t + noise", where z_t is
the track's predicted variant effect divided by its standard deviation over the element (not centred,
so the sign keeps the alt-vs-ref direction). With Gaussian priors on a and b, the measured panel y_S
under track t is distributed as
    y_S ~ N(0, s^2 [(1 - beta) I + alpha 11' + beta z_tS z_tS'])        (s^2 = mean of y_S^2)
and the posterior over tracks is proportional to that marginal likelihood times the prior. The chosen
proxy is the MAP hypothesis; the posterior-weighted average of all predictive means is reported too.

Two optional refinements (off by default: on the development elements neither improved the
posterior-weighted prediction, see docs/internal/pilot_A1_results_zh.md):
  ensemble hypotheses  the all-track mean (and, when the element's own cell type is profiled, the mean of
                       its metadata-matched tracks) are added as hypotheses holding --ensemble-prior of the
                       prior mass, so with few or noisy labels the answer falls back to the average.
  measurement noise    each measured variant carries its MPRA standard error. In Kircher 2019 it is
                       c / sqrt(barcodes) with one constant c per element, so it is known before a variant
                       is measured. The covariance becomes s^2 [...] + diag(se^2) with s^2 the signal part
                       of mean(y_S^2), and the sequential panel divides track disagreement by the variant's
                       expected noise (residual misfit + se^2).

Panels (k measured variants, chosen from a random pool of variants that "could be assayed"):
  random        k variants at random
  high_effect   largest |all-track mean prediction| (the usual "pick strong predicted hits")
  disagree      largest unweighted variance across tracks (static committee disagreement)
  disc          sequential: each next variant maximises the posterior-weighted variance of the
                hypotheses' predictive means over its expected noise (Box-Hill style model discrimination)
  disc_stop     disc, stopping once the recommended proxy has stopped changing: the MAP track's
                predictions correlate >= --stop-r with those of the MAP tracks of the previous
                --stop-patience panels (uses no held-out labels)
All panels use the same track-choice rule, so they differ only in which variants get measured.
With --regression, few-label regression baselines (the AlphaGenome CAGI5 protocol fits LASSO over all tracks
with every label) are also fitted on the random and disc panels of the first --regression-repeats repeats: ridge (penalty by exact leave-one-out) and
LASSO (scikit-learn LassoCV, leave-one-out), both over all single tracks. They add rows only; every other
number is unchanged.
Label-free references: matched (mean of metadata-matched tracks), all_mean, oracle (best single raw track
using every label; a reference, not a method and not an upper bound: the fitted slope may be negative and a
posterior average can beat every single track).
With --common-eval, every strategy and reference is scored on the variants outside the repeat's pool.
Method v2 (docs/internal/method_v2_zh.md): --anchor-prior adds the mean of the matched tracks as a hypothesis holding
that prior mass, with its slope prior centred on +sqrt(beta)*s; --family-prior splits the tracks' prior equally between
prediction clusters (no labels used; duplicating a track changes nothing); --controls adds same-budget metadata
controls (+matched_cal, +meta_ridge) on the random and disc panels.
With --comparators (post hoc, docs/internal/acquisition_comparators_zh.md), on the first --comparator-repeats repeats:
  disc_static   one round: the k variants with the largest prior-weighted disagreement (disc's first score)
  batch_kb      one round: disc with each chosen variant's label replaced by the current posterior-average
                prediction (kriging believer), so the panel is fixed before any measurement
  eig           sequential: expected information gain about the hypothesis index, with each hypothesis'
                calibration uncertainty in its predictive variance
  disc_stop_bma@r  disc, stopping once the posterior-average predictions of all variants correlate >= r with
                those of the previous --stop-patience panels
Scenarios: all tracks, and hide_matched (tracks of the element's own cell type removed, which mimics
a cell type the model does not profile).

usage:
  python proxy_select.py <pilot_dir>                                  # every scored element
  python proxy_select.py <pilot_dir> --elements PKLR-48h --repeats 50
  python proxy_select.py <pilot_dir> --ensemble --noise --tag full          # ablation: both refinements
Output: <pilot_dir>/results/proxy_select[_<tag>]/{runs.tsv.gz, summary.tsv, tracks/<element>.tsv}
"""
import os, glob, argparse, time, warnings
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.linear_model import LassoCV
from sklearn.model_selection import LeaveOneOut

ap = argparse.ArgumentParser()
ap.add_argument("pilot_dir")
ap.add_argument("--elements", default="", help="comma-separated; default every element with a prediction file")
ap.add_argument("--pred-dir", default="predictions/alphagenome_all", help="relative to pilot_dir")
ap.add_argument("--variants", default="data/variants_grch38.tsv", help="relative to pilot_dir; the Siraj table is "
                "data/siraj/proxy_variants.tsv (with --pred-dir predictions/siraj_all)")
ap.add_argument("--version", default="ALL_FOLDS")
ap.add_argument("--outputs", default="CAGE,DNASE,ATAC")
ap.add_argument("--budgets", default="2,4,8,12,16")
ap.add_argument("--repeats", type=int, default=200)
ap.add_argument("--pool-frac", type=float, default=0.5, help="fraction of variants available for a panel")
ap.add_argument("--alpha", type=float, default=0.1, help="prior variance of the intercept (units of s^2)")
ap.add_argument("--beta", type=float, default=0.5, help="prior fraction of variance a track explains")
ap.add_argument("--ensemble-prior", type=float, default=0.5, help="prior mass shared by the ensemble hypotheses")
ap.add_argument("--ensemble", action="store_true", help="add ensemble hypotheses (ablation)")
ap.add_argument("--noise", action="store_true", help="use per-variant measurement error (ablation)")
ap.add_argument("--tag", default="", help="suffix of the output directory")
ap.add_argument("--stop-r", type=float, default=0.95)
ap.add_argument("--stop-patience", type=int, default=2)
ap.add_argument("--min-labels", type=int, default=3, help="disc_stop never stops before this many labels")
ap.add_argument("--cluster-r", type=float, default=0.9, help="tracks correlated above this share a cluster (reported only)")
ap.add_argument("--regression", action="store_true", help="also fit few-label ridge/LASSO baselines")
ap.add_argument("--regression-repeats", type=int, default=50, help="regression baselines on the first N repeats (slow)")
ap.add_argument("--common-eval", action="store_true",
                help="score every strategy and label-free reference on the same variants: those outside the "
                     "repeat's pool (added 2026-10-29 after review; the frozen v1 scores each panel on its own "
                     "non-panel variants and the references on all variants)")
ap.add_argument("--anchor-prior", type=float, default=0.0,
                help="method v2: prior mass of the metadata anchor hypothesis (mean of the matched tracks, positive "
                     "slope prior); 0 = no anchor (v1). See docs/internal/method_v2_zh.md")
ap.add_argument("--anchor-slope-sd", type=float, default=0.5,
                help="sd of the anchor's slope prior as a fraction of its mean sqrt(beta)*s")
ap.add_argument("--family-prior", action="store_true",
                help="method v2: equal prior mass per prediction cluster (--cluster-r), split within the cluster")
ap.add_argument("--controls", action="store_true",
                help="same-budget metadata controls where matched tracks exist: +matched_cal and +meta_ridge")
ap.add_argument("--comparators", action="store_true",
                help="post hoc acquisition comparators on the first --comparator-repeats repeats "
                     "(docs/internal/acquisition_comparators_zh.md): disc_static, batch_kb, eig and disc_stop_bma@r")
ap.add_argument("--comparator-repeats", type=int, default=50)
ap.add_argument("--eig-top", type=int, default=100, help="eig: hypotheses kept (largest posterior weight)")
ap.add_argument("--eig-grid", type=int, default=96, help="eig: quadrature points for the mixture entropy")
ap.add_argument("--stop-bma-grid", default="0.95,0.98,0.99,0.995",
                help="disc_stop_bma thresholds: stop once consecutive posterior-average predictions correlate >= r")
ap.add_argument("--scenarios", default="all,hide_exact,hide_matched")
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
assert not (a.ensemble and a.anchor_prior > 0), "--ensemble and --anchor-prior are alternatives"

budgets = [int(k) for k in a.budgets.split(",")]
stop_bma_grid = [float(r) for r in a.stop_bma_grid.split(",")] if a.comparators else []
kmax = max(budgets)
outputs = a.outputs.split(",")
pred_dir = os.path.join(a.pilot_dir, a.pred_dir)
out_dir = os.path.join(a.pilot_dir, "results", "proxy_select" + (f"_{a.tag}" if a.tag else ""))
os.makedirs(os.path.join(out_dir, "tracks"), exist_ok=True)

var = pd.read_csv(os.path.join(a.pilot_dir, a.variants), sep="\t", low_memory=False)
var["key"] = var.chrom + ":" + var.pos.astype(str) + ":" + var.ref + ">" + var.alt
if "pred_element" not in var:   # Kircher: one prediction file per experiment
    var["pred_element"] = var.element
scored = {os.path.basename(f).split("__")[0] for f in glob.glob(os.path.join(pred_dir, f"*__{a.version}.tsv.gz"))}
elements = (a.elements.split(",") if a.elements else
            sorted(var.element[var.pred_element.isin(scored)].unique()))


def load(element):
    """Measured effects, their standard errors, and the variant x track matrix of predictions, scaled
    per track."""
    lab = var[var.element == element]
    cols = ["variant_id", "output_type", "track_name", "track_strand", "biosample_name", "biosample_type", "raw_score"]
    p = pd.read_csv(os.path.join(pred_dir, f"{lab.pred_element.iloc[0]}__{a.version}.tsv.gz"), sep="\t", usecols=cols,
                    low_memory=False)
    p = p[p.output_type.isin(outputs)]
    p["track"] = p.output_type + "|" + p.track_name + "|" + p.track_strand.fillna(".").astype(str)
    X = p.pivot_table(index="variant_id", columns="track", values="raw_score", aggfunc="mean")
    lab = lab[lab.key.isin(X.index)]
    X = X.loc[lab.key]
    sd = X.std(ddof=0)
    X = X.loc[:, (sd > 0) & X.notna().all()]
    meta = p.drop_duplicates("track").set_index("track").loc[X.columns, ["output_type", "track_name", "track_strand",
                                                                          "biosample_name", "biosample_type"]]
    kws = [k for k in lab.model_biosample_keywords.iloc[0].split(";") if k]
    text = meta.biosample_name.fillna("").astype(str) + " " + meta.track_name.astype(str)
    meta["matched"] = np.any([text.str.contains(k, case=False, regex=False) for k in kws], axis=0)
    # related biosamples (exact names; second dataset only): hidden with the matched ones in hide_matched, and the
    # anchor of hide_exact, where only the matched cell line's own tracks are hidden
    rel = lab.related_biosamples.iloc[0] if "related_biosamples" in lab else ""
    meta["related"] = meta.biosample_name.isin([k for k in str(rel).split(";") if k]) & ~meta.matched
    # standard error c / sqrt(barcodes); c from variants whose p-value pins it down (p < 0.5, effect != 0)
    # (second dataset: the published standard errors)
    if "barcodes" in lab:
        ok = (lab.pvalue < 0.5) & (lab.effect != 0)
        c = np.median(lab.effect[ok].abs() / norm.isf(lab.pvalue[ok].clip(lower=1e-300) / 2) * np.sqrt(lab.barcodes[ok]))
        se = c / np.sqrt(lab.barcodes.to_numpy(float))
    else:
        se = lab.se.to_numpy(float)
    if not a.noise:
        se = np.zeros_like(se)
    return (lab.effect.to_numpy(float), se, lab.truth_class.to_numpy(int),
            (X / X.std(ddof=0)).to_numpy(float), meta)


def rank_corr(M, y):
    """Spearman correlation of each column of M (n, T) with y (n,)."""
    R = rankdata(M, axis=0) if M.ndim == 2 else rankdata(M)[:, None]
    r = rankdata(y)
    R = R - R.mean(0)
    r = r - r.mean()
    return (R * r[:, None]).sum(0) / np.sqrt((R ** 2).sum(0) * (r ** 2).sum())


def posterior(Zs, ys, ses, logprior, slope=None):
    """Log posterior (unnormalised) of every hypothesis for panel (Zs (k, T), ys (k,), standard errors
    ses (k,)), the posterior-mean intercept and slope of each hypothesis' map from scaled prediction to
    the true effect, and the signal variance s^2. slope = (mean (T,), variance (T,)) of each hypothesis'
    slope prior in units of s and s^2; default N(0, beta s^2) for all (v1)."""
    k, T = Zs.shape
    m2 = float(np.mean(ys ** 2))
    s2 = max(m2 - float(np.mean(ses ** 2)), 0.1 * m2, 1e-6)
    zt = Zs.T
    if slope is None:
        S = s2 * ((1 - a.beta) * np.eye(k) + a.alpha + a.beta * zt[:, :, None] * zt[:, None, :]) + np.diag(ses ** 2)
        L = np.linalg.cholesky(S)
        c = np.linalg.solve(S, np.broadcast_to(ys, (T, k))[..., None])[..., 0]
        ll = -0.5 * ((ys * c).sum(1) + 2 * np.log(np.diagonal(L, axis1=1, axis2=2)).sum(1)) + logprior
        return ll, s2 * a.alpha * c.sum(1), s2 * a.beta * (zt * c).sum(1), s2
    mb, vb = slope
    s = np.sqrt(s2)
    S = s2 * ((1 - a.beta) * np.eye(k) + a.alpha + vb[:, None, None] * zt[:, :, None] * zt[:, None, :]) + np.diag(ses ** 2)
    L = np.linalg.cholesky(S)
    r = ys - s * mb[:, None] * zt                       # (T, k) residual from each prior-mean line
    c = np.linalg.solve(S, r[..., None])[..., 0]
    ll = -0.5 * ((r * c).sum(1) + 2 * np.log(np.diagonal(L, axis1=1, axis2=2)).sum(1)) + logprior
    return ll, s2 * a.alpha * c.sum(1), s * mb + s2 * vb * (zt * c).sum(1), s2


def softmax(ll):
    w = np.exp(ll - ll.max())
    return w / w.sum()


def evaluate(Z, y, se, cls, panel, track_rho, logprior, test=None, slope=None):
    """Fit the posterior on the panel and score the MAP hypothesis and the posterior average on `test`
    (default: every variant outside the panel)."""
    ll, A, B, _ = posterior(Z[panel], y[panel], se[panel], logprior, slope)
    w = softmax(ll)
    if test is None:
        test = np.ones(len(y), bool)
        test[panel] = False
    t = int(np.argmax(ll))
    pick = A[t] + B[t] * Z[test, t]
    bma = (A + Z[test] * B) @ w
    sig = cls[test] != 0
    return dict(pick=t, w=w, rho_pick=rank_corr(pick, y[test])[0], rho_bma=rank_corr(bma, y[test])[0],
                acc_pick=np.mean(np.sign(pick[sig]) == cls[test][sig]),
                acc_bma=np.mean(np.sign(bma[sig]) == cls[test][sig]),
                # diagnostic (docs/internal/direction_accuracy_zh.md): the same average without the fitted offsets
                acc_bma0=np.mean(np.sign((Z[test] * B) @ w)[sig] == cls[test][sig]),
                pick_pct=np.mean(track_rho <= track_rho[t]))


def regress(Z, y, panel):
    """Few-label ridge and LASSO over all single tracks, fitted on the panel; predictions for every variant."""
    out = {"ridge": ridge(Z, y, panel)}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lasso = LassoCV(alphas=20, cv=LeaveOneOut(), max_iter=5000).fit(Z[panel], y[panel])
    out["lasso"] = lasso.predict(Z)
    return out


def ridge(Z, y, panel):
    """Ridge over all single tracks fitted on the panel, penalty by exact leave-one-out; predictions for
    every variant."""
    Xs, ys = Z[panel], y[panel]
    xm, ym = Xs.mean(0), ys.mean()
    Xc, yc = Xs - xm, ys - ym
    K = Xc @ Xc.T
    best = None
    for lam in np.trace(K) / len(panel) * np.logspace(-3, 3, 13):
        Hm = K @ np.linalg.inv(K + lam * np.eye(len(panel)))      # hat matrix of the centred fit
        loo = np.mean(((yc - Hm @ yc) / np.maximum(1 - np.diag(Hm) - 1 / len(panel), 1e-6)) ** 2)
        if best is None or loo < best[0]:
            best = (loo, lam)
    w = Xc.T @ np.linalg.solve(K + best[1] * np.eye(len(panel)), yc)
    return ym + (Z - xm) @ w


def disc_sequence(Z, y, se, pool, C, logprior, slope=None, fantasy=False):
    """Sequential discriminating panel. Returns the panel order, the first panel size at which the
    MAP track correlates >= --stop-r with the MAP tracks of the previous --stop-patience panels (kmax if
    never), and for each --stop-bma-grid threshold r the first panel size at which the posterior-average
    predictions of all variants correlate >= r with those of each of the previous --stop-patience panels.
    C is the track x track correlation of predictions. fantasy=True (batch_kb) never reads y: each chosen
    variant gets the current posterior-average prediction as its label ("kriging believer"), so the
    whole panel is fixed before any measurement."""
    w = softmax(logprior)
    mu = Z[pool]                       # before any label: each hypothesis' own prediction, positive slope
    noise = np.median(se[pool] ** 2) + se[pool] ** 2   # before any label, misfit assumed ~ typical noise
    left = np.ones(len(pool), bool)
    yy = np.zeros_like(y) if fantasy else y
    order, k_stop, maps, bmas = [], kmax, [], []
    k_bma = {r: kmax for r in stop_bma_grid}
    for step in range(kmax):
        mbar = mu @ w
        score = ((mu - mbar[:, None]) ** 2) @ w / np.maximum(noise, 1e-12)
        score[~left] = -np.inf
        j = int(np.argmax(score))
        left[j] = False
        order.append(pool[j])
        if fantasy:
            yy[pool[j]] = mbar[j]
        ll, A, B, s2 = posterior(Z[order], yy[order], se[order], logprior, slope)
        w = softmax(ll)
        mu = A + Z[pool] * B
        noise = s2 * (1 - a.beta) + se[pool] ** 2
        if fantasy:
            continue
        maps.append(int(np.argmax(ll)))
        prev = maps[-1 - a.stop_patience:-1]
        if (k_stop == kmax and len(order) >= a.min_labels and len(prev) == a.stop_patience
                and all(C[maps[-1], t] >= a.stop_r for t in prev)):
            k_stop = len(order)
        if not k_bma:
            continue
        bmas.append(A @ w + Z @ (B * w))
        prev = bmas[-1 - a.stop_patience:-1]
        if len(order) >= a.min_labels and len(prev) == a.stop_patience:
            r_min = min(np.corrcoef(bmas[-1], p)[0, 1] for p in prev)
            for r in k_bma:
                if k_bma[r] == kmax and r_min >= r:
                    k_bma[r] = len(order)
    return order, k_stop, k_bma


def eig_sequence(Z, y, se, pool, logprior, slope=None):
    """Sequential panel by expected information gain about the hypothesis index (the quantity Box-Hill's
    disagreement approximates). Under hypothesis t the next label at x is Gaussian with the posterior-mean
    line A_t + B_t z_t(x) and variance noise + [1, z_t(x)] Cov_t(a, b) [1, z_t(x)]'; the information gain
    is the entropy of the weighted mixture minus the mean entropy of its components, by quadrature on
    --eig-grid points, over the --eig-top hypotheses with the largest weight. Before any label s = 1
    (the gain does not depend on the unit when every line passes through 0)."""
    T = Z.shape[1]
    mb, vb = (np.zeros(T), np.full(T, a.beta)) if slope is None else slope
    w, A, B, s2 = softmax(logprior), np.zeros(T), mb.copy(), 1.0
    left = np.ones(len(pool), bool)
    order = []
    for step in range(kmax):
        sig2 = s2 * (1 - a.beta)
        p11, p22 = np.full(T, 1 / (a.alpha * s2)), 1 / (vb * s2)        # posterior precision of (a, b)
        p12 = np.zeros(T)
        if order:
            d = 1 / (sig2 + se[order] ** 2)
            zo = Z[order]
            p11, p12, p22 = p11 + d.sum(), d @ zo, p22 + d @ zo ** 2
        det = p11 * p22 - p12 ** 2
        top = np.argsort(-w)[:a.eig_top]
        wt = w[top] / w[top].sum()
        zc = Z[pool][:, top]
        m = A[top] + zc * B[top]
        v = sig2 + se[pool, None] ** 2 + (p22[top] - 2 * p12[top] * zc + p11[top] * zc ** 2) / det[top]
        sd = np.sqrt(v)
        gain = np.full(len(pool), -np.inf)
        for c in range(0, len(pool), 64):
            mc, sc, vc = m[c:c + 64], sd[c:c + 64], v[c:c + 64]
            lo, hi = (mc - 6 * sc).min(1), (mc + 6 * sc).max(1)
            g = lo[:, None] + (hi - lo)[:, None] * np.linspace(0, 1, a.eig_grid)[None, :]
            dens = (np.exp(-0.5 * (g[:, :, None] - mc[:, None, :]) ** 2 / vc[:, None, :]) / sc[:, None, :]) @ wt
            dens /= np.sqrt(2 * np.pi)
            h_mix = -(dens * np.log(np.maximum(dens, 1e-300))).sum(1) * (hi - lo) / (a.eig_grid - 1)
            gain[c:c + 64] = h_mix - 0.5 * np.log(2 * np.pi * np.e * vc) @ wt
        gain[~left] = -np.inf
        j = int(np.argmax(gain))
        left[j] = False
        order.append(pool[j])
        ll, A, B, s2 = posterior(Z[order], y[order], se[order], logprior, slope)
        w = softmax(ll)
    return order


runs, t0 = [], time.time()
for element in elements:
    y, se, cls, Z_all, meta_all = load(element)
    n = len(y)
    scenarios = {"all": np.ones(len(meta_all), bool)}
    if meta_all.matched.any():
        scenarios["hide_matched"] = ~(meta_all.matched | meta_all.related).to_numpy()
    if meta_all.matched.any() and meta_all.related.any():
        scenarios["hide_exact"] = ~meta_all.matched.to_numpy()
    scenarios = {k: v for k, v in scenarios.items() if k in a.scenarios.split(",")}
    track_rho_all = rank_corr(Z_all, y)
    clusters_all = fcluster(linkage(Z_all.T, method="average", metric="correlation"),
                            t=1 - a.cluster_r, criterion="distance")
    sig = cls != 0
    tracks = meta_all.assign(spearman=track_rho_all, dir_acc_sig=(np.sign(Z_all[sig]) == cls[sig, None]).mean(0),
                             cluster=clusters_all)
    tracks["cluster_size"] = tracks.groupby("cluster").cluster.transform("size")
    tracks.sort_values("spearman", ascending=False).to_csv(os.path.join(out_dir, "tracks", f"{element}.tsv"), sep="\t")
    print(f"{element}: {n} SNVs ({sig.sum()} significant), {Z_all.shape[1]} tracks, "
          f"{meta_all.matched.sum()} matched, {len(np.unique(clusters_all))} clusters", flush=True)

    for scen, keep in scenarios.items():
        Z, meta, track_rho = Z_all[:, keep], meta_all[keep], track_rho_all[keep]
        clusters = np.unique(clusters_all[keep], return_inverse=True)[1]
        base = dict(element=element, scenario=scen, n_tracks=int(keep.sum()))
        best = int(np.argmax(track_rho))
        rho_oracle = track_rho[best]
        # hypotheses: every single track, plus the ensemble averages (scaled like the tracks)
        ens = {"ensemble|all_mean": Z.mean(1)} if a.ensemble else {}
        if ens and scen == "all" and meta.matched.any():
            ens["ensemble|matched"] = Z[:, meta.matched.to_numpy()].mean(1)
        # metadata anchor: the matched tracks when kept, the related ones when only the matched are hidden
        anchor_mask = (meta.matched if scen == "all" else meta.related if scen == "hide_exact"
                       else pd.Series(False, index=meta.index)).to_numpy()
        has_matched = bool(anchor_mask.any())
        if a.anchor_prior > 0 and has_matched:
            ens = {"anchor|matched": Z[:, anchor_mask].mean(1)}
        T = Z.shape[1]
        H = np.column_stack([Z] + [v / v.std() for v in ens.values()])
        hmeta = pd.concat([meta, pd.DataFrame(dict(biosample_name=list(ens), matched=[k.endswith("matched") for k in ens]),
                                              index=list(ens))])
        h_rho = np.concatenate([track_rho, rank_corr(H[:, T:], y)]) if ens else track_rho
        h_clusters = np.concatenate([clusters, clusters.max() + 1 + np.arange(len(ens))])
        logprior = np.full(H.shape[1], np.log((1 - a.ensemble_prior) / T) if ens else -np.log(T))
        if ens:
            logprior[T:] = np.log(a.ensemble_prior / len(ens))
        slope = None
        anchor_slope = (np.array([np.sqrt(a.beta)]), np.array([(a.anchor_slope_sd * np.sqrt(a.beta)) ** 2]))
        if a.family_prior or a.anchor_prior > 0:
            # method v2: the tracks' mass is split equally between prediction clusters, then within each;
            # the anchor holds --anchor-prior and its slope prior is centred on a positive value
            lp = (-np.log(clusters.max() + 1) - np.log(np.bincount(clusters)[clusters]) if a.family_prior
                  else np.full(T, -np.log(T)))
            mass = a.anchor_prior if ens else 0.0
            logprior = np.concatenate([lp + np.log1p(-mass), np.full(len(ens), np.log(max(mass, 1e-300)))])
            if ens:
                slope = (np.concatenate([np.zeros(T), anchor_slope[0]]), np.concatenate([np.full(T, a.beta), anchor_slope[1]]))
        C = np.corrcoef(H.T)

        def label_free(name, pred, m=None, rep=-1):
            m = np.ones(n, bool) if m is None else m
            runs.append(dict(base, strategy=name, k=0, repeat=rep, labels=0, rho_pick=rank_corr(pred[m], y[m])[0],
                             acc_pick=np.mean(np.sign(pred[m & sig]) == cls[m & sig]), rho_oracle=rho_oracle))
        refs = {"oracle": Z[:, best], "all_mean": Z.mean(1)}
        if has_matched:
            refs["matched" if scen == "all" else "related"] = Z[:, anchor_mask].mean(1)
        if not a.common_eval:
            for name, pred in refs.items():
                label_free(name, pred)

        mean_pred = np.abs(Z.mean(1))
        spread = Z.var(1)
        for rep in range(a.repeats):
            rng = np.random.default_rng([a.seed, rep])
            pool = rng.choice(n, int(a.pool_frac * n), replace=False)
            outside = None
            if a.common_eval:
                outside = np.ones(n, bool)
                outside[pool] = False
                for name, pred in refs.items():
                    label_free(name, pred, outside, rep)
            order, k_stop, k_bma = disc_sequence(H, y, se, pool, C, logprior, slope)
            panels = {"random": rng.permutation(pool),
                      "high_effect": pool[np.argsort(-mean_pred[pool], kind="stable")],
                      "disagree": pool[np.argsort(-spread[pool], kind="stable")],
                      "disc": np.array(order)}
            todo = [(s, k, p[:k]) for s, p in panels.items() for k in budgets]
            todo.append(("disc_stop", kmax, np.array(order[:k_stop])))
            if a.comparators and rep < a.comparator_repeats:
                # one-round panels (no label seen while choosing) and the information-gain panel; the same
                # inference scores all of them (no random draws here, so every other row is unchanged)
                w0, mu0 = softmax(logprior), H[pool]
                static = ((mu0 - (mu0 @ w0)[:, None]) ** 2) @ w0
                extra = {"disc_static": pool[np.argsort(-static, kind="stable")],
                         "batch_kb": np.array(disc_sequence(H, y, se, pool, C, logprior, slope, fantasy=True)[0]),
                         "eig": np.array(eig_sequence(H, y, se, pool, logprior, slope))}
                todo += [(s, k, p[:k]) for s, p in extra.items() for k in budgets]
                todo += [(f"disc_stop_bma@{r:g}", kmax, np.array(order[:k])) for r, k in k_bma.items()]
            for strategy, k, panel in todo:
                r = evaluate(H, y, se, cls, panel, h_rho, logprior, outside, slope)
                t = r["pick"]
                runs.append(dict(base, strategy=strategy, k=k, repeat=rep, labels=len(panel),
                                 pick_track=hmeta.index[t], pick_biosample=hmeta.biosample_name.iloc[t],
                                 pick_matched=bool(hmeta.matched.iloc[t]), pick_ensemble=t >= T,
                                 ensemble_mass=float(r["w"][T:].sum()), pick_cluster_mass=float(
                                     r["w"][h_clusters == h_clusters[t]].sum()), exp_loss=float(1 - r["w"] @ C[:, t]),
                                 rho_pick=r["rho_pick"], rho_bma=r["rho_bma"], acc_pick=r["acc_pick"],
                                 acc_bma=r["acc_bma"], acc_bma0=r["acc_bma0"], pick_pct=r["pick_pct"], rho_oracle=rho_oracle,
                                 regret=rho_oracle - h_rho[t]))
            if a.regression and rep < a.regression_repeats:
                for pname in ["random", "disc"]:
                    for k in [4, 8, 16]:
                        panel = panels[pname][:k]
                        if outside is None:
                            test = np.ones(n, bool)
                            test[panel] = False
                        else:
                            test = outside
                        sig_t = cls[test] != 0
                        for name, pred in regress(Z, y, panel).items():
                            runs.append(dict(base, strategy=f"{pname}+{name}", k=k, repeat=rep, labels=k,
                                             rho_bma=rank_corr(pred[test], y[test])[0],
                                             acc_bma=np.mean(np.sign(pred[test][sig_t]) == cls[test][sig_t]),
                                             rho_oracle=rho_oracle))
            if a.controls and has_matched:
                # same-budget metadata controls: the matched mean recalibrated with a symmetric slope prior
                # (v1 inference on that one hypothesis), and metadata-centred ridge (matched mean calibrated
                # with the anchor's positive slope prior, then ridge over all tracks on the residuals)
                m = Z[:, anchor_mask].mean(1)
                m = (m / m.std())[:, None]
                for pname in ["random", "disc"]:
                    for k in [4, 8, 16]:
                        panel = panels[pname][:k]
                        if outside is None:
                            test = np.ones(n, bool)
                            test[panel] = False
                        else:
                            test = outside
                        sig_t = cls[test] != 0
                        _, A, B, _ = posterior(m[panel], y[panel], se[panel], np.zeros(1))
                        cal = A[0] + B[0] * m[:, 0]
                        _, A, B, _ = posterior(m[panel], y[panel], se[panel], np.zeros(1), anchor_slope)
                        off = A[0] + B[0] * m[:, 0]
                        for name, pred in {"matched_cal": cal, "meta_ridge": off + ridge(Z, y - off, panel)}.items():
                            runs.append(dict(base, strategy=f"{pname}+{name}", k=k, repeat=rep, labels=k,
                                             rho_bma=rank_corr(pred[test], y[test])[0],
                                             acc_bma=np.mean(np.sign(pred[test][sig_t]) == cls[test][sig_t]),
                                             rho_oracle=rho_oracle))
        print(f"  {scen}: {a.repeats} repeats done, {time.time() - t0:.0f}s", flush=True)

runs = pd.DataFrame(runs)
for col in ["pick_matched", "pick_ensemble"]:
    runs[col] = runs[col].astype(float)  # NaN for the label-free references
runs.to_csv(os.path.join(out_dir, "runs.tsv.gz"), sep="\t", index=False)
g = runs.groupby(["element", "scenario", "strategy", "k"], sort=False)
summary = g.agg(repeats=("rho_pick", "size"), labels=("labels", "mean"), rho_pick=("rho_pick", "mean"),
                rho_pick_q05=("rho_pick", lambda s: s.quantile(0.05)), rho_bma=("rho_bma", "mean"),
                regret=("regret", "mean"), near_oracle=("regret", lambda s: (s <= 0.05).mean()),
                acc_pick=("acc_pick", "mean"), pick_matched=("pick_matched", "mean"),
                pick_ensemble=("pick_ensemble", "mean"), ensemble_mass=("ensemble_mass", "mean"),
                rho_oracle=("rho_oracle", "first")).reset_index()
summary.to_csv(os.path.join(out_dir, "summary.tsv"), sep="\t", index=False, float_format="%.4f")

# console view: label-free references, then every panel strategy at k = 8 and the stopping rule
show = summary[(summary.k == 0) | (summary.k == 8) | (summary.strategy == "disc_stop")]
cols = ["element", "scenario", "strategy", "k", "labels", "rho_pick", "rho_pick_q05", "rho_bma", "regret",
        "near_oracle", "acc_pick", "pick_ensemble"]
with pd.option_context("display.width", 200, "display.max_rows", 500):
    print(show[cols].round(3).to_string(index=False))
print(f"wrote {out_dir}")
