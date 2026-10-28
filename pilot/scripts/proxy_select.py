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
Label-free references: matched (mean of metadata-matched tracks), all_mean, oracle (best single track
using every label; an upper bound, not a method).
Scenarios: all tracks, and hide_matched (tracks of the element's own cell type removed, which mimics
a cell type the model does not profile).

usage:
  python proxy_select.py <pilot_dir>                                  # every scored element
  python proxy_select.py <pilot_dir> --elements PKLR-48h --repeats 50
  python proxy_select.py <pilot_dir> --ensemble --noise --tag full          # ablation: both refinements
Output: <pilot_dir>/results/proxy_select[_<tag>]/{runs.tsv.gz, summary.tsv, tracks/<element>.tsv}
"""
import os, glob, argparse, time
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm
from scipy.cluster.hierarchy import linkage, fcluster

ap = argparse.ArgumentParser()
ap.add_argument("pilot_dir")
ap.add_argument("--elements", default="", help="comma-separated; default every element with a prediction file")
ap.add_argument("--pred-dir", default="predictions/alphagenome_all", help="relative to pilot_dir")
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
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()

budgets = [int(k) for k in a.budgets.split(",")]
kmax = max(budgets)
outputs = a.outputs.split(",")
pred_dir = os.path.join(a.pilot_dir, a.pred_dir)
out_dir = os.path.join(a.pilot_dir, "results", "proxy_select" + (f"_{a.tag}" if a.tag else ""))
os.makedirs(os.path.join(out_dir, "tracks"), exist_ok=True)

var = pd.read_csv(os.path.join(a.pilot_dir, "data", "variants_grch38.tsv"), sep="\t", low_memory=False)
var["key"] = var.chrom + ":" + var.pos.astype(str) + ":" + var.ref + ">" + var.alt
elements = (a.elements.split(",") if a.elements else
            sorted(os.path.basename(f).split("__")[0] for f in glob.glob(os.path.join(pred_dir, f"*__{a.version}.tsv.gz"))))


def load(element):
    """Measured effects, their standard errors, and the variant x track matrix of predictions, scaled
    per track."""
    lab = var[var.element == element]
    cols = ["variant_id", "output_type", "track_name", "track_strand", "biosample_name", "biosample_type", "raw_score"]
    p = pd.read_csv(os.path.join(pred_dir, f"{element}__{a.version}.tsv.gz"), sep="\t", usecols=cols, low_memory=False)
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
    # standard error c / sqrt(barcodes); c from variants whose p-value pins it down (p < 0.5, effect != 0)
    ok = (lab.pvalue < 0.5) & (lab.effect != 0)
    c = np.median(lab.effect[ok].abs() / norm.isf(lab.pvalue[ok].clip(lower=1e-300) / 2) * np.sqrt(lab.barcodes[ok]))
    se = c / np.sqrt(lab.barcodes.to_numpy(float))
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


def posterior(Zs, ys, ses, logprior):
    """Log posterior (unnormalised) of every hypothesis for panel (Zs (k, T), ys (k,), standard errors
    ses (k,)), the posterior-mean intercept and slope of each hypothesis' map from scaled prediction to
    the true effect, and the signal variance s^2."""
    k, T = Zs.shape
    m2 = float(np.mean(ys ** 2))
    s2 = max(m2 - float(np.mean(ses ** 2)), 0.1 * m2, 1e-6)
    zt = Zs.T
    S = s2 * ((1 - a.beta) * np.eye(k) + a.alpha + a.beta * zt[:, :, None] * zt[:, None, :]) + np.diag(ses ** 2)
    L = np.linalg.cholesky(S)
    c = np.linalg.solve(S, np.broadcast_to(ys, (T, k))[..., None])[..., 0]
    ll = -0.5 * ((ys * c).sum(1) + 2 * np.log(np.diagonal(L, axis1=1, axis2=2)).sum(1)) + logprior
    return ll, s2 * a.alpha * c.sum(1), s2 * a.beta * (zt * c).sum(1), s2


def softmax(ll):
    w = np.exp(ll - ll.max())
    return w / w.sum()


def evaluate(Z, y, se, cls, panel, track_rho, logprior):
    """Fit the posterior on the panel and score the MAP hypothesis and the posterior average on every
    variant outside the panel."""
    ll, A, B, _ = posterior(Z[panel], y[panel], se[panel], logprior)
    w = softmax(ll)
    test = np.ones(len(y), bool)
    test[panel] = False
    t = int(np.argmax(ll))
    pick = A[t] + B[t] * Z[test, t]
    bma = (A + Z[test] * B) @ w
    sig = cls[test] != 0
    return dict(pick=t, w=w, rho_pick=rank_corr(pick, y[test])[0], rho_bma=rank_corr(bma, y[test])[0],
                acc_pick=np.mean(np.sign(pick[sig]) == cls[test][sig]),
                acc_bma=np.mean(np.sign(bma[sig]) == cls[test][sig]),
                pick_pct=np.mean(track_rho <= track_rho[t]))


def disc_sequence(Z, y, se, pool, C, logprior):
    """Sequential discriminating panel. Returns the panel order and the first panel size at which the
    MAP track correlates >= --stop-r with the MAP tracks of the previous --stop-patience panels (kmax if
    never). C is the track x track correlation of predictions."""
    w = softmax(logprior)
    mu = Z[pool]                       # before any label: each hypothesis' own prediction, positive slope
    noise = np.median(se[pool] ** 2) + se[pool] ** 2   # before any label, misfit assumed ~ typical noise
    left = np.ones(len(pool), bool)
    order, k_stop, maps = [], kmax, []
    for step in range(kmax):
        mbar = mu @ w
        score = ((mu - mbar[:, None]) ** 2) @ w / np.maximum(noise, 1e-12)
        score[~left] = -np.inf
        j = int(np.argmax(score))
        left[j] = False
        order.append(pool[j])
        ll, A, B, s2 = posterior(Z[order], y[order], se[order], logprior)
        w = softmax(ll)
        mu = A + Z[pool] * B
        noise = s2 * (1 - a.beta) + se[pool] ** 2
        maps.append(int(np.argmax(ll)))
        prev = maps[-1 - a.stop_patience:-1]
        if (k_stop == kmax and len(order) >= a.min_labels and len(prev) == a.stop_patience
                and all(C[maps[-1], t] >= a.stop_r for t in prev)):
            k_stop = len(order)
    return order, k_stop


runs, t0 = [], time.time()
for element in elements:
    y, se, cls, Z_all, meta_all = load(element)
    n = len(y)
    scenarios = {"all": np.ones(len(meta_all), bool)}
    if meta_all.matched.any():
        scenarios["hide_matched"] = ~meta_all.matched.to_numpy()
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
        T = Z.shape[1]
        H = np.column_stack([Z] + [v / v.std() for v in ens.values()])
        hmeta = pd.concat([meta, pd.DataFrame(dict(biosample_name=list(ens), matched=[k.endswith("matched") for k in ens]),
                                              index=list(ens))])
        h_rho = np.concatenate([track_rho, rank_corr(H[:, T:], y)]) if ens else track_rho
        h_clusters = np.concatenate([clusters, clusters.max() + 1 + np.arange(len(ens))])
        logprior = np.full(H.shape[1], np.log((1 - a.ensemble_prior) / T) if ens else -np.log(T))
        if ens:
            logprior[T:] = np.log(a.ensemble_prior / len(ens))
        C = np.corrcoef(H.T)

        def label_free(name, pred):
            runs.append(dict(base, strategy=name, k=0, repeat=-1, labels=0, rho_pick=rank_corr(pred, y)[0],
                             acc_pick=np.mean(np.sign(pred[sig]) == cls[sig]), rho_oracle=rho_oracle))
        label_free("oracle", Z[:, best])
        label_free("all_mean", Z.mean(1))
        if scen == "all" and meta.matched.any():
            label_free("matched", Z[:, meta.matched.to_numpy()].mean(1))

        mean_pred = np.abs(Z.mean(1))
        spread = Z.var(1)
        for rep in range(a.repeats):
            rng = np.random.default_rng([a.seed, rep])
            pool = rng.choice(n, int(a.pool_frac * n), replace=False)
            order, k_stop = disc_sequence(H, y, se, pool, C, logprior)
            panels = {"random": rng.permutation(pool),
                      "high_effect": pool[np.argsort(-mean_pred[pool], kind="stable")],
                      "disagree": pool[np.argsort(-spread[pool], kind="stable")],
                      "disc": np.array(order)}
            todo = [(s, k, p[:k]) for s, p in panels.items() for k in budgets]
            todo.append(("disc_stop", kmax, np.array(order[:k_stop])))
            for strategy, k, panel in todo:
                r = evaluate(H, y, se, cls, panel, h_rho, logprior)
                t = r["pick"]
                runs.append(dict(base, strategy=strategy, k=k, repeat=rep, labels=len(panel),
                                 pick_track=hmeta.index[t], pick_biosample=hmeta.biosample_name.iloc[t],
                                 pick_matched=bool(hmeta.matched.iloc[t]), pick_ensemble=t >= T,
                                 ensemble_mass=float(r["w"][T:].sum()), pick_cluster_mass=float(
                                     r["w"][h_clusters == h_clusters[t]].sum()), exp_loss=float(1 - r["w"] @ C[:, t]),
                                 rho_pick=r["rho_pick"], rho_bma=r["rho_bma"], acc_pick=r["acc_pick"],
                                 acc_bma=r["acc_bma"], pick_pct=r["pick_pct"], rho_oracle=rho_oracle,
                                 regret=rho_oracle - h_rho[t]))
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
