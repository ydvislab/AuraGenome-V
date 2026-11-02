"""A one-file HTML report of a TrackPilot fit (no external resources)."""
from __future__ import annotations

import html

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .core import ANCHOR

CSS = """
:root { --bg: #fff; --fg: #1d1d1f; --muted: #6e6e73; --line: #e5e5ea; --bar: #2a6fdb; --anchor: #c2570c; }
@media (prefers-color-scheme: dark) { :root { --bg: #161618; --fg: #f2f2f7; --muted: #a1a1a6; --line: #333;
  --bar: #6ea1ff; --anchor: #f29b5c; } }
body { background: var(--bg); color: var(--fg); font: 14px/1.5 -apple-system, system-ui, sans-serif;
  max-width: 920px; margin: 32px auto; padding: 0 16px; }
h1 { font-size: 22px; } h2 { font-size: 16px; margin-top: 28px; } .muted { color: var(--muted); }
table { border-collapse: collapse; width: 100%; font-variant-numeric: tabular-nums; }
td, th { border-bottom: 1px solid var(--line); padding: 4px 8px; text-align: left; }
td.n, th.n { text-align: right; } .bar { height: 10px; background: var(--bar); border-radius: 0 3px 3px 0; }
.bar.a { background: var(--anchor); } .stats { display: flex; gap: 32px; flex-wrap: wrap; }
.stats div b { display: block; font-size: 22px; }
"""


def _rows(df, cols):
    out = []
    for _, r in df.iterrows():
        tds = []
        for c, fmt in cols:
            v = r[c]
            if fmt == "bar":
                cls = "bar a" if r["hypothesis"] == ANCHOR else "bar"
                tds.append(f'<td style="width:30%"><div class="{cls}" style="width:{100 * v:.1f}%"></div></td>')
            elif fmt:
                tds.append(f'<td class="n">{v:{fmt}}</td>')
            else:
                tds.append(f"<td>{html.escape(str(v))}</td>")
        out.append("<tr>" + "".join(tds) + "</tr>")
    return "\n".join(out)


def write_report(path, tp, fit, post, preds):
    fam = (post.groupby("family").agg(weight=("weight", "sum"), n=("weight", "size"),
                                      top=("hypothesis", "first"), biosample=("biosample_name", "first"))
           .sort_values("weight", ascending=False))
    n90 = int(np.searchsorted(fam.weight.cumsum().to_numpy(), 0.9)) + 1
    m = preds.measured.notna()
    fit_rho = spearmanr(preds.prediction[m], preds.measured[m])[0] if m.sum() > 2 else float("nan")
    anchor_w = float(post.weight[post.hypothesis == ANCHOR].sum())
    top = post.head(15).assign(label=lambda d: d.biosample_name.fillna("").astype(str))
    body = f"""<h1>TrackPilot report</h1>
<p class="muted">{len(tp.variants)} variants, {len(tp.tracks)} tracks,
{int(m.sum())} measured. Predictions are the posterior average over all hypotheses.</p>
<div class="stats"><div><b>{html.escape(str(fam.biosample.iloc[0]))}</b>leading track family
({fam.weight.iloc[0]:.2f} of the posterior, {fam.n.iloc[0]} tracks)</div>
<div><b>{n90}</b>track families hold 90% of the posterior</div>
<div><b>{anchor_w:.2f}</b>anchor weight{'' if ANCHOR in tp.hypotheses else ' (no anchor)'}</div>
<div><b>{fit_rho:.2f}</b>Spearman on the measured variants (in-sample)</div></div>
<h2>Track families</h2>
<table><tr><th>#</th><th>biosample of the top track</th><th>weight</th><th class="n">p</th><th class="n">tracks</th></tr>
{_rows(fam.head(10).reset_index(drop=True).assign(hypothesis=lambda d: d.top, rank=lambda d: d.index + 1),
       [("rank", ""), ("biosample", ""), ("weight", "bar"), ("weight", ".3f"), ("n", "d")])}</table>
<h2>Hypotheses</h2>
<table><tr><th>hypothesis</th><th>biosample</th><th>weight</th><th class="n">p</th><th class="n">slope</th></tr>
{_rows(top, [("hypothesis", ""), ("label", ""), ("weight", "bar"), ("weight", ".3f"), ("slope", ".3g")])}</table>
<p class="muted">A small in-sample Spearman says the measured variants are not explained by any track; a posterior
spread over many families says more measurements are needed before trusting one track.</p>"""
    with open(path, "w") as f:
        f.write(f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
                f'<meta name="viewport" content="width=device-width, initial-scale=1"><title>TrackPilot report</title>'
                f"<style>{CSS}</style></head><body>{body}</body></html>")
