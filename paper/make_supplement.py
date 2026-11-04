"""Supplementary Tables S1-S6 from the result tables under results/a1 (run from the repository root):
    python paper/make_supplement.py && cd paper && tectonic supplement.tex
Writes paper/supp_tables.tex, which paper/supplement.tex inputs. No number is typed by hand."""
import re
import numpy as np
import pandas as pd

R = "results/a1"


def tsv(path):
    return pd.read_csv(f"{R}/{path}", sep="\t")


def f(x, d=3, sign=True):
    return "--" if pd.isna(x) else (f"{x:+.{d}f}" if sign else f"{x:.{d}f}").replace("-", "$-$")


def ci(lo, hi):
    return f"{f(lo)} to {f(hi)}"


def table(label, caption, cols, rows, spec=None):
    spec = spec or "l" * len(cols)
    body = "\n".join(" & ".join(str(c) for c in r) + r" \\" for r in rows)
    if len(rows) > 30:   # spans pages
        head = " & ".join(cols) + r" \\"
        return (f"{{\\small\\begin{{longtable}}{{{spec}}}\n\\caption{{{caption}}}\\label{{{label}}}\\\\\n"
                f"\\toprule\n{head}\n\\midrule\n\\endfirsthead\n\\toprule\n{head}\n\\midrule\n\\endhead\n"
                f"{body}\n\\bottomrule\n\\end{{longtable}}}}\n")
    return (f"\\begin{{table}}[!ht]\n\\caption{{{caption}}}\\label{{{label}}}\n\\centering\\small\n"
            f"\\begin{{tabular}}{{{spec}}}\n\\toprule\n{' & '.join(cols)} \\\\\n\\midrule\n{body}\n"
            f"\\bottomrule\n\\end{{tabular}}\n\\end{{table}}\n")


def decision(path):
    """P1-P4 lines of a siraj_decision.txt -> {name: (estimate, test statistic)}"""
    out = {}
    for line in open(f"{R}/{path}"):
        m = re.match(r"(P\d(?: \(S4, secondary\))?) (?:\S+ )?k=(\d+): (.*?) -> (\w+)", line)
        if m:
            name, k, text, verdict = m.groups()
            num = re.findall(r"[-+]?\d+\.\d+(?:e-?\d+)?|\d+/\d+", text)
            out[f"{name[:2]} $k={k}$"] = (num, verdict)
    return out


parts = []

# S1: the pre-registered Siraj test with labels relative to the reference oligo and with the published centring
ref, pub = decision("siraj_ref/siraj_decision.txt"), decision("siraj_decision.txt")
what = {"P1": "mean D (lower bound)", "P2": "harmed (upper bound)", "P3": "v2 $-$ all-track mean ($P$)",
        "P4": "v2 $-$ related mean ($P$)"}
rows = []
for key in ref:
    def cell(n):
        return f"{f(float(n[0]))} ({f(float(n[1]))})" if key[:2] == "P1" else \
            f"{n[1]} ({float(n[-2]):.3f})" if key[:2] == "P2" else f"{f(float(n[0]))} ($P$={float(n[-1]):.0e})"
    rows.append([key, what[key[:2]], cell(ref[key][0]), cell(pub[key][0])])
s_ref, s_pub = tsv("siraj_ref/siraj_secondary.tsv"), tsv("siraj_secondary.tsv")
for _, r in s_ref.iterrows():
    p = s_pub[s_pub.comparison == r.comparison]
    name = r.comparison.split(":")[0]
    rows.append([name, r.comparison.split(": ")[1].replace("_", r"\_").replace(" - ", " $-$ "),
                 f"{f(r['mean'])} [{ci(r.ci_low, r.ci_high)}]",
                 f"{f(p['mean'].iloc[0])} [{ci(p.ci_low.iloc[0], p.ci_high.iloc[0])}]" if len(p) else "--"])
parts.append(table("tab:s1", "Siraj test (125 loci): effects relative to the measured reference oligo (primary) and "
                   "centred on all oligos of the element as published (original analysis). P1--P4 are the "
                   "pre-registered hypotheses; all held under both definitions. Secondary rows give the mean "
                   "difference over loci and its 95\\% bootstrap interval. S3b was added with the reference-relative "
                   "labels.", ["", "Quantity", "Reference-relative", "Published centring"], rows, "llll"))

# S2: directional accuracy
d = tsv("direction/direction_contrasts.tsv")
names = {"siraj": "Siraj", "kircher_v1": "Kircher (v1)", "kircher_v2": "Kircher (v2)",
         "kircher_v2_no_offset": "Kircher (v2, no offset)", "kircher_v2_with_offset_c3": "Kircher (v2)"}
rows = [[names[r.set], r.contrast, r.first.replace("_", r"\_"), r.second.replace("_", r"\_"), r.units,
         f(r.acc_first, sign=False), f(r.acc_second, sign=False), f(r["mean"]), ci(r.ci_low, r.ci_high),
         f"{r.positive}/{r.negative}"] for _, r in d.iterrows()]
parts.append(table("tab:s2", "Directional accuracy (post hoc): fraction of significant held-out variants "
                   "($P<10^{-5}$, $|\\text{effect}|\\ge0.1$) whose predicted sign is correct. Labelled methods use "
                   "the posterior-average prediction; \\emph{no offset} drops the fitted offsets (diagnostic). "
                   "Interval: 95\\% bootstrap over units; last column, units better/worse.",
                   ["Set", "Contrast", "First", "Second", "$n$", "Acc.\\,1", "Acc.\\,2", "Diff.", "95\\% CI", "+/$-$"],
                   rows, "llllrrrrll"))

# S3: Siraj with the held-out-fold checkpoints
c1 = decision("siraj_ref/c1/siraj_decision.txt")
delta = tsv("siraj_ref/c1/c1_siraj_delta.tsv")
lab = {"C1_D_k8": "P1 $k=8$", "C1_D_k16": "P1 $k=16$", "C3_D_k8": "P3 $k=8$", "C2_D_k16": "P4 $k=16$"}
rows = []
for _, r in delta.iterrows():
    n, verdict = c1[lab[r.contrast]]
    rows.append([lab[r.contrast], f(r.all_folds), f(r.held_out), f"{f(r.change)} [{ci(r.ci_low, r.ci_high)}]",
                 verdict.lower()])
for key in ["P2 $k=8$", "P2 $k=16$"]:
    rows.append([key, f"{decision('siraj_ref/siraj_decision.txt')[key][0][1]} harmed", f"{c1[key][0][1]} harmed",
                 "--", c1[key][1].lower()])
parts.append(table("tab:s3", "Siraj test repeated with, for each locus, the AlphaGenome checkpoint that held out its "
                   "genomic fold (labels relative to the reference oligo). Columns give the pre-registered "
                   "contrast with the distilled checkpoint (ALL\\_FOLDS) and with the held-out checkpoint, and the "
                   "paired change with its 95\\% bootstrap interval over loci.",
                   ["Hypothesis", "ALL\\_FOLDS", "Held-out fold", "Change [95\\% CI]", "Verdict"], rows, "lrrll"))

# S4: acquisition comparators
rows = []
for name, path in [("Siraj", "siraj_ref/cmp/cmp_siraj_contrasts.tsv"), ("Kircher", "cmp/cmp_test_contrasts.tsv")]:
    t = tsv(path)
    for _, r in t.iterrows():
        rows.append([name, r.view.replace("_", r"\_"), r.k, r.contrast.replace("_", r"\_").replace(" - ", " $-$ "),
                     r.units, f(r["mean"]), ci(r.ci_low, r.ci_high), r.reading.replace("within 0.02", "")])
parts.append(table("tab:s4", "Acquisition comparators (post hoc; frozen v2 inference, 50 repeats). Siraj views C1--C3; "
                   "Kircher test regions with matched tracks kept (deploy) or without an anchor (no\\_anchor). "
                   "Reading: 95\\% interval within $\\pm0.02$ = equivalent; above or below 0 = first or second better.",
                   ["Data", "View", "$k$", "Contrast", "$n$", "Mean", "95\\% CI", "Reading"], rows, "llrlrrll"))

# S5: CODA / MODEL SELECTOR style baselines
rows = []
for name, path in [("Siraj", "coda/coda_siraj_contrasts.tsv"), ("Kircher", "coda/coda_test_contrasts.tsv")]:
    t = tsv(path)
    view = {"all": "C1", "hide_exact": "C2", "hide_matched": "C3"}
    for _, r in t.iterrows():
        rows.append([name, view[r.view], r.k, r.contrast.replace("_", r"\_").replace(" - ", " $-$ "), r.units,
                     f(r["mean"]), ci(r.ci_low, r.ci_high), r.positive, r.reading.replace("within 0.02", "")])
parts.append(table("tab:s5", "CODA- and MODEL SELECTOR-style baselines (post hoc; 50 repeats). Both use the calibrated "
                   "likelihood, no metadata and the information-gain panel; CODA uses a consensus prior "
                   "($\\kappa=10$, fixed on development elements), MODEL SELECTOR a uniform prior ($\\kappa=0$). "
                   "\\emph{Single track}: the most probable track, the methods' own output; \\emph{average}: their "
                   "posterior average. Positive: units where the first is better.",
                   ["Data", "View", "$k$", "Contrast", "$n$", "Mean", "95\\% CI", "+", "Reading"], rows,
                   "llrlrrlrl"))

# S6: EIG convergence
cov = tsv("eigconv/eigconv_siraj_cover.tsv").assign(data="Siraj")
cov = pd.concat([cov, tsv("eigconv/eigconv_test_cover.tsv").assign(data="Kircher")])
cov = cov[(cov.k == 8) & (cov.quantity.isin(["eig_cover_first", "eig_cover_min"]))]
rows = []
for (data, view), g in cov.groupby(["data", "view"], sort=False):
    g = g.set_index(["setting", "quantity"])
    rows.append([data, view.replace("_", r"\_")] +
                [f"{g.loc[(s, q), 'median']:.2f} ({g.loc[(s, q), 'minimum']:.2f})"
                 for s in ["E0", "E1"] for q in ["eig_cover_first", "eig_cover_min"]])
parts.append("\\renewcommand{\\thetable}{S6a}\n")
parts.append(table("tab:s6a", "EIG convergence check (post hoc), coverage: posterior mass held by the hypotheses kept, "
                   "median over units (minimum), at the first step and the lowest over the eight steps of a $k=8$ "
                   "panel. E0: the 100 most probable hypotheses, 96 quadrature points (as in the main text); E1: 300 "
                   "hypotheses; E2 (not shown) keeps $\\ge$99\\% of the mass at every step with 192 points.",
                   ["Data", "View", "E0 first", "E0 lowest", "E1 first", "E1 lowest"], rows, "llrrrr"))
rows = []
for name, path in [("Siraj", "eigconv/eigconv_siraj_contrasts.tsv"), ("Kircher", "eigconv/eigconv_test_contrasts.tsv")]:
    t = tsv(path)
    for _, r in t.iterrows():
        rows.append([name, r.view.replace("_", r"\_"), r.k, r.contrast.replace(" - ", " $-$ "), r.units,
                     f(r["mean"], 4), ci(r.ci_low, r.ci_high), r.reading.replace("within 0.02", "")])
parts.append("\\addtocounter{table}{-1}\\renewcommand{\\thetable}{S6b}\n")
parts.append(table("tab:s6b", "EIG convergence check, accuracy: differences in Spearman correlation of the "
                   "posterior-average prediction (same pools and 50 repeats; 95\\% bootstrap interval over units).",
                   ["Data", "View", "$k$", "Contrast", "$n$", "Mean", "95\\% CI", "Reading"], rows, "llrlrrll"))

# S7: label efficiency and where the gain is (post hoc, descriptive; pilot/scripts/reframe_report.py)
parts.append("\\renewcommand{\\thetable}{S\\arabic{table}}\n")
view = {"all": "C1 profiled", "hide_exact": "C2 related only", "hide_matched": "C3 unprofiled"}
eff, gap = tsv("reframe/label_efficiency.tsv"), tsv("reframe/gain_by_gap.tsv")
rows = [[view[r.scenario], int(r.k), f(r.pilot, sign=False), f(r.random, sign=False), f(r.default, sign=False),
         f(r.oracle, sign=False), f"{100 * r.gap_closed:.0f}\\%"] for _, r in eff.iterrows()]
parts.append(table("tab:s7", "Label efficiency on Siraj (post hoc, descriptive; 125 loci, 200 repeats). Mean Spearman "
                   "correlation of the posterior-average prediction after $k$ pilot or random measurements, the "
                   "label-free default and the best single track chosen with all labels (oracle). Gap closed: share "
                   "of the default-to-oracle gap recovered by the pilot panel.",
                   ["View", "$k$", "Pilot", "Random", "Default", "Oracle", "Gap closed"], rows, "lrrrrrr"))
rows = [[view[r.scenario], int(r.quartile), int(r.loci), f(r.gap, sign=False), f(r.gain_k8),
         f"{100 * r.share_closed:.0f}\\%", int(r.worse_by_005)] for _, r in gap.iterrows()]
parts.append(table("tab:s8", "Where the gain is (post hoc, descriptive): Siraj loci split into quartiles of the gap "
                   "between the label-free default and the oracle track (which needs all labels), gain of the pilot "
                   "panel over the default at $k=8$, and loci worse than the default by more than 0.05.",
                   ["View", "Quartile", "Loci", "Gap", "Gain", "Share closed", "Worse $>0.05$"], rows, "lrrrrrr"))
# S9: cross-locus prior (pre-registered extension; pilot/scripts/cross_locus.py)
dec = open(f"{R}/cross_locus/decision.txt").read().splitlines()
nums = [re.findall(r"[-+]?\d+\.\d+(?:e-?\d+)?|\d+/\d+", l) for l in dec[1:3]]
pfmt = lambda p: "{:.0f}\\times10^{{{:d}}}".format(p / 10 ** int(np.floor(np.log10(p))), int(np.floor(np.log10(p))))
x1x2 = (f"X1 (C3, $k=8$, cross-locus minus single locus at least 0.02): {nums[0][0]}, {nums[0][2]} loci positive, "
        f"one-sided Wilcoxon $P={pfmt(float(nums[0][3]))}$; not met. X2 (C1, $k=8$, non-inferior to the matched mean): "
        f"{nums[1][0]}, lower bound {nums[1][1]}, {nums[1][3]} loci worse by more than 0.05 (upper bound "
        f"{nums[1][4]}); met.")
sec = tsv("cross_locus/secondary.tsv")
rows = []
for _, r in sec.iterrows():
    m = "10" if pd.isna(r.m) else str(int(r.m))
    if pd.isna(r.m) or r.m != 10:
        rows.append([view[r.scenario], int(r.k), m, "--" if pd.isna(r.single) else f(r.single, sign=False),
                     f(r.cross, sign=False), f(r.cross_minus_single, 4),
                     "--" if pd.isna(r.positive) else f"{int(r.positive)}/{int(r.loci)}"])
parts.append(table("tab:s9", "Cross-locus prior (pre-registered extension; 85 test loci, 20 repeats, $\\lambda=0.75$ "
                   "chosen on 40 development loci). The prior of a new locus mixes the frozen prior with the mean "
                   "posterior of $m$ other loci in the same cell type after their own 8-variant pilot panels. "
                   + x1x2 + " ",
                   ["View", "$k$", "$m$", "Single locus", "Cross-locus", "Difference", "Loci $>0$"], rows, "lrrrrrr"))

# S10: neuron MPRA (pre-registered third dataset; pilot/scripts/kosicki_report.py)
lev = {m: tsv(f"kosicki_{m}/levels.tsv").set_index(["scenario", "k"]) for m in ["alphagenome", "enformer"]}
verdicts = []
for m, name in [("alphagenome", "AlphaGenome"), ("enformer", "Enformer")]:
    lines = open(f"{R}/kosicki_{m}/decision.txt").read().splitlines()[1:5]
    n = [re.findall(r"[-+]?\d+\.\d+(?:e-?\d+)?|\d+/\d+", l) for l in lines]
    held = ["held" if "HOLDS" in l else "not met" for l in lines]
    verdicts.append(f"{name}: N1 {f(float(n[0][0]))} ({held[0]}), N2 {f(float(n[1][0]))} ($P={float(n[1][-1]):.2g}$, "
                    f"{held[1]}), N3 {f(float(n[2][0]))}, lower bound {f(float(n[2][1]))}, {n[2][3]} regions worse "
                    f"by more than 0.05 ({held[2]}), N4 {f(float(n[3][0]))} ({held[3]})")
rows = [[view[s], int(k)] + [f(lev[m].loc[(s, k), c], sign=False) for m in lev for c in ["pilot", "random", "default"]]
        for s, k in lev["alphagenome"].index]
orc = " and ".join(f"{lev[m].oracle.mean():.2f}" for m in lev)
parts.append(table("tab:s10", "Neuron MPRA (pre-registered third dataset; WTC11-NGN2 induced excitatory neurons, 159 "
                   "regions, 200 repeats). Mean Spearman correlation of TrackPilot after $k$ pilot or random measurements "
                   "and of the label-free default (matched, related or all-track mean). "
                   + "; ".join(verdicts) + ". Decision (N1 and N3): not supported for either model. The best single "
                   f"track chosen with all labels reached {orc} (AlphaGenome, Enformer), mostly selection among about "
                   "1,000 tracks on about 34 held-out variants per tile.",
                   ["View", "$k$", "AG pilot", "AG random", "AG default", "Enf.\\ pilot", "Enf.\\ random",
                    "Enf.\\ default"], rows, "lrrrrrrr"))

# S11: applicability check (pre-registered extension; pilot/scripts/applicability_check.py)
dec = open(f"{R}/applicability/decision.txt").read().splitlines()
sec = tsv("applicability/secondary.tsv")
name = {"kircher": "Kircher (11 regions)", "siraj": "Siraj (85 loci)", "neuron": "Neurons (109 regions)"}
rows = [[name[r.dataset], view[r.scenario], f(r.rho_tp, sign=False), f(r.gated, sign=False), f(r.rho_default, sign=False),
         f"{100 * r.share_no_signal:.0f}\\%"] for _, r in sec.iterrows()]
parts.append(table("tab:s11", "Applicability check (pre-registered extension; test units, $k=8$, 20 repeats, threshold "
                   "$P>0.5$ chosen on development units). After the pilot panel is measured, a permutation test (200 "
                   "shuffles of the eight labels) compares the mixture's log evidence with that of shuffled labels; if "
                   "the labels are explained no better, the label-free default is returned. Mean Spearman correlation "
                   "of TrackPilot without and with the check, of the default, and the share of runs judged to carry no "
                   "trackable signal. A1 (neurons: checked $-$ unchecked $>0$; checked non-inferior to the default, "
                   "margin 0.02) and A2 (Siraj: gains over the defaults kept; checked non-inferior to unchecked, "
                   "margin 0.01) both held; see the repository for the exact statistics. Secondary, neuron test regions with "
                   "Enformer predictions: checked $-$ unchecked $+0.019$ ($P=2\\times10^{-4}$), checked $-$ default $-0.015$ "
                   "(lower bound $-0.020$).",
                   ["Data", "View", "Unchecked", "Checked", "Default", "No signal"], rows, "llrrrr"))

open("paper/supp_tables.tex", "w").write("\n".join(parts))
print("wrote paper/supp_tables.tex")
