"""Pilot step 1: variant table with measured effects (ground truth) from Kircher et al. 2019.

Input : data/raw/elements.tsv.gz  (official data of the MPRA portal, github.com/kircherlab/MPRA_SaturationMutagenesis)
Output: data/variants_grch38.tsv  one row per SNV, GRCh38, with
        element, element_type, cell line used in the MPRA, whether that cell line is human,
        model biosample keywords to look up in model outputs (human cells only; proxies flagged),
        measured log2 effect, p-value, and a ground-truth class:
            class = +1 / -1 if p < P_SIG and |effect| >= MIN_EFFECT (direction of the measured effect)
                     0      otherwise (no detectable effect)
        P_SIG and MIN_EFFECT are parameters; defaults follow a CAGI5-style convention and must be
        confirmed against the CAGI5 regulation-saturation challenge description before use.
Deletions ('-') are kept in a separate file because not every model scorer handles them alike.

Data use: cite Kircher M, et al. (2019) Nat Commun 10:3583, doi:10.1038/s41467-019-11526-w.

usage: python build_variant_table.py <pilot_dir> [--p-sig 1e-5] [--min-effect 0.1]
"""
import sys, os, argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("pilot_dir")
ap.add_argument("--p-sig", type=float, default=1e-5)
ap.add_argument("--min-effect", type=float, default=0.1)
a = ap.parse_args()
RAW = os.path.join(a.pilot_dir, "data", "raw")

# element -> (type, MPRA cell line, human?, keywords to match model biosample names, proxy?)
# From the portal's element tables. Cell lines that models do not cover are mapped to the closest
# human biosample and flagged as proxy; mouse cell lines are flagged and excluded from the
# cell-type-matched analysis (they stay in the table for the cell-type-agnostic analysis).
ELEMENTS = {
    "F9": ("promoter", "HepG2", True, ["HepG2"], False),
    "FOXE1": ("promoter", "HeLa", True, ["HeLa"], False),
    "GP1BA": ("promoter", "HEL 92.1.7", True, ["K562"], True),
    "HBB": ("promoter", "HEL 92.1.7", True, ["K562"], True),
    "HBG1": ("promoter", "HEL 92.1.7", True, ["K562"], True),
    "HNF4A": ("promoter", "HEK293T", True, ["HEK293"], False),
    "LDLR": ("promoter", "HepG2", True, ["HepG2"], False),
    "LDLR.2": ("promoter", "HepG2", True, ["HepG2"], False),
    "MSMB": ("promoter", "HEK293T", True, ["HEK293"], False),
    "PKLR-24h": ("promoter", "K562", True, ["K562"], False),
    "PKLR-48h": ("promoter", "K562", True, ["K562"], False),
    "TERT-HEK": ("promoter", "HEK293T", True, ["HEK293"], False),
    "TERT-GBM": ("promoter", "SF7996 (glioblastoma)", True, ["SK-N-SH", "astrocyte"], True),
    "TERT-GAa": ("promoter", "glioblastoma (alt.)", True, ["SK-N-SH", "astrocyte"], True),
    "TERT-GSc": ("promoter", "glioblastoma (alt.)", True, ["SK-N-SH", "astrocyte"], True),
    "BCL11A": ("enhancer", "HEL 92.1.7", True, ["K562"], True),
    "IRF4": ("enhancer", "SK-MEL-28", True, ["melanocyte", "SK-MEL"], True),
    "IRF6": ("enhancer", "HaCaT", True, ["keratinocyte"], True),
    "MYCrs11986220": ("enhancer", "LNCaP + DHT", True, ["LNCaP"], False),
    "MYCrs6983267": ("enhancer", "HEK293T", True, ["HEK293"], False),
    "SORT1": ("enhancer", "HepG2", True, ["HepG2"], False),
    "SORT1.2": ("enhancer", "HepG2", True, ["HepG2"], False),
    "SORT1-flip": ("enhancer", "HepG2", True, ["HepG2"], False),
    "RET": ("enhancer", "Neuro-2a (mouse)", False, [], True),
    "UC88": ("enhancer", "Neuro-2a (mouse)", False, [], True),
    "TCF7L2": ("enhancer", "MIN6 (mouse)", False, [], True),
    "ZFAND3": ("enhancer", "MIN6 (mouse)", False, [], True),
    "ZRSh-13": ("enhancer", "NIH/3T3 (mouse)", False, [], True),
    "ZRSh-13h2": ("enhancer", "NIH/3T3 (mouse)", False, [], True),
}

d = pd.read_csv(os.path.join(RAW, "elements.tsv.gz"), sep="\t", low_memory=False)
d = d[d.Release == "GRCh38"].copy()
missing = sorted(set(d.Element) - set(ELEMENTS))
assert not missing, f"elements without metadata: {missing}"
meta = pd.DataFrame([(k, *v) for k, v in ELEMENTS.items()],
                    columns=["element", "element_type", "mpra_cell_line", "human_cells", "model_biosample_keywords", "biosample_is_proxy"])
d = d.rename(columns={"Chrom": "chrom", "Pos": "pos", "Ref": "ref", "Alt": "alt", "Barcodes": "barcodes",
                      "DNA": "dna", "RNA": "rna", "Coefficient": "effect", "pValue": "pvalue", "Element": "element"})
# a few rows have no GRCh38 coordinates in the source file (11 MSMB rows); they cannot be scored
no_coord = d.chrom.isna() | d.pos.isna()
if no_coord.any():
    print(f"dropping {int(no_coord.sum())} rows without GRCh38 coordinates: {d[no_coord].element.value_counts().to_dict()}")
d = d[~no_coord].copy()
d["pos"] = d.pos.astype(int)
d["chrom"] = "chr" + d.chrom.astype(str).str.replace(r"\.0$", "", regex=True)
d = d.merge(meta, on="element", how="left")
d["model_biosample_keywords"] = d.model_biosample_keywords.map(lambda k: ";".join(k))
sig = (d.pvalue < a.p_sig) & (d.effect.abs() >= a.min_effect)
d["truth_class"] = 0
d.loc[sig & (d.effect > 0), "truth_class"] = 1
d.loc[sig & (d.effect < 0), "truth_class"] = -1
d["variant_id"] = d.chrom + ":" + d.pos.astype(str) + ":" + d.ref + ">" + d.alt + ":" + d.element
cols = ["variant_id", "chrom", "pos", "ref", "alt", "element", "element_type", "mpra_cell_line", "human_cells",
        "model_biosample_keywords", "biosample_is_proxy", "effect", "pvalue", "truth_class", "barcodes", "dna", "rna"]
snv, dele = d[d.alt != "-"][cols], d[d.alt == "-"][cols]
out = os.path.join(a.pilot_dir, "data")
snv.to_csv(os.path.join(out, "variants_grch38.tsv"), sep="\t", index=False)
dele.to_csv(os.path.join(out, "deletions_grch38.tsv"), sep="\t", index=False)

# summary for the pilot README
s = snv.groupby(["element", "element_type", "mpra_cell_line", "human_cells", "biosample_is_proxy"]).agg(
    n=("variant_id", "size"), up=("truth_class", lambda x: int((x == 1).sum())), down=("truth_class", lambda x: int((x == -1).sum()))).reset_index()
s.to_csv(os.path.join(out, "element_summary.tsv"), sep="\t", index=False)
print(f"SNVs {len(snv)} | deletions {len(dele)} | elements {snv.element.nunique()}")
print(f"truth classes (p<{a.p_sig}, |effect|>={a.min_effect}): {snv.truth_class.value_counts().to_dict()}")
print(f"human-cell, non-proxy SNVs (cell-type-matched analysis): {int((snv.human_cells & ~snv.biosample_is_proxy).sum())}")
print(s.to_string(index=False))
