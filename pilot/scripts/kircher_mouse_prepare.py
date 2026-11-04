"""Fourth test (docs/internal/kircher_mouse_zh.md): the six Kircher et al. experiments measured in mouse cell lines, never
analysed before. Human regulatory elements in Neuro-2a (RET, UC88), MIN6 (TCF7L2, ZFAND3) and NIH/3T3 (ZRS, two
experiments); no human track profiles these cells. A proxy is named from track metadata only: human neuroblastoma
lines for Neuro-2a, endocrine-pancreas precursors for MIN6, fibroblasts for NIH/3T3.

Reads pilot/data/variants_grch38.tsv and writes pilot/data/kircher_mouse/proxy_variants_{alphagenome,enformer}.tsv in the
format of data/kosicki (element = experiment, pred_element = element). No effect is compared with a prediction here.

usage (from the repository root): python pilot/scripts/kircher_mouse_prepare.py
"""
import os
import pandas as pd

OUT = "pilot/data/kircher_mouse"
PROXY = {"alphagenome": {"Neuro-2a": ["BE2C", "SK-N-DZ", "SK-N-SH"],
                         "MIN6": ["islet precursor cell", "progenitor cell of endocrine pancreas"],
                         "NIH/3T3": ["fibroblast"]},
         "enformer": {"Neuro-2a": ["BE2C", "SK-N-DZ", "SK-N-SH", "neuroblastoma cell line"],
                      "MIN6": ["islet precursor cell"],
                      "NIH/3T3": ["fibroblast"]}}
os.makedirs(OUT, exist_ok=True)
v = pd.read_csv("pilot/data/variants_grch38.tsv", sep="\t", low_memory=False)
v = v[~v.human_cells].copy()
v["cell_type"] = v.mpra_cell_line.str.split(" ").str[0]
for model, proxy in PROXY.items():
    out = v.assign(pred_element=v.element, se=0.0, related_biosamples="",
                   model_biosample_keywords=v.cell_type.map(lambda c: ";".join(proxy[c])))
    cols = ["element", "pred_element", "cell_type", "chrom", "pos", "ref", "alt", "effect", "se", "pvalue",
            "truth_class", "model_biosample_keywords", "related_biosamples"]
    out[cols].sort_values(["element", "pos", "alt"]).to_csv(f"{OUT}/proxy_variants_{model}.tsv", sep="\t", index=False)
print(v.groupby(["element", "cell_type"]).size().to_string())
