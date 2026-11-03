"""Third dataset (docs/internal/neuron_mpra_zh.md): lentiMPRA in WTC11-NGN2 induced excitatory neurons
(Kosicki et al. 2025, Nat Commun; ENCODE ENCSR548AQS, CC BY 4.0), synthetic single-base transversions every fourth
base of 270-bp tiles.

Reads pilot/data/raw/kosicki2025/ (ENCFF958ZLG.tsv element quantifications, ENCFF752NFF.fasta.gz designed sequences)
and writes pilot/data/kosicki/proxy_variants_{alphagenome,enformer}.tsv in the format of data/siraj_ref:
  effect      = ratio(mutant tile) - ratio(reference tile), the reference tile measured in the same experiment
  element     = the tile (KOS_<tile id>); pred_element = the region of overlapping tiles (KOSR<nnn>), scored once
  matched     = model_biosample_keywords (substring of biosample + track name): neuron tracks
  related     = related_biosamples (exact names): neural and brain biosamples other than the matched ones
Rules fixed before any prediction exists: '+'-strand tiles only (the '-' tiles hold 172 SNVs whose positions do not
match the reference sequence); tiles with at least 30 quantified SNVs; reference base checked against the designed
sequence. Only names, positions, sequences and track metadata are read here; no effect is compared with a prediction.

usage (from the repository root): python pilot/scripts/kosicki_prepare.py
"""
import os, re, gzip, glob
import numpy as np
import pandas as pd

RAW, OUT = "pilot/data/raw/kosicki2025", "pilot/data/kosicki"
MIN_SNV = 30
os.makedirs(OUT, exist_ok=True)

q = pd.read_csv(f"{RAW}/ENCFF958ZLG.tsv", sep="\t")
p = q.name.str.split("|", expand=True)
c = p[1].str.extract(r"^(chr\w+)-(\d+)([+-])(\d+)$")
seq, name = {}, None
with gzip.open(f"{RAW}/ENCFF752NFF.fasta.gz", "rt") as f:
    for line in f:
        if line.startswith(">"):
            name = line[1:].strip()
            seq[name] = []
        else:
            seq[name].append(line.strip())
seq = {k: "".join(v).upper() for k, v in seq.items()}

ref_ratio = pd.Series(q.ratio[p[2] == "reference"].to_numpy(), index=p[1][p[2] == "reference"].to_numpy())
ref_ratio = ref_ratio[~ref_ratio.index.duplicated()]
snv = p[2].str.match(r"^[ACGT]-[ACGT]-\d+$") & (c[2] == "+")
v = pd.DataFrame({"tile": p[0][snv], "coord": p[1][snv], "chrom": c[0][snv], "start": c[1][snv].astype(int),
                  "ratio": q.ratio[snv]})
v[["ref", "alt", "offset"]] = p[2][snv].str.split("-", expand=True)
v["offset"] = v.offset.astype(int)
v["pos"] = v.start + v.offset - 1
refseq = {k.split("|")[1]: s for k, s in seq.items() if k.endswith("|reference")}
v = v[np.array([refseq[t][o - 1] == r for t, o, r in zip(v.coord, v.offset, v.ref)]) & (v.ref != v.alt).to_numpy()]
v = v[v.coord.isin(ref_ratio.index)]
n = v.groupby("coord").size()
v = v[v.coord.isin(n[n >= MIN_SNV].index)].copy()
v["effect"] = v.ratio - v.coord.map(ref_ratio)

# regions of overlapping tiles, scored once per region
tiles = v.groupby("coord").agg(chrom=("chrom", "first"), start=("start", "first")).reset_index()
tiles["end"] = tiles.start + 269
tiles = tiles.sort_values(["chrom", "start"])
region, rid, last = {}, 0, (None, -1)
for t in tiles.itertuples():
    if t.chrom != last[0] or t.start > last[1]:
        rid += 1
        last = (t.chrom, t.end)
    else:
        last = (t.chrom, max(last[1], t.end))
    region[t.coord] = f"KOSR{rid:03d}"
v["pred_element"] = v.coord.map(region)
v["element"] = "KOS_" + v.tile

# matched and related tracks, from each model's track metadata only
INC = re.compile(r"neur|brain|cortex|cerebel|hippocamp|WTC11|glutamat|spinal cord|medulla|putamen|caudate|thalam|pons|"
                 r"substantia", re.I)
EXC = re.compile(r"renal|kidney|adrenal|astrocyte|pericyte|endothelial|perineur|smooth muscle|crest|glioblastoma|glioma|"
                 r"neuroblastoma|tumor|neuroepithelioma|neurofibroma", re.I)
MATCHED = {"alphagenome": ["CL:0000540", "glutamatergic neuron", "motor neuron"],   # 'neuron' (CAGE), glutamatergic
           "enformer": ["bipolar neuron", "Neurons,"]}                              # neuron, motor neuron; iPSC/primary
META = {"alphagenome": glob.glob("pilot/predictions/siraj_all/*ALL_FOLDS*")[0],
        "enformer": glob.glob("pilot/predictions/enformer_siraj/*")[0]}
for model, path in META.items():
    m = pd.read_csv(path, sep="\t", usecols=["track_name", "biosample_name"]).drop_duplicates()
    text = m.biosample_name.fillna("") + " " + m.track_name
    is_matched = np.any([text.str.contains(k, case=False, regex=False) for k in MATCHED[model]], axis=0)
    b = m.biosample_name.fillna("")
    related = sorted(set(b[b.str.contains(INC) & ~b.str.contains(EXC) & ~is_matched]) - set(b[is_matched]))
    out = v.assign(cell_type="iN", se=0.0, pvalue=np.nan, truth_class=0,
                   model_biosample_keywords=";".join(MATCHED[model]), related_biosamples=";".join(related))
    cols = ["element", "pred_element", "cell_type", "chrom", "pos", "ref", "alt", "effect", "se", "pvalue",
            "truth_class", "model_biosample_keywords", "related_biosamples"]
    out[cols].sort_values(["pred_element", "element", "pos"]).to_csv(f"{OUT}/proxy_variants_{model}.tsv", sep="\t",
                                                                     index=False)
    print(f"{model}: {int(is_matched.sum())} matched track rows, {len(related)} related biosamples")
print(f"{v.element.nunique()} tiles, {v.pred_element.nunique()} regions, {len(v)} SNV measurements, "
      f"{v.drop_duplicates(['chrom', 'pos', 'alt']).shape[0]} unique variants")
