"""Enformer variant-effect scores for every MPRA SNV (second foundation model; docs/internal/enformer_generalisation_zh.md).

The output has the tidy format of alphagenome_score.py, so proxy_select.py reads it with --version ENFORMER.
Model: enformer-pytorch with the EleutherAI/enformer-official-rough weights, human head, CAGE, DNase and ATAC tracks
(pilot/data/enformer/targets_human.txt). Input: 196,608 bp of GRCh38 centred on the element (UCSC API, cached under
pilot/data/enformer/seq/). One forward pass per sequence (batch 1 on MPS: batches > 1 give wrong outputs there; batch 1 matches the CPU
to 0.005), no reverse-complement or shift averaging. Score of a variant
on a track: log2(1 + sum of alt predictions) - log2(1 + sum of ref predictions) over the 128-bp bins that overlap the
501 bp centred on the variant, the analogue of AlphaGenome's CenterMaskScorer(width=501, DIFF_LOG2_SUM).

usage (from the repository root, with .venv-enformer):
  python pilot/scripts/enformer_score.py pilot --variants data/siraj_ref/proxy_variants.tsv \
      --out predictions/enformer_siraj --elements <pred_element,...>
  python pilot/scripts/enformer_score.py pilot --out predictions/enformer_kircher --elements IRF4,LDLR,...
Output: <pilot_dir>/<out>/<element>__ENFORMER.tsv.gz
"""
import os, sys, json, time, argparse, urllib.request
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("pilot_dir")
ap.add_argument("--variants", default="data/variants_grch38.tsv", help="relative to pilot_dir")
ap.add_argument("--out", required=True, help="relative to pilot_dir")
ap.add_argument("--elements", default="", help="comma-separated prediction elements; default all human-cell ones")
ap.add_argument("--batch", type=int, default=1, help="MPS gives wrong outputs for batches > 1 (torch 2.x)")
ap.add_argument("--device", default="mps")
ap.add_argument("--limit", type=int, default=0, help="smoke test: score only the first N variants per region")
a = ap.parse_args()
assert a.batch == 1 or a.device != "mps", "MPS gives wrong outputs for batches > 1"

SEQ, CROP, BIN, HALF = 196608, 40960, 128, 250
var = pd.read_csv(os.path.join(a.pilot_dir, a.variants), sep="\t", low_memory=False)
if "pred_element" not in var:          # Kircher: one prediction file per experiment
    var = var[var.human_cells].assign(pred_element=lambda d: d.element)
var["variant_id"] = var.chrom + ":" + var.pos.astype(str) + ":" + var.ref + ">" + var.alt
if a.elements:
    var = var[var.pred_element.isin(a.elements.split(","))]
var = var.drop_duplicates(["pred_element", "variant_id"])
reg = var.groupby("pred_element").agg(chrom=("chrom", "first"), lo=("pos", "min"), hi=("pos", "max")).reset_index()
reg["region"] = reg.chrom + ":" + reg.lo.astype(str) + "-" + reg.hi.astype(str)
out_dir = os.path.join(a.pilot_dir, a.out)
seq_dir = os.path.join(a.pilot_dir, "data", "enformer", "seq")
os.makedirs(out_dir, exist_ok=True)
os.makedirs(seq_dir, exist_ok=True)

tg = pd.read_csv(os.path.join(a.pilot_dir, "data", "enformer", "targets_human.txt"), sep="\t")
tg["output_type"] = tg.description.str.split(":").str[0]
tg = tg[tg.output_type.isin(["CAGE", "DNASE", "ATAC"])]
tg["biosample_name"] = tg.description.str.split(":", n=1).str[1].str.strip()
tg["track_name"] = tg.identifier + " " + tg.description             # descriptions repeat across replicates
idx = tg.index.to_numpy()
print(f"{len(reg)} elements in {reg.region.nunique()} regions, {len(var)} SNVs, {len(tg)} tracks", flush=True)


def fetch(chrom, start):
    """GRCh38 [start, start + SEQ), 0-based, upper case, cached."""
    f = os.path.join(seq_dir, f"{chrom}_{start}.txt")
    if not os.path.exists(f):
        url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={start + SEQ}"
        dna = json.load(urllib.request.urlopen(url, timeout=120))["dna"].upper()
        assert len(dna) == SEQ, (chrom, start, len(dna))
        open(f, "w").write(dna)
    return open(f).read()


import torch  # noqa: E402
from enformer_pytorch import from_pretrained  # noqa: E402
model = from_pretrained("EleutherAI/enformer-official-rough", use_tf_gamma=False).to(a.device).eval()
LUT = np.full(256, 4, dtype=np.int64)       # enformer-pytorch index encoding: ACGT = 0..3, N = 4
for i, b in enumerate("ACGT"):
    LUT[ord(b)] = i


def predict(seqs):
    x = torch.from_numpy(np.stack([LUT[np.frombuffer(s.encode(), np.uint8)] for s in seqs])).to(a.device)
    with torch.inference_mode():
        return model(x, head="human")[:, :, idx].float().cpu().numpy()


# experiments that measured the same region (e.g. the SORT1 or TERT experiments) share one scoring run
for region, grp in reg.groupby("region"):
    todo = [e for e in grp.pred_element if not os.path.exists(os.path.join(out_dir, f"{e}__ENFORMER.tsv.gz"))]
    if not todo:
        continue
    t0 = time.time()
    r = grp.iloc[0]
    start = (r.lo + r.hi) // 2 - SEQ // 2                       # 0-based start of the input window
    ref_seq = fetch(r.chrom, start)
    v = var[var.pred_element.isin(grp.pred_element)].drop_duplicates("variant_id")
    if a.limit:
        v = v.head(a.limit)
    off = (v.pos - 1 - start).to_numpy()
    bad = [i for i, (o, b) in enumerate(zip(off, v.ref)) if ref_seq[o] != b]
    assert not bad, f"{region}: reference mismatch at {len(bad)} variants"
    lo_bin = np.maximum((off - HALF - CROP) // BIN, 0)
    hi_bin = (off + HALF - CROP) // BIN
    ref = predict([ref_seq])[0]
    csum = np.vstack([np.zeros((1, ref.shape[1])), np.cumsum(ref, 0)])
    ref_sum = csum[hi_bin + 1] - csum[lo_bin]
    alt_sum = np.empty_like(ref_sum)
    alts = list(v.alt)
    for s in range(0, len(v), a.batch):
        seqs = [ref_seq[:o] + b + ref_seq[o + 1:] for o, b in zip(off[s:s + a.batch], alts[s:s + a.batch])]
        p = predict(seqs)
        for j in range(len(seqs)):
            alt_sum[s + j] = p[j, lo_bin[s + j]:hi_bin[s + j] + 1].sum(0)
    score = np.log2(1 + alt_sum) - np.log2(1 + ref_sum)
    df = pd.DataFrame(score, columns=tg.track_name.to_numpy(), index=v.variant_id.to_numpy()).stack().rename(
        "raw_score").reset_index().rename(columns={"level_0": "variant_id", "level_1": "track_name"})
    meta = tg.set_index("track_name")
    df["output_type"] = df.track_name.map(meta.output_type)
    df["biosample_name"] = df.track_name.map(meta.biosample_name)
    df["track_strand"], df["biosample_type"] = ".", ""
    for e in todo:
        path = os.path.join(out_dir, f"{e}__ENFORMER.tsv.gz")
        ids = set(var.variant_id[var.pred_element == e])
        df[df.variant_id.isin(ids)].to_csv(path + ".tmp", sep="\t", index=False, compression="gzip")
        os.replace(path + ".tmp", path)
    print(f"{region} {todo}: {len(v)} SNVs, {time.time() - t0:.0f}s", flush=True)
print("done")
