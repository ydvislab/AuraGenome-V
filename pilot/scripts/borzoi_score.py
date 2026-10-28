"""Pilot step 2b: Borzoi variant-effect scores from the four released replicates (rep0-rep3).

Second model family for cross-model disagreement, and a second replicate-disagreement
estimate (research plan §4.2). For every SNV: predict reference and alternative sequence in a
524,288 bp window centred on the element, then take, for each matched track, the
log2((alt_sum + 1) / (ref_sum + 1)) over the output bins covering a 501 bp window around the
variant - the same idea as AlphaGenome's CenterMaskScorer(width=501, DIFF_LOG2_SUM), so the two
models' scores are comparable in direction.

Requires: Python >= 3.10, a CUDA GPU (CPU works but is very slow), `pip install borzoi-pytorch pyfaidx`,
hg38 FASTA (UCSC hg38.fa + .fai), and Borzoi's human target table
(https://raw.githubusercontent.com/calico/borzoi/main/examples/targets_human.txt).
Weights: Hugging Face johahi/borzoi-replicate-{0..3} (ported with permission; check the licence
of the port before redistribution).
NOT yet executed: the borzoi-pytorch call pattern follows the package README and its example
notebook (input (batch, 4, 524288) one-hot ACGT; output (batch, tracks, bins)). Verify the
output bin size and cropping (`--bin-size`, `--crop-bins`) on the first run with --check.

usage:
  python borzoi_score.py <pilot_dir> --fasta hg38.fa --targets targets_human.txt --check
  python borzoi_score.py <pilot_dir> --fasta hg38.fa --targets targets_human.txt --elements PKLR-48h --replicates 0,1,2,3
Output: <pilot_dir>/predictions/borzoi/<element>__rep<k>.tsv.gz
        columns: variant_id, track_index, identifier, description, raw_score
"""
import os, sys, argparse, time
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("pilot_dir")
ap.add_argument("--fasta", required=True)
ap.add_argument("--targets", required=True)
ap.add_argument("--elements", default="")
ap.add_argument("--replicates", default="0,1,2,3")
ap.add_argument("--seq-len", type=int, default=524288)
ap.add_argument("--bin-size", type=int, default=32)
ap.add_argument("--crop-bins", type=int, default=5120, help="bins cropped from each side before the output; verify with --check")
ap.add_argument("--width", type=int, default=501)
ap.add_argument("--batch", type=int, default=4)
ap.add_argument("--check", action="store_true", help="run one reference prediction and print output shape, then exit")
a = ap.parse_args()

import torch  # noqa: E402
from pyfaidx import Fasta  # noqa: E402
from borzoi_pytorch import Borzoi  # noqa: E402

dev = "cuda" if torch.cuda.is_available() else "cpu"
var = pd.read_csv(os.path.join(a.pilot_dir, "data", "variants_grch38.tsv"), sep="\t", low_memory=False)
var = var[var.human_cells]
if a.elements:
    var = var[var.element.isin(a.elements.split(","))]
tg = pd.read_csv(a.targets, sep="\t", index_col=0)
fa = Fasta(a.fasta)
out_dir = os.path.join(a.pilot_dir, "predictions", "borzoi")
os.makedirs(out_dir, exist_ok=True)
LUT = {"A": 0, "C": 1, "G": 2, "T": 3}

def one_hot(seq):
    x = np.zeros((4, len(seq)), np.float32)
    for i, b in enumerate(seq.upper()):
        if b in LUT:
            x[LUT[b], i] = 1.0
    return x

def matched_tracks(keywords):
    kws = [k for k in keywords.split(";") if k]
    m = np.zeros(len(tg), bool)
    for k in kws:
        m |= tg.description.str.contains(k, case=False, regex=False).values
    return np.where(m)[0]

def window(chrom, centre):
    start = centre - a.seq_len // 2
    return start, str(fa[chrom][start:start + a.seq_len])

for rep in [int(r) for r in a.replicates.split(",")]:
    model = Borzoi.from_pretrained(f"johahi/borzoi-replicate-{rep}").to(dev).eval()
    for element, ev in var.groupby("element"):
        out = os.path.join(out_dir, f"{element}__rep{rep}.tsv.gz")
        if os.path.exists(out) and not a.check:
            continue
        tracks = matched_tracks(ev.model_biosample_keywords.iloc[0])
        chrom = ev.chrom.iloc[0]
        centre = int((ev.pos.min() + ev.pos.max()) // 2)
        wstart, ref_seq = window(chrom, centre)
        ref_x = one_hot(ref_seq)
        with torch.no_grad():
            ref_pred = model(torch.as_tensor(ref_x)[None].to(dev)).float().cpu().numpy()[0]  # (tracks, bins)
        if a.check:
            print(f"element {element}: input {ref_x.shape}, output {ref_pred.shape}, matched tracks {len(tracks)}")
            print(tg.iloc[tracks[:10]].description.to_string())
            sys.exit(0)
        out_start = wstart + a.crop_bins * a.bin_size  # genomic start of the first output bin
        rows, t0 = [], time.time()
        for i in range(0, len(ev), a.batch):
            chunk = ev.iloc[i:i + a.batch]
            xs = []
            for v in chunk.itertuples():
                off = v.pos - 1 - wstart
                assert ref_seq[off].upper() == v.ref, f"reference mismatch at {v.variant_id}"
                x = ref_x.copy(); x[:, off] = 0; x[LUT[v.alt], off] = 1; xs.append(x)
            with torch.no_grad():
                alt_pred = model(torch.as_tensor(np.stack(xs)).to(dev)).float().cpu().numpy()
            for v, ap_ in zip(chunk.itertuples(), alt_pred):
                b0 = max(0, (v.pos - 1 - a.width // 2 - out_start) // a.bin_size)
                b1 = (v.pos - 1 + a.width // 2 - out_start) // a.bin_size + 1
                ref_sum = ref_pred[tracks, b0:b1].sum(1); alt_sum = ap_[tracks, b0:b1].sum(1)
                score = np.log2(alt_sum + 1) - np.log2(ref_sum + 1)
                rows += [(v.variant_id, int(t), tg.identifier.iloc[t], tg.description.iloc[t], float(s)) for t, s in zip(tracks, score)]
        pd.DataFrame(rows, columns=["variant_id", "track_index", "identifier", "description", "raw_score"]).to_csv(out, sep="\t", index=False)
        print(f"rep{rep} {element}: {len(ev)} variants x {len(tracks)} tracks in {time.time() - t0:.0f}s", flush=True)
print("done")
