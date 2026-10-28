"""Pilot step 2a: AlphaGenome variant-effect scores for every MPRA SNV, from five model versions.

Model versions (alphagenome.models.v1.dna_client.ModelVersion): ALL_FOLDS (distilled, default)
and FOLD_0..FOLD_3 (cross-validation checkpoints). Disagreement between versions is the
replicate-disagreement signal of the reliability profile (research plan §4.2).

Strategy: one in-silico-mutagenesis call per element region and model version
(DnaClient.score_ism_variants scores all three alternative bases at every position), instead
of one call per variant. Only tracks whose biosample matches the element's keywords are kept,
plus every track of the requested output types for a small set of always-kept biosamples
(--keep-biosamples), so outputs stay small.

Requires: Python >= 3.10, `pip install alphagenome` (written against alphagenome 0.10.0 source;
NOT yet executed - run with --dry-run first), and an API key in the environment variable
ALPHAGENOME_API_KEY (non-commercial terms of use apply).

usage:
  python alphagenome_score.py <pilot_dir> --dry-run
  python alphagenome_score.py <pilot_dir> --elements PKLR-48h,F9 --versions ALL_FOLDS,FOLD_0
  python alphagenome_score.py <pilot_dir>                     # everything, resumable
Output: <pilot_dir>/predictions/alphagenome/<element>__<version>.tsv.gz (tidy scores, filtered)
"""
import os, sys, argparse, time
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("pilot_dir")
ap.add_argument("--elements", default="", help="comma-separated element names; default all human-cell elements")
ap.add_argument("--versions", default="ALL_FOLDS,FOLD_0,FOLD_1,FOLD_2,FOLD_3")
ap.add_argument("--outputs", default="CAGE,RNA_SEQ,DNASE,ATAC,CHIP_HISTONE",
                help="keys of variant_scorers.RECOMMENDED_VARIANT_SCORERS")
ap.add_argument("--seq-len", default="SEQUENCE_LENGTH_1MB", help="key of dna_client.SUPPORTED_SEQUENCE_LENGTHS")
ap.add_argument("--keep-biosamples", default="K562;HepG2", help="always keep these biosamples (cross-cell-type contrast)")
ap.add_argument("--include-mouse-cell-elements", action="store_true")
ap.add_argument("--chunk-bp", type=int, default=100, help="ISM positions per API call (memory bound)")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()

var = pd.read_csv(os.path.join(a.pilot_dir, "data", "variants_grch38.tsv"), sep="\t", low_memory=False)
if not a.include_mouse_cell_elements:
    var = var[var.human_cells]
if a.elements:
    var = var[var.element.isin(a.elements.split(","))]

# element regions: identical genomic regions measured in several experiments (e.g. PKLR-24h/48h)
# are scored once and shared
regions = (var.groupby("element").agg(chrom=("chrom", "first"), start=("pos", "min"), end=("pos", "max"),
                                      keywords=("model_biosample_keywords", "first"), n=("variant_id", "size"))
           .reset_index())
regions["region"] = regions.chrom + ":" + regions.start.astype(str) + "-" + regions.end.astype(str)
versions = a.versions.split(",")
out_dir = os.path.join(a.pilot_dir, "predictions", "alphagenome")
os.makedirs(out_dir, exist_ok=True)

print(f"elements {len(regions)}, unique regions {regions.region.nunique()}, SNVs {len(var)}, "
      f"model versions {len(versions)} -> ISM calls {regions.region.nunique() * len(versions)}")
print(regions[["element", "region", "n", "keywords"]].to_string(index=False))
if a.dry_run:
    sys.exit(0)

from alphagenome.data import genome  # noqa: E402
from alphagenome.models.v1 import dna_client, variant_scorers  # noqa: E402

key = os.environ.get("ALPHAGENOME_API_KEY")
if not key:
    sys.exit("set ALPHAGENOME_API_KEY (see https://github.com/google-deepmind/alphagenome)")
scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[k] for k in a.outputs.split(",")]
seq_len = dna_client.SUPPORTED_SEQUENCE_LENGTHS[a.seq_len]
always = [k for k in a.keep_biosamples.split(";") if k]

def keep(df, keywords):
    kws = [k for k in keywords.split(";") if k] + always
    name = df.get("biosample_name", pd.Series("", index=df.index)).fillna("").astype(str)
    track = df["track_name"].fillna("").astype(str)
    hit = pd.Series(False, index=df.index)
    for k in kws:
        hit |= name.str.contains(k, case=False, regex=False) | track.str.contains(k, case=False, regex=False)
    return df[hit]

def attach_scorers(scores):
    """alphagenome 0.10.0: score_ism_variants returns one AnnData per scorer for each variant but,
    unlike score_variant, does not record which scorer produced it, so tidy_scores rejects it.
    Attach scorers by position (the order score_variant itself relies on) and check that every
    position carries the same track set across variants, so a different order would fail loudly."""
    tracks_at = {}
    for per_variant in scores:
        assert len(per_variant) == len(scorers), f"expected {len(scorers)} scores per variant, got {len(per_variant)}"
        for k, ad in enumerate(per_variant):
            ad.uns["variant_scorer"] = scorers[k]
            names = tuple(ad.var["name"]) if "name" in ad.var else tuple(ad.var.index)
            assert tracks_at.setdefault(k, names) == names, f"scorer {k}: track set differs between variants"

def score_region(client, chrom, start, end, keywords):
    """ISM over [start, end] (1-based, inclusive) in chunks, keeping matched tracks only.
    The prediction window is centred on the whole region so all chunks share one context."""
    region = genome.Interval(chromosome=chrom, start=start - 1, end=end)  # 0-based, end-exclusive
    interval = region.resize(seq_len)
    parts = []
    for s in range(start, end + 1, a.chunk_bp):
        e = min(s + a.chunk_bp - 1, end)
        ism = genome.Interval(chromosome=chrom, start=s - 1, end=e)
        t0 = time.time()
        scores = client.score_ism_variants(interval=interval, ism_interval=ism, variant_scorers=scorers,
                                           organism=dna_client.Organism.HOMO_SAPIENS, progress_bar=False)
        attach_scorers(scores)
        df = keep(variant_scorers.tidy_scores(scores), keywords)
        print(f"  {chrom}:{s}-{e}: {len(df)} kept rows, {time.time() - t0:.0f}s", flush=True)
        parts.append(df)
    return pd.concat(parts, ignore_index=True)

# experiments that measured the same region (e.g. PKLR-24h / PKLR-48h) share one scoring run,
# filtered with the union of their biosample keywords and then split per experiment
for version in versions:
    client = dna_client.create(key, model_version=dna_client.ModelVersion[version])
    for region, grp in regions.groupby("region"):
        todo = [r for r in grp.itertuples() if not os.path.exists(os.path.join(out_dir, f"{r.element}__{version}.tsv.gz"))]
        if not todo:
            continue
        union = ";".join(sorted({k for r in todo for k in r.keywords.split(";") if k}))
        print(f"{region} {version} -> {[r.element for r in todo]}", flush=True)
        df = score_region(client, todo[0].chrom, int(todo[0].start), int(todo[0].end), union)
        for r in todo:
            keep(df, r.keywords).to_csv(os.path.join(out_dir, f"{r.element}__{version}.tsv.gz"), sep="\t", index=False)
print("done")
