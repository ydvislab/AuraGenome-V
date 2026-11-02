"""Optional: score every substitution of a region with AlphaGenome by in-silico mutagenesis.

Needs `pip install trackpilot[alphagenome]` and an API key in the environment variable ALPHAGENOME_API_KEY (the
key is read from the environment only; AlphaGenome's terms of use apply). The calls are those of
pilot/scripts/alphagenome_score.py, which produced the predictions in the paper: 1-Mb context centred on the region,
the recommended CAGE, DNase and ATAC scorers, all human biosamples, default (distilled) checkpoint.
"""
from __future__ import annotations

import os
import time

import pandas as pd


def score_region(region, outputs=("CAGE", "DNASE", "ATAC"), chunk_bp=100, model_version="ALL_FOLDS",
                 seq_len="SEQUENCE_LENGTH_1MB"):
    """Tidy scores (one row per variant x track) for every substitution in `region` ('chrom:start-end', 1-based,
    inclusive, GRCh38)."""
    from alphagenome.data import genome
    from alphagenome.models.v1 import dna_client, variant_scorers

    key = os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        raise SystemExit("set ALPHAGENOME_API_KEY (see https://github.com/google-deepmind/alphagenome)")
    chrom, span = region.split(":")
    start, end = (int(x.replace(",", "")) for x in span.split("-"))
    scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[k] for k in outputs]
    client = dna_client.create(key, model_version=dna_client.ModelVersion[model_version])
    interval = genome.Interval(chromosome=chrom, start=start - 1, end=end).resize(
        dna_client.SUPPORTED_SEQUENCE_LENGTHS[seq_len])
    parts = []
    for s in range(start, end + 1, chunk_bp):
        e = min(s + chunk_bp - 1, end)
        t0 = time.time()
        scores = client.score_ism_variants(interval=interval, ism_interval=genome.Interval(chromosome=chrom,
                                           start=s - 1, end=e), variant_scorers=scorers,
                                           organism=dna_client.Organism.HOMO_SAPIENS, progress_bar=False)
        # score_ism_variants does not record which scorer produced each AnnData (alphagenome 0.10.0): attach them
        # by position and check that each position keeps one track set
        names_at = {}
        for per_variant in scores:
            assert len(per_variant) == len(scorers)
            for i, ad in enumerate(per_variant):
                ad.uns["variant_scorer"] = scorers[i]
                names = tuple(ad.var["name"]) if "name" in ad.var else tuple(ad.var.index)
                assert names_at.setdefault(i, names) == names, f"scorer {i}: track set differs between variants"
        parts.append(variant_scorers.tidy_scores(scores))
        print(f"  {chrom}:{s}-{e}: {time.time() - t0:.0f}s", flush=True)
    df = pd.concat(parts, ignore_index=True)
    df["variant_id"] = df.variant_id.astype(str)          # e.g. chr11:5227002:G>A
    return df
