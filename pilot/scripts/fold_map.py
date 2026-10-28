"""C1 robustness, step 1: which AlphaGenome checkpoint never saw each MPRA element during training.

AlphaGenome's FOLD_k checkpoints each hold out a validation and a test fold of the Borzoi genome split
(alphagenome.data.fold_intervals): FOLD_0 -> fold0/fold1, FOLD_1 -> fold3/fold4, FOLD_2 -> fold2/fold5,
FOLD_3 -> fold6/fold7; ALL_FOLDS (the default, distilled) saw every fold. For each element we list the
folds of the training windows that overlap its region and the checkpoint(s) that held all of them out.

usage: python fold_map.py <pilot_dir>
Output: <pilot_dir>/data/element_folds.tsv, plus the scoring commands grouped by checkpoint
"""
import os, sys
import pandas as pd
from alphagenome.data import fold_intervals as fi
from alphagenome.models.v1 import dna_client

pilot_dir = sys.argv[1]
var = pd.read_csv(os.path.join(pilot_dir, "data", "variants_grch38.tsv"), sep="\t", low_memory=False)
var = var[var.human_cells]
el = var.groupby("element").agg(chrom=("chrom", "first"), start=("pos", "min"), end=("pos", "max")).reset_index()

held_out = {}   # fold name -> checkpoint that never trained on it
for v in [dna_client.ModelVersion.FOLD_0, dna_client.ModelVersion.FOLD_1, dna_client.ModelVersion.FOLD_2,
          dna_client.ModelVersion.FOLD_3]:
    for f in fi.get_fold_names(v, fi.Subset.VALID) + fi.get_fold_names(v, fi.Subset.TEST):
        held_out[f] = v.name
# every fold's intervals (the ALL_FOLDS training set lists all eight)
iv = fi.get_fold_intervals(dna_client.ModelVersion.ALL_FOLDS, dna_client.Organism.HOMO_SAPIENS, fi.Subset.TRAIN)
cols = {c.lower(): c for c in iv.columns}
chrom_c, start_c, end_c = cols.get("chromosome", cols.get("chrom")), cols["start"], cols["end"]
fold_c = next(c for c in iv.columns if c.lower() in ("fold", "name", "subset", "split"))

rows = []
for r in el.itertuples():
    hit = iv[(iv[chrom_c] == r.chrom) & (iv[start_c] < r.end) & (iv[end_c] > r.start - 1)]
    folds = sorted(set(hit[fold_c]))
    versions = sorted({held_out[f] for f in folds})
    rows.append(dict(element=r.element, region=f"{r.chrom}:{r.start}-{r.end}", folds=",".join(folds),
                     held_out_checkpoint=versions[0] if len(versions) == 1 else "",
                     note="" if len(versions) == 1 else f"windows from several folds: {versions}"))
out = pd.DataFrame(rows)
out.to_csv(os.path.join(pilot_dir, "data", "element_folds.tsv"), sep="\t", index=False)
print(out.to_string(index=False))
print("\nscoring commands (held-out checkpoint per element):")
for v, g in out[out.held_out_checkpoint != ""].groupby("held_out_checkpoint"):
    print(f".venv/bin/python pilot/scripts/alphagenome_score.py pilot --all-biosamples --outputs CAGE,DNASE,ATAC "
          f"--versions {v} --elements {','.join(g.element)}")
