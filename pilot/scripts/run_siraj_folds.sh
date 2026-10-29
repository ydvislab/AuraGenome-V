#!/bin/sh
# Score each Siraj locus with the AlphaGenome checkpoint that held its genomic fold out (C1 on Siraj).
# Fold map: pilot/data/siraj/element_folds.tsv from fold_map.py. Needs ALPHAGENOME_API_KEY in this shell.
# Resumable: finished loci are skipped. Outputs <locus>__FOLD_k.tsv.gz next to the ALL_FOLDS files.
# usage (from the repository root): sh pilot/scripts/run_siraj_folds.sh [FOLD_0 FOLD_1 ...]   (default all four)
PY=${PY:-.venv/bin/python}
MAP=pilot/data/siraj/element_folds.tsv
[ $# -gt 0 ] || set -- FOLD_0 FOLD_1 FOLD_2 FOLD_3
for ck in "$@"; do
  els=$(awk -F'\t' -v ck="$ck" 'NR > 1 && $4 == ck {printf "%s%s", sep, $1; sep = ","}' "$MAP")
  [ -n "$els" ] || { echo "no loci for $ck"; continue; }
  echo "== $ck: $(echo "$els" | tr ',' '\n' | wc -l | tr -d ' ') loci"
  "$PY" pilot/scripts/alphagenome_score.py pilot --variants data/siraj/variants_grch38.tsv --out predictions/siraj_all \
    --versions "$ck" --outputs CAGE,DNASE,ATAC --all-biosamples --elements "$els" || exit 1
done
