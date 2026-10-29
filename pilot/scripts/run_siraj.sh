#!/bin/sh
# Pre-registered Siraj SatMut run of the frozen method v2 (docs/internal/preregistration_v2_siraj_zh.md, section 7):
# one process per scored locus, the main configuration and the uniform-prior ablation (S5), then the decision.
# usage (from the repository root): sh pilot/scripts/run_siraj.sh [parallel jobs, default 8]
# C1 on Siraj (docs/internal/c1_robustness_zh.md): sh pilot/scripts/run_siraj.sh [jobs] fold
#   reruns the main configuration with each locus's held-out-fold checkpoint (pilot/data/siraj/element_folds.tsv)
#   and writes the same decision tables to results/a1/c1_siraj/.
PY=${PY:-.venv/bin/python}
LOG=pilot/results/siraj_logs
MAP=pilot/data/siraj/element_folds.tsv
if [ "$1" = "--one" ]; then
  l=$2
  "$PY" pilot/scripts/proxy_select.py pilot --variants data/siraj/proxy_variants.tsv --pred-dir predictions/siraj_all \
    --common-eval --anchor-prior 0.95 --family-prior --controls --elements "$l.K562,$l.HEPG2" --tag "siraj/$l" \
    > "$LOG/$l.log" 2>&1 || echo "FAIL main $l"
  "$PY" pilot/scripts/proxy_select.py pilot --variants data/siraj/proxy_variants.tsv --pred-dir predictions/siraj_all \
    --common-eval --anchor-prior 0.95 --scenarios all --elements "$l.K562,$l.HEPG2" --tag "siraj_uni/$l" \
    > "$LOG/$l.uni.log" 2>&1 || echo "FAIL uni $l"
  exit 0
fi
if [ "$1" = "--fold" ]; then
  l=$2
  ck=$(awk -F'\t' -v l="$l" 'NR > 1 && $1 == l {print $4}' "$MAP")
  "$PY" pilot/scripts/proxy_select.py pilot --variants data/siraj/proxy_variants.tsv --pred-dir predictions/siraj_all \
    --version "$ck" --common-eval --anchor-prior 0.95 --family-prior --controls --elements "$l.K562,$l.HEPG2" \
    --tag "siraj_fold/$l" > "$LOG/$l.fold.log" 2>&1 || echo "FAIL fold $l"
  exit 0
fi
set -e
mkdir -p "$LOG"
if [ "$2" = "fold" ]; then
  awk -F'\t' 'NR > 1 {print $1 "__" $4 ".tsv.gz"}' "$MAP" | while read -r f; do
    [ -f "pilot/predictions/siraj_all/$f" ] && echo "${f%%__*}"
  done > "$LOG/loci_fold.txt"
  echo "$(wc -l < "$LOG/loci_fold.txt" | tr -d ' ') loci with a held-out-fold prediction"
  xargs -P "${1:-8}" -n 1 sh "$0" --fold < "$LOG/loci_fold.txt"
  "$PY" pilot/scripts/siraj_report.py pilot/results/proxy_select_siraj_fold results/a1/c1_siraj
  exit 0
fi
ls pilot/predictions/siraj_all/*__ALL_FOLDS.tsv.gz | xargs -n 1 basename | sed 's/__ALL_FOLDS.tsv.gz//' > "$LOG/loci.txt"
echo "$(wc -l < "$LOG/loci.txt" | tr -d ' ') scored loci"
xargs -P "${1:-8}" -n 1 sh "$0" --one < "$LOG/loci.txt"
"$PY" pilot/scripts/siraj_report.py pilot/results/proxy_select_siraj results/a1
