#!/bin/sh
# Pre-registered Siraj SatMut run of the frozen method v2 (docs/internal/preregistration_v2_siraj_zh.md, section 7):
# one process per scored locus, the main configuration and the uniform-prior ablation (S5), then the decision.
# usage (from the repository root): sh pilot/scripts/run_siraj.sh [parallel jobs, default 8]
PY=${PY:-.venv/bin/python}
LOG=pilot/results/siraj_logs
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
set -e
mkdir -p "$LOG"
ls pilot/predictions/siraj_all/*__ALL_FOLDS.tsv.gz | xargs -n 1 basename | sed 's/__ALL_FOLDS.tsv.gz//' > "$LOG/loci.txt"
echo "$(wc -l < "$LOG/loci.txt" | tr -d ' ') scored loci"
xargs -P "${1:-8}" -n 1 sh "$0" --one < "$LOG/loci.txt"
"$PY" pilot/scripts/siraj_report.py pilot/results/proxy_select_siraj results/a1
