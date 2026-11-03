#!/bin/sh
# Neuron MPRA (docs/internal/neuron_mpra_zh.md): frozen TrackPilot (proxy_select.py, v2 configuration) on every region,
# then the pre-registered decision. One process per region (all tiles of the region).
# usage (from the repository root): sh pilot/scripts/run_kosicki.sh <alphagenome|enformer> [parallel jobs, default 8]
PY=${PY:-.venv/bin/python}
if [ "$1" = "--one" ]; then
  m=$2; r=$3
  if [ "$m" = alphagenome ]; then ver=ALL_FOLDS; else ver=ENFORMER; fi
  els=$(awk -F'\t' -v r="$r" 'NR > 1 && $2 == r {print $1}' "pilot/data/kosicki/proxy_variants_$m.tsv" | sort -u | paste -sd, -)
  "$PY" pilot/scripts/proxy_select.py pilot --variants "data/kosicki/proxy_variants_$m.tsv" \
    --pred-dir "predictions/kosicki_$m" --version "$ver" --common-eval --anchor-prior 0.95 --family-prior --controls \
    --elements "$els" --tag "kosicki_$m/$r" > "pilot/results/kosicki_logs/$m.$r.log" 2>&1 || echo "FAIL $m $r"
  exit 0
fi
set -e
m=$1
mkdir -p pilot/results/kosicki_logs
awk -F'\t' 'NR > 1 {print $2}' "pilot/data/kosicki/proxy_variants_$m.tsv" | sort -u | while read -r r; do echo "$m $r"; done \
  | xargs -P "${2:-8}" -n 2 sh "$0" --one
"$PY" pilot/scripts/kosicki_report.py "$m"
