#!/bin/sh
# Kircher mouse-cell experiments (docs/internal/kircher_mouse_zh.md): frozen TrackPilot (proxy_select.py, v2
# configuration) on every experiment, then the pre-registered decision.
# usage (from the repository root): sh pilot/scripts/run_kircher_mouse.sh <alphagenome|enformer>
PY=${PY:-.venv/bin/python}
m=$1
if [ "$m" = alphagenome ]; then ver=ALL_FOLDS; pred=predictions/alphagenome_all; else ver=ENFORMER; pred=predictions/enformer_kircher_mouse; fi
mkdir -p pilot/results/kircher_mouse_logs
for e in RET UC88 TCF7L2 ZFAND3 ZRSh-13 ZRSh-13h2; do
  "$PY" pilot/scripts/proxy_select.py pilot --variants "data/kircher_mouse/proxy_variants_$m.tsv" --pred-dir "$pred" \
    --version "$ver" --common-eval --anchor-prior 0.95 --family-prior --controls --elements "$e" \
    --tag "kircher_mouse_$m/$e" > "pilot/results/kircher_mouse_logs/$m.$e.log" 2>&1 &
done
wait
"$PY" pilot/scripts/kircher_mouse_report.py "$m"
