#!/bin/sh
# C1 robustness (docs/internal/c1_robustness_zh.md): rerun v1 and frozen v2 on each Kircher element with the
# checkpoint that held its genomic fold out, next to ALL_FOLDS, then summarise.
# usage (from the repository root): sh pilot/scripts/run_c1.sh [parallel jobs, default 8]
PY=${PY:-.venv/bin/python}
LOG=pilot/results/c1_logs
V2="--anchor-prior 0.95 --family-prior --controls --scenarios all,hide_matched"
if [ "$1" = "--one" ]; then
  cfg=$2; ver=$3; el=$4
  if [ "$cfg" = v1 ]; then flags=""; else flags=$V2; fi
  if [ "$ver" = ALL_FOLDS ]; then side=all; else side=fold; fi
  "$PY" pilot/scripts/proxy_select.py pilot --common-eval --version "$ver" --elements "$el" $flags \
    --tag "c1/${cfg}_${side}/$el" > "$LOG/${cfg}_${side}_$el.log" 2>&1 || echo "FAIL $cfg $ver $el"
  exit 0
fi
set -e
mkdir -p "$LOG"
: > "$LOG/jobs.txt"
for pair in FOLD_0:F9,MSMB,SORT1,SORT1-flip,SORT1.2 FOLD_1:BCL11A,FOXE1,MYCrs11986220,MYCrs6983267 \
    FOLD_2:GP1BA,HBB,HBG1,IRF4,PKLR-48h,TERT-GAa,TERT-GBM,TERT-GSc,TERT-HEK FOLD_3:HNF4A,IRF6,LDLR,LDLR.2; do
  fold=${pair%%:*}
  for el in $(echo "${pair#*:}" | tr ',' ' '); do
    for cfg in v1 v2; do
      echo "$cfg $fold $el" >> "$LOG/jobs.txt"
      echo "$cfg ALL_FOLDS $el" >> "$LOG/jobs.txt"
    done
  done
done
echo "$(wc -l < "$LOG/jobs.txt" | tr -d ' ') jobs"
xargs -P "${1:-8}" -L 1 sh "$0" --one < "$LOG/jobs.txt"
"$PY" pilot/scripts/c1_report.py pilot results/a1
