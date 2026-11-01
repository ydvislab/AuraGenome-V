#!/bin/sh
# Post hoc EIG convergence check (docs/internal/eig_convergence_zh.md): frozen v2 plus --eig-check on the first 50
# repeats (the same pools as the original comparator runs), then the summary.
# usage (from the repository root):
#   sh pilot/scripts/run_eig_check.sh test  [jobs]   Kircher test elements
#   sh pilot/scripts/run_eig_check.sh siraj [jobs]   the 125 Siraj loci, labels relative to the reference oligo
PY=${PY:-.venv/bin/python}
LOG=pilot/results/eigconv_logs
V2="--common-eval --anchor-prior 0.95 --family-prior --controls --repeats 50 --comparator-repeats 50"
EIG="--eig-check 100:96,300:96,0.99:192"
if [ "$1" = "--kircher" ]; then
  "$PY" pilot/scripts/proxy_select.py pilot $V2 $EIG --scenarios all,hide_matched --elements "$2" \
    --tag "eigconv/$2" > "$LOG/$2.log" 2>&1 || echo "FAIL $2"
  exit 0
fi
if [ "$1" = "--siraj" ]; then
  "$PY" pilot/scripts/proxy_select.py pilot --variants data/siraj_ref/proxy_variants.tsv --pred-dir predictions/siraj_all \
    $V2 $EIG --elements "$2.K562,$2.HEPG2" --tag "siraj_ref_eigconv/$2" > "$LOG/siraj_$2.log" 2>&1 || echo "FAIL $2"
  exit 0
fi
set -e
mkdir -p "$LOG"
TEST="FOXE1 GP1BA HBB HBG1 IRF4 LDLR LDLR.2 MSMB MYCrs11986220 MYCrs6983267 SORT1 SORT1-flip SORT1.2 TERT-GAa TERT-GBM TERT-GSc TERT-HEK"
case "$1" in
  test)  printf '%s\n' $TEST | xargs -P "${2:-12}" -n 1 sh "$0" --kircher
         "$PY" pilot/scripts/eig_check_report.py kircher pilot pilot/results/proxy_select_eigconv results/a1/eigconv ;;
  siraj) ls pilot/predictions/siraj_all/*__ALL_FOLDS.tsv.gz | xargs -n 1 basename | sed 's/__ALL_FOLDS.tsv.gz//' \
           | xargs -P "${2:-12}" -n 1 sh "$0" --siraj
         "$PY" pilot/scripts/eig_check_report.py siraj pilot pilot/results/proxy_select_siraj_ref_eigconv results/a1/eigconv ;;
  *)     echo "usage: sh pilot/scripts/run_eig_check.sh test|siraj [jobs]"; exit 1 ;;
esac
