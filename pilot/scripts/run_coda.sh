#!/bin/sh
# Post hoc CODA / MODEL SELECTOR style baseline (docs/internal/coda_baseline_zh.md): frozen v2 plus --coda on the
# first 50 repeats (the same pools as the original runs), then the summary.
# usage (from the repository root):
#   sh pilot/scripts/run_coda.sh dev   [jobs]           Kircher development elements, kappa grid 0,3,10,30; chooses kappa
#   sh pilot/scripts/run_coda.sh test  <kappa> [jobs]   Kircher test elements, kappa 0 and the chosen kappa
#   sh pilot/scripts/run_coda.sh siraj <kappa> [jobs]   the 125 Siraj loci, labels relative to the reference oligo
PY=${PY:-.venv/bin/python}
LOG=pilot/results/coda_logs
V2="--common-eval --anchor-prior 0.95 --family-prior --controls --coda --repeats 50 --comparator-repeats 50"
if [ "$1" = "--kircher" ]; then
  "$PY" pilot/scripts/proxy_select.py pilot $V2 --coda-kappa "$KAPPA" --scenarios all,hide_matched --elements "$2" \
    --tag "coda/$2" > "$LOG/$2.log" 2>&1 || echo "FAIL $2"
  exit 0
fi
if [ "$1" = "--siraj" ]; then
  "$PY" pilot/scripts/proxy_select.py pilot --variants data/siraj_ref/proxy_variants.tsv --pred-dir predictions/siraj_all \
    $V2 --coda-kappa "$KAPPA" --elements "$2.K562,$2.HEPG2" --tag "siraj_ref_coda/$2" > "$LOG/siraj_$2.log" 2>&1 \
    || echo "FAIL $2"
  exit 0
fi
set -e
mkdir -p "$LOG"
DEV="PKLR-48h F9 HNF4A IRF6 BCL11A"
TEST="FOXE1 GP1BA HBB HBG1 IRF4 LDLR LDLR.2 MSMB MYCrs11986220 MYCrs6983267 SORT1 SORT1-flip SORT1.2 TERT-GAa TERT-GBM TERT-GSc TERT-HEK"
case "$1" in
  dev)   export KAPPA=0,3,10,30
         printf '%s\n' $DEV | xargs -P "${2:-8}" -n 1 sh "$0" --kircher
         "$PY" pilot/scripts/coda_report.py dev pilot/results/proxy_select_coda results/a1/coda ;;
  test)  export KAPPA=0,$2
         printf '%s\n' $TEST | xargs -P "${3:-8}" -n 1 sh "$0" --kircher
         "$PY" pilot/scripts/coda_report.py test pilot/results/proxy_select_coda results/a1/coda "$2" ;;
  siraj) export KAPPA=0,$2
         ls pilot/predictions/siraj_all/*__ALL_FOLDS.tsv.gz | xargs -n 1 basename | sed 's/__ALL_FOLDS.tsv.gz//' \
           | xargs -P "${3:-8}" -n 1 sh "$0" --siraj
         "$PY" pilot/scripts/coda_report.py siraj pilot/results/proxy_select_siraj_ref_coda results/a1/coda "$2" ;;
  *)     echo "usage: sh pilot/scripts/run_coda.sh dev [jobs] | test <kappa> [jobs] | siraj <kappa> [jobs]"; exit 1 ;;
esac
