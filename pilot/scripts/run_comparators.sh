#!/bin/sh
# Post hoc acquisition comparators (docs/internal/acquisition_comparators_zh.md): frozen v2 with --comparators on
# the first 50 repeats (the same pools and labels as the original runs), then the summary.
# usage (from the repository root):
#   sh pilot/scripts/run_comparators.sh dev   [jobs]   Kircher development elements; freezes the stopping threshold
#   sh pilot/scripts/run_comparators.sh test  [jobs]   Kircher test elements
#   sh pilot/scripts/run_comparators.sh siraj [jobs]   the 125 Siraj loci
#   SET=siraj_ref ... siraj: labels relative to the measured reference oligo (docs/internal/siraj_refnorm_zh.md)
PY=${PY:-.venv/bin/python}
LOG=pilot/results/cmp_logs
SET=${SET:-siraj}
CMP="--common-eval --anchor-prior 0.95 --family-prior --controls --comparators --repeats 50 --comparator-repeats 50"
if [ "$1" = "--kircher" ]; then
  "$PY" pilot/scripts/proxy_select.py pilot $CMP --scenarios all,hide_matched --elements "$2" --tag "cmp/$2" \
    > "$LOG/$2.log" 2>&1 || echo "FAIL $2"
  exit 0
fi
if [ "$1" = "--siraj" ]; then
  "$PY" pilot/scripts/proxy_select.py pilot --variants "data/$SET/proxy_variants.tsv" --pred-dir predictions/siraj_all \
    $CMP --elements "$2.K562,$2.HEPG2" --tag "${SET}_cmp/$2" > "$LOG/${SET}_$2.log" 2>&1 || echo "FAIL $2"
  exit 0
fi
set -e
mkdir -p "$LOG"
DEV="PKLR-48h F9 HNF4A IRF6 BCL11A"
TEST="FOXE1 GP1BA HBB HBG1 IRF4 LDLR LDLR.2 MSMB MYCrs11986220 MYCrs6983267 SORT1 SORT1-flip SORT1.2 TERT-GAa TERT-GBM TERT-GSc TERT-HEK"
case "$1" in
  dev)   printf '%s\n' $DEV | xargs -P "${2:-8}" -n 1 sh "$0" --kircher
         "$PY" pilot/scripts/cmp_report.py kircher pilot pilot/results/proxy_select_cmp results/a1/cmp dev ;;
  test)  printf '%s\n' $TEST | xargs -P "${2:-8}" -n 1 sh "$0" --kircher
         "$PY" pilot/scripts/cmp_report.py kircher pilot pilot/results/proxy_select_cmp results/a1/cmp test ;;
  siraj) ls pilot/predictions/siraj_all/*__ALL_FOLDS.tsv.gz | xargs -n 1 basename | sed 's/__ALL_FOLDS.tsv.gz//' \
           | xargs -P "${2:-8}" -n 1 sh "$0" --siraj
         if [ "$SET" = siraj ]; then o=results/a1/cmp; else o=results/a1/$SET/cmp; fi
         "$PY" pilot/scripts/cmp_report.py siraj pilot "pilot/results/proxy_select_${SET}_cmp" "$o" ;;
  *)     echo "usage: sh pilot/scripts/run_comparators.sh dev|test|siraj [jobs]"; exit 1 ;;
esac
