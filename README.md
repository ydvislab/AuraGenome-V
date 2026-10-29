# AuraGenome-V

**Choosing which foundation-model track to trust for a cell type the model does not profile, from a handful of measurements.**

Genomic foundation models such as AlphaGenome predict the effect of any variant on thousands of functional
tracks (cell types × assays). When the cell type of interest is not among them, a researcher has to pick a proxy
track, and that choice can flip conclusions. AuraGenome-V (working name) treats every track as a hypothesis,
chooses the few variants whose measurement best discriminates between the tracks (sequential Bayesian model
discrimination), and returns a posterior-weighted prediction with a stopping rule; a visual-analytics interface
will show how the posterior spreads over track clusters. Validation uses saturation-mutagenesis MPRA
(Kircher et al. 2019) with labels masked.

> **Pre-registered Siraj test (2026-10-29): all four primary hypotheses hold; the safety-net claim is supported.**
> Frozen v2 on 125 new loci of the Siraj 2026 saturation MPRA (K562 and HepG2, pre-registered in commit 239d903):
> - with matched tracks kept, 8 labels are non-inferior to the matched mean (+0.010, one-sided lower bound +0.005;
>   3/125 loci worse by more than 0.05);
> - with own and related cell tracks hidden, 8 labels beat the all-track mean (+0.048, p = 8e-13);
> - with only own cell tracks hidden, 16 labels beat the related-cell mean (+0.041, p = 1e-19).
>
> The discrimination panel is essential: a random panel falls below the label-free mean. Results are in
> `results/a1/siraj_*.tsv` and section 8 of docs/internal/preregistration_v2_siraj_zh.md.

> **Method v2 (2026-10-29, docs/internal/method_v2_zh.md):** the matched-track mean enters as a hypothesis with
> prior mass 0.95 and a positive slope prior (chosen on the development set by a rule declared before the runs, then
> frozen and pushed before the test run). On the Kircher test regions, with matched tracks kept, v2 ties the matched
> mean at k = 8 (+0.001; 1/11 regions harmed vs 7/11 for v1) and is ahead at k = 16 (+0.018, none harmed); rescue is
> slower (LDLR 0.25 → 0.29 at k = 8, 0.41 at k = 16). Its inference ties a metadata-centred ridge on the same panel, so
> the increment comes from the panel choice; 8-label affine recalibration of the matched mean is harmful (sign flips).

> **Status (2026-10-29): pilot passed; method v2 needed.** On the pre-registered held-out test (11 regions, 17 MPRA
> experiments, tracks of the element's own and related cell types hidden), 8 discrimination-chosen measurements beat
> random panels (+0.08 Spearman, Wilcoxon p = 0.003) and the label-free all-track mean (+0.07, 8/11 regions); the
> advantage over picking variants with the largest predicted effect is not significant (+0.03, p = 0.14). With the
> same 8 labels it also beats LASSO over all tracks, the protocol AlphaGenome uses on this benchmark (+0.11 to +0.18
> depending on how degenerate LASSO fits are scored, p = 0.002; added after the pre-registration).
> **Limitation found after the test:** when metadata-matched or related-cell tracks are kept, simply averaging them
> is as good as the method (method − matched mean −0.03, 3/11 regions positive); the method helps where the matched
> proxy is poor (e.g. LDLR 0.24 → 0.51). The next version puts metadata into the prior and will be tested on a
> second dataset. Re-scoring every method and reference on one common evaluation set (variants outside the pool)
> leaves all conclusions unchanged (`results/a1/common_eval_*`). Cross-model novelty check: PROCEED WITH CAUTION
> (7/10); cross-model review: PROCEED WITH CAUTION, readiness 5/10.
> Collaborators: start with **[docs/STATUS_zh.md](docs/STATUS_zh.md)** (project overview, in Chinese) and the
> figures in [results/a1/](results/a1/).

![Held-out test learning curve](results/a1/fig3_test_learning_curve.png)

## Repository layout

```
pilot/
  scripts/
    download_data.py        ground truth: Kircher et al. 2019 saturation-mutagenesis MPRA (pinned, checksummed)
    build_variant_table.py  SNV table (GRCh38) with measured effects and ground-truth classes
    alphagenome_score.py    AlphaGenome in-silico mutagenesis per element (--all-biosamples keeps every track)
    proxy_select.py         A1 method (frozen) and baselines: panel strategies x budgets x repeats
    prereg_test.py          pre-registered held-out test and decision rule
    a1_report.py            collects tracked tables and figures into results/a1/
    compute_signals.py      superseded "reliability profile" pilot (kept for reference)
    borzoi_score.py         Borzoi replicate scores (GPU; planned robustness analysis)
  data/, predictions/, results/   downloaded and derived data (not tracked)
results/a1/                 tracked summary tables and figures
docs/STATUS_zh.md           project overview for collaborators (Chinese)
docs/internal/              pre-registration, pilot notes, idea-selection report; older plans marked superseded
```

## Reproduce

Requires Python ≥ 3.10 and an AlphaGenome API key (non-commercial terms; each person uses their own key,
exported in the shell, never committed).

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python pilot/scripts/download_data.py pilot
python pilot/scripts/build_variant_table.py pilot
export ALPHAGENOME_API_KEY=...
# development set (5 elements, ~25 min of API time)
python pilot/scripts/alphagenome_score.py pilot --all-biosamples --outputs CAGE,DNASE,ATAC --versions ALL_FOLDS --elements PKLR-48h,F9,HNF4A,IRF6,BCL11A
python pilot/scripts/proxy_select.py pilot
# ablations
python pilot/scripts/proxy_select.py pilot --noise --tag noise
python pilot/scripts/proxy_select.py pilot --ensemble --tag ensemble
python pilot/scripts/proxy_select.py pilot --ensemble --noise --tag full
# held-out test set (pre-registered: docs/internal/preregistration_A1_zh.md)
TEST=FOXE1,GP1BA,HBB,HBG1,IRF4,LDLR,LDLR.2,MSMB,MYCrs11986220,MYCrs6983267,SORT1,SORT1-flip,SORT1.2,TERT-GAa,TERT-GBM,TERT-GSc,TERT-HEK
python pilot/scripts/alphagenome_score.py pilot --all-biosamples --outputs CAGE,DNASE,ATAC --versions ALL_FOLDS --elements $TEST
python pilot/scripts/proxy_select.py pilot --elements $TEST --tag test
python pilot/scripts/prereg_test.py pilot
python pilot/scripts/a1_report.py pilot results/a1
```

## Data and model terms

- **Ground truth:** Kircher M, *et al.* (2019) *Nat Commun* 10:3583, doi:10.1038/s41467-019-11526-w. The data are
  downloaded from the authors' portal repository and are not redistributed here; cite the publication for any use.
- **AlphaGenome:** API and weights are available for non-commercial use under Google DeepMind's terms of use.
- **Borzoi:** Linder J, *et al.* (2025) *Nat Genet*; PyTorch weights via `borzoi-pytorch` (check the licence of
  the port you use).

## Licence

To be decided before the first public release.
