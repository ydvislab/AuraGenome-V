# AuraGenome-V

**Verifying genomic foundation-model predictions with multi-source evidence.**

Genomic foundation models such as AlphaGenome and Borzoi predict the effect of any variant on thousands of
functional genomic tracks, but they give no indication of *where* a prediction can be trusted. AuraGenome-V
(working name) builds a reliability profile for every prediction, places the prediction next to independent
experimental evidence, and presents both in a multi-scale visualization with an evidence-constrained LLM
assistant, so that researchers can separate reliable from unreliable predictions during variant interpretation.

> **Status: pilot.** The repository currently contains the feasibility pilot (research plan, weeks 3–6):
> can model-internal reliability signals tell correct from incorrect variant-effect predictions, measured
> against saturation-mutagenesis MPRA ground truth? The visualization and LLM components come after the pilot.

## Repository layout

```
pilot/
  scripts/
    download_data.py        ground truth: Kircher et al. 2019 saturation-mutagenesis MPRA (pinned, checksummed)
    build_variant_table.py  SNV table (GRCh38) with measured effects and ground-truth classes
    alphagenome_score.py    AlphaGenome scores from ALL_FOLDS and FOLD_0-3 (in-silico mutagenesis per element)
    borzoi_score.py         Borzoi scores from replicates rep0-rep3 (GPU)
    compute_signals.py      reliability signals, AUROC, element-level bootstrap, pre-registered go/no-go rule
    proxy_select.py         proxy-track choice from a few measured variants (sequential discrimination panel vs baselines)
  data/                     downloaded and derived data (not tracked)
docs/internal/              research plan, novelty check and pilot notes (Chinese; internal working documents)
```

## Quick start (pilot)

Requires Python ≥ 3.10.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python pilot/scripts/download_data.py pilot
python pilot/scripts/build_variant_table.py pilot
python pilot/scripts/compute_signals.py pilot --simulate        # pipeline test on synthetic predictions
```

Real predictions:

```bash
export ALPHAGENOME_API_KEY=...                                  # non-commercial terms of use apply
python pilot/scripts/alphagenome_score.py pilot --dry-run
python pilot/scripts/alphagenome_score.py pilot --elements PKLR-48h --versions ALL_FOLDS
python pilot/scripts/alphagenome_score.py pilot
python pilot/scripts/compute_signals.py pilot
# proxy-track selection: score every biosample's tracks, then compare panel strategies
python pilot/scripts/alphagenome_score.py pilot --all-biosamples --outputs CAGE,DNASE,ATAC --versions ALL_FOLDS --elements PKLR-48h,F9,HNF4A,IRF6,BCL11A
python pilot/scripts/proxy_select.py pilot
# optional, GPU: pip install -r requirements-borzoi.txt
python pilot/scripts/borzoi_score.py pilot --fasta hg38.fa --targets targets_human.txt --check
```

**Pre-registered decision rule.** The pilot is a GO if the element-bootstrap 95% CI of the combined signal's
AUROC (direction correct vs. wrong) has a lower bound ≥ 0.6 and the top reliability tertile is at least 15
percentage points more accurate than the bottom tertile.

## Data and model terms

- **Ground truth:** Kircher M, *et al.* (2019) *Nat Commun* 10:3583, doi:10.1038/s41467-019-11526-w. The data are
  downloaded from the authors' portal repository and are not redistributed here; cite the publication for any use.
- **AlphaGenome:** API and weights are available for non-commercial use under Google DeepMind's terms of use.
- **Borzoi:** Linder J, *et al.* (2025) *Nat Genet*; PyTorch weights via `borzoi-pytorch` (check the licence of
  the port you use).

## Licence

To be decided before the first public release.
