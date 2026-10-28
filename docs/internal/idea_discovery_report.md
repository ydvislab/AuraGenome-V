# Idea Discovery Report

> Copied from the ARIS run directory on 2026-10-29. Cross-model novelty verdict for A1: PROCEED WITH CAUTION (7/10); cross-model research review: PROCEED WITH CAUTION, readiness 5/10 (see "External Critical Review"). Current project state: [`docs/STATUS_zh.md`](../STATUS_zh.md).

**Direction**: concrete, practical problems where CS (visualization, uncertainty quantification, LLMs) advances bioinformatics analysis of genomic foundation-model variant-effect predictions. Target: *Bioinformatics* (Original Paper).
**Date**: 2026-10-28
**Pipeline**: research-lit → idea-creator → novelty-check → research-review → research-refine-pipeline
**Run id**: genomic-fm-vep-verification-20261028

## Executive Summary
(filled at Phase 5)

## Literature Landscape

Sources contributed: arxiv (ARIS fetcher), web (WebSearch). Gemini unavailable (CLI not installed); Zotero/Obsidian MCP not configured.
WARN: local contributed nothing — no PDFs found in papers/, literature/, or a configured paper library.
Verification (verify_papers.py): 10 verified, 8 verify_pending (bioRxiv/medRxiv preprints without DOI in the query; their pages were read earlier in this session), 0 hallucinated.

| Paper | Venue | Method | Key result | Relevance | Status |
|---|---|---|---|---|---|
| Avsec et al., AlphaGenome | Nature 2026 | 1 Mb sequence → 5,930 human tracks | SOTA on 25/26 VEP benchmarks; **CAGI5 saturation MPRA (= Kircher 2019) is one of its benchmarks**; cell-type-matched DNase r=0.57; LASSO over tracks makes matching unnecessary when labels exist | Our ground truth is already a model benchmark → "how accurate is the model" is not our contribution | ✅ crossref |
| Corgi (context-aware S2F) | Nat Commun 2026 | sequence + trans-regulator expression → predictions in **held-out cell types** | learns cell-type regulators without training on them | Closest work to unseen-cell-type VEP; requires a new model, not post-hoc use of AlphaGenome | ✅ crossref |
| AlphaGenome commentary | Cell Res 2026 | — | S2F models "implicitly assume a fixed cellular state"; may only extend to related cell types | States the unseen-cell-type gap | ✅ crossref |
| DNACipher | bioRxiv 2025 | contrasts predicted signals to impute unseen cell-type/assay combinations | most contexts are imputable | Competing angle on unseen contexts | … pending |
| He, Palamuttam & Danko | bioRxiv 2024 / 2025 | training on personal genomes | VEP retained when fine-tuning to unseen contexts | fine-tuning route | … pending |
| Kircher et al. 2019 | Nat Commun | saturation MPRA, 20 elements | per-base effects; several elements measured in cell lines models do not profile | our ground truth | ✅ crossref |
| Shigaki et al. 2019 (CAGI5 regulation) | Hum Mutat | assessment of challenge submissions | integrating several epigenomic marks helps; discusses *choosing among cell types* for training | prior evidence that cell-type choice matters | ✅ crossref |
| Bajwa et al. 2024 | MLCB (PMLR 240) | 5 Basenji2 replicates | replicates disagree in >50% of eQTL cases; overconfident on reference | uncertainty background | … pending |
| RegTrust-XAI (Trust-Aware S2F) | bioRxiv 2026 | ensemble consensus × attribution coherence × applicability domain | single-sequence trust badge (HF demo); PKLR MPRA ρ≈0.23 with a small CNN | nearest "trust" competitor | … pending |
| Dark Regulome | arXiv 2606.06834 | residualization/permutation diagnostics | DNA-LM ISM tracks sequence predictability; models' top lists disjoint | cautions on scores | ✅ arxiv |
| Position: beyond anecdotal interpretability | arXiv 2606.07607 | benchmark of attribution methods | methods contradict each other, miss known motifs | explanation validity gap | ✅ arxiv |
| EVEE | bioRxiv 2026 | Evo 2 probes + LLM narrative explanations, 4.2M ClinVar variants | explanations judged only by LLM-as-judge (3.8/5); no faithfulness check, no user study | LLM-explanation gap | … pending |
| AlphaGenome Atlas | medRxiv 2026 | precomputed AVI scores, attributions, motif maps, GUI | no uncertainty or disagreement shown; "only part of the evidence chain" | portal competitor | … pending |
| tangermeme | Nat Methods 2026 | attribution/ablation toolkit | fast ISM/DeepLIFT | tooling | … pending |
| Archaic introgression with S2F | bioRxiv 2026 | AlphaGenome on introgressed variants | accessibility predictions match experiments; expression no better than chance | modality-dependent reliability | … pending |
| CellPrism | arXiv 2608.01669 (VIS 2026 Best Paper) | VA for virtual-cell perturbation predictions | future work: show prediction uncertainty | VIS competitor (different model class) | ✅ arxiv |
| CAGI overview | arXiv 2205.05897 | community assessment | — | evaluation tradition | ✅ arxiv |
| K562 CRISPRi enhancer benchmark (blog) | genomicsxai 2026 | zero-shot S2F on CRISPRi | AlphaGenome r=0.67/0.45; all models under-predict distal effects | distance-dependent weakness | web |

**Landscape synthesis.**
(1) Model accuracy on standard benchmarks is high and already published, including on our planned ground truth (CAGI5 = Kircher MPRA); our own pilot agrees (PKLR in K562: 96% direction accuracy). A paper whose premise is "models are often wrong" would be weak.
(2) The recognised weak spots are *contextual*: cell types absent from training (commentary, Corgi, DNACipher), expression vs accessibility readouts (introgression study), distal vs proximal effects (CRISPRi benchmark). The practitioner-facing version of (2) — *I have a cell type the model does not cover; which of its 5,930 tracks should I read, and how much does that choice change my answer?* — is handled either by retraining (Corgi, fine-tuning) or by supervised LASSO with labels (AlphaGenome CAGI5 protocol), not by a label-free, post-hoc procedure.
(3) Explanations: attribution methods disagree and LLM narratives (EVEE) are unchecked; saturation MPRA gives per-base ground truth that could verify mechanistic claims, but nobody has used it to audit LLM explanations.
(4) Trust signals exist only as research diagnostics (Bajwa, RegTrust) or a single badge; no published tool ties them to a concrete decision.

## Ranked Ideas

Generation: 12 ideas (gpt-5.6-sol) + 8 ideas (gpt-5.6-terra) → 15 after mechanical dedup, all budget-feasible (API-only). gpt-6-astra (configured reviewer) was at capacity; gpt-5.5 / gpt-6 / gpt-6.1-sol are not available on this ChatGPT account; both generator and judge are therefore GPT-5-family (cross-family with the Claude executor). Traces: `.aris/traces/idea-creator/2026-10-28_run01/`.

Cross-model triage (gpt-5.6-sol, thread 01a18187…): **TOP: A1, C1, E1.**

| Rank | Id | Idea | Reviewer's best case | Named risk |
|---|---|---|---|---|
| 1 | A1 | **Eight-variant proxy picker** — measure only the few variants on which candidate tracks disagree most, use them to pick the proxy track for an unprofiled cell type, predict the rest | a citable protocol: adapt AlphaGenome to a new cell type with ~8 measurements | retrospective masking may overstate prospective value; best proxy may vary by variant ("no stable best proxy") |
| 2 | C1 | **Locus-exposure delta (FoldSwitch)** — compare ALL_FOLDS with the checkpoint that held the locus out | first per-locus rule for checkpoint choice; either outcome is informative | ALL_FOLDS vs FOLD differ in distillation too → diagnostic, not causal; signal may vanish |
| 3 | E1 | **Explanation compiler** — explanations must emit executable claims (factor, span, bases, direction, evidence track) tested against the saturation-MPRA footprint | makes genomic explanations falsifiable; best fit for team's LLM + vis identity | MPRA cannot prove TF occupancy; many claims may be undecidable |
| 4 | D1 | Context-delta test (TERT cell lines, PKLR time points) | tests "matters more in context A than B" | proxy confounding; noisy interaction truth |
| 5 | A3 | Track archetypes / proxy medoids | compresses 5,930 tracks into a few distinct conclusions | generic clustering; locus-specific clusters |
| 6 | B1 | Proxy-unanimity gate | simple robust vs context-sensitive call | arbitrary consensus set; consensus without discrimination |
| 7 | A4 | Expression-nearest proxy retrieval | label-free proxy from RNA-seq | batch effects; close to Corgi conceptually |
| 8 | E2 | Three-allele motif triangulation | stringent use of all three alleles | narrow; allele resolution |
| 9 | D2 | Contrast ladder | sanity benchmark | incommensurable contrasts |
| 10 | A5 | One-locus context adapter | cross-locus transfer | few cell-line × element combinations |
| 11 | A2 | Activity-gated aggregation | one-line default | too obvious |
| 12 | B2 | Rank envelopes | conservative shortlists | generic |
| 13 | B3 | Four-checkpoint voting | cheap stability badge | not a valid ensemble; Bajwa 2024 precedent |
| 14 | F1 | Mutational-footprint segmenter | compact intervals | wrong representation assumption |
| 15 | F2 | Local-rank calibration | common unit | not calibration |

## Novelty Verification

**Status (2026-10-29): A1 verdict obtained — PROCEED WITH CAUTION 7/10 (see below). C1 and E1 were not re-checked (deferred by the user).** Earlier status: BLOCKED — cross-model verdict not yet obtained. The Codex reviewer returned "usage limit reached, try again at 8:49 PM" (2026-10-28). Phase B (multi-source search) is complete; Phase C (cross-model verdict) must be re-run. Dossier: `idea-stage/NOVELTY_DOSSIER.md`. The notes below are search findings, not verdicts.

| Idea | Closest prior work found | What stays different (to be verified by reviewer) |
|---|---|---|
| A1 disagreement-selected calibration panel | AlphaGenome CAGI5 protocol (LASSO over tracks with *all* labels); GenoNet / semi-supervised MPRA models (Bioinformatics 2021, limited labels, train a new model); active learning for tumour variants (bioRxiv 2021, small logistic model); neuronal MPRA 2026 (averages Enformer DNase tracks); Genome Biology 2024 (models degrade in cell-type-specific regions) | selecting *which existing foundation-model track to trust* for an unprofiled cell type from a handful of disagreement-chosen measurements — no paper found that does this |
| C1 locus-exposure delta | AlphaGenome docs/paper: FOLD_k used for aggregate out-of-fold evaluation; Flashzoi notes train–test similarity concerns; no per-locus ALL_FOLDS-vs-held-out comparison found | per-locus, user-facing diagnostic of training exposure; confound: distillation |
| E1 explanation compiler | Multi-ISM (bioRxiv 2026-09) benchmarks *attribution maps* against CAGI5 saturation MPRA; GeneAgent (Nat Methods 2025) verifies LLM gene-set claims against *databases*; EVEE (LLM-as-judge only); position paper 2026 on anecdotal interpretability | executable, structured explanation claims (LLM or attribution) falsified against *per-base experimental footprints* — no paper found; Multi-ISM is the nearest race and must be cited |

### Update 2026-10-28 20:55 — reframed A1 only (verdict still BLOCKED)

Codex returned "usage limit … try again at Oct 29th, 2026 1:52 AM" for both gpt-6-astra and gpt-5.6-sol. The
dossier for the reframed A1 is `idea-stage/NOVELTY_DOSSIER_A1.md` (method frozen, pre-registered held-out test GO,
few-label LASSO/ridge baselines). Additional Phase-B finding, executor search only, not a verdict:

| Prior work | Overlap | What stays different (to be verified by the reviewer) |
|---|---|---|
| MODEL SELECTOR (Okanovic et al., arXiv 2410.13609; AISTATS 2025) | label few examples to identify the best of >1,500 pretrained classifiers | regression with an unknown per-candidate affine map, marginalised; posterior-weighted prediction rather than a single pick; thousands of correlated output heads of one genomic model; unprofiled-cell-type application with a pre-registered held-out test |
| Model Picker (Karimi et al. 2021), CAMS (NeurIPS 2024) | streaming active model selection among pretrained classifiers | same as above; pool-based, ≤16 labels |
| Box & Hill 1967; Bayesian design for model discrimination | sequential discrimination design | scale (≈1,000 correlated candidates) and the genomic setting |
| AlphaGenome (CAGI5 LASSO over tracks) | uses labels to combine tracks | with 8 labels our method beats it by +0.117 Spearman on the held-out regions (post hoc) |

Implication for the paper: MODEL SELECTOR must be cited and, ideally, adapted to regression as a baseline.

### Cross-model verdict 2026-10-29 01:56 — A1: PROCEED WITH CAUTION, 7/10

Reviewer: Codex MCP, requested model gpt-6-astra (xhigh, read-only); the reply did not self-report a model id.
Thread 01a183a7-289f-7440-96fe-a0ce59f0e22c; trace `.aris/traces/novelty-check/2026-10-29_run01/`.

- **Clears the check**: no published paper contains the combination of small-panel calibration, predictive
  averaging across genomic tracks, and the validation.
- **Delta (verifiable)**: "Unlike pretrained-classifier selection and AlphaGenome's supervised CAGI5 aggregation, A1
  uses a tiny, adaptively acquired panel of quantitative variant effects to marginalize affine calibration across
  frozen genomic tracks and produce locus-specific predictions for unmeasured variants."
- **Why caution**: the acquisition score is posterior-weighted pairwise disagreement (Box–Hill mean separation), and
  its edge over simple targeted panels is unproven (H3 failed; disc − static disagreement only +0.010, 7/11).
- **New prior work to cite**: CODA (Kay et al., ICCV 2025), SELECT-LLM (arXiv 2510.09418v2), Gardner et al. NIPS
  2015, Olofsson et al. ICML 2018, Sugiyama & Rubens, Sawade et al. 2008/2012, Enformer Celltyping, Borzoi PEFT.
- **Correction**: AlphaGenome's CAGI5 LASSO uses the challenge training labels, not all labels.
- **Audit findings to fix**: (1) `hide_matched` also removes legitimate related proxies — add the deployment
  comparison with them retained (IRF4: metadata proxy 0.652 vs disc 0.575); (2) undefined LASSO Spearman values
  silently dropped (221/850 random+LASSO at test k = 8); (3) "oracle" is a reference, not an upper bound;
  (4) score all strategies on a common outside-pool set.
- **Second dataset**: Siraj SatMut is public (Zenodo 15297965, `data_satmut.zip`).
- **Defensible claim now**: at eight measurements, the Bayesian adapter beats random-panel adaptation and
  indiscriminate track averaging; the value of sequential acquisition over simpler targeted panels is unresolved.

## External Critical Review

**2026-10-29, Codex MCP (requested gpt-6-astra, xhigh; the reply did not self-report a model id), thread
01a184bc-aa73-7960-9842-81142cc6d582, trace `.aris/traces/research-review/2026-10-29_run01/`.**
Bundle: `idea-stage/REVIEW_BUNDLE_A1.md`, which includes the post-test deployment finding (v1 − matched-track mean
−0.032, 3/11 regions).

**Bottom line: PROCEED WITH CAUTION, not abandon.** Submission readiness is 5/10 now and is expected to reach 7/10
only if the plan below succeeds on frozen external data.

**Defensible claim today:** with metadata-selected tracks withheld, 8 measured variants improve prediction over
random-panel adaptation and over indiscriminate track averaging. Advantages over simple targeted panels and over
available metadata proxies are not established. Losses against the matched mean exceed 0.05 in 7/11 regions.

**What the "safety net" claim needs:** one frozen procedure that shows rescue and preservation together on
independent data, with:
- per-locus paired differences on identical evaluation variants;
- a noninferiority margin and a harm-rate bound, both pre-specified;
- "poor proxy" defined without the outcome.

**Required before more tuning:** a common outside-pool evaluation set for all methods and references. Corrected v1
is reported separately and not retuned.

**v2 design constraints:**
- fixed or very small grid, with leave-one-element-out sensitivity;
- an unchanged-proxy control and a metadata-centred shrinkage baseline (the affine fit can flip the matched mean);
- a family prior defined without labels and invariant to duplicated tracks;
- Siraj stays untouched until v2 is frozen.

**Scope of the Siraj data:** it is external-locus validation, not unprofiled-cell validation. Handle sequence
backgrounds and locus grouping.

**Mandatory baselines:**
- unchanged matched mean and all-track mean;
- matched-mean calibration and metadata shrinkage at the same budget;
- the acquisition strategies under identical inference;
- one-round batch panel;
- budget-matched ridge/LASSO;
- ablations;
- continuous EIG if active design stays central.

Citing MODEL SELECTOR and CODA is mandatory; implementing them is conditional.

**Order of work:** evaluation fix → minimal metadata-aware adapter and controls → freeze and pre-register on Siraj
→ Siraj and C1 → decide on manuscript → interface.

**Stop rule:** drop the safety-net claim if v2 cannot beat simple metadata-centred adaptation, or cannot deliver
rescue with controlled harm.

## User decision (2026-10-28, before the cross-model verdict)

- **Selected direction: A1 only**, reframed so the core contribution is a CS method: *label-efficient discrimination among thousands of correlated foundation-model output tracks for a cell type the model does not profile* (choose the few variants whose measurements best separate candidate tracks; stopping rule; visual-analytics interface), validated on real genomic data (Kircher 2019 MPRA). Target: *Bioinformatics* Original Paper.
- C1 → robustness analysis inside the A1 paper (does the method hold on the held-out checkpoint?). E1 → deferred (later paper).
- Supporting evidence (pilot, PKLR-48h, ALL_FOLDS): replacing matched K562 tracks with HepG2 tracks lowers direction accuracy on significant SNVs from 0.95–0.96 to 0.69–0.85 (ATAC 0.69, RNA-seq 0.77, CAGE/DNase 0.85) and Spearman from ~0.90 to 0.30–0.70; K562 and HepG2 tracks disagree in direction on 13–32% of significant SNVs.
- Open for the reviewer: is the reframed algorithmic problem a recognised CS contribution beyond uncertainty sampling, and what is the closest prior method (optimal experimental design for model discrimination, active model selection)?
