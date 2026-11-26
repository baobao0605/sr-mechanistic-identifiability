# Kern/Turner independent SR separation validation

Analysis date: 2026-11-24. Primary narrative: `report.md` (Chinese).

## Scope

This branch reconstructs 4,560 frozen scores from two real repeated-measurement subsets (Kern DSMR101 MEG, Turner ID13 EEG; one participant each), calculates continuous SR margins, compares deterministic error budgets with the matching decision rule, and audits the separate THINGS record-disjoint transfer. Computational conditions are reused and are not independent participants.

The nearest-family ablation refits the target. It is a different decision problem from frozen-training held-out prediction. Under the same observation-class rule, its 77 multi-class sufficient-condition cases had 39 visible direction recoveries, 37 current-only fits, and 1 symmetric fit under the old rule; the matching target-refit rule uniquely recovered all 77. The original candidate-label bookkeeping had 37 ties that all collapsed to current. The 38 conditions without direction recovery must not be described as 38 forward/backward direction errors. This is not a held-out prediction improvement claim.

## Reproduction order

Run from the parent project directory with NumPy, pandas, SciPy, scikit-learn and matplotlib installed. Set `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` for deterministic, efficient small-matrix calculations.

1. `python independent/build_inputs.py` — reconstruct inputs and all frozen scores. Requires the restored original source caches under `prior/restored/complete` and old code/results under `prior/sr_round2_20261124/math`.
2. `python independent/analyze_margins.py` — 48 folds, 912 conditions, 2,304 continuous geometry records. Uses `repeated_inputs.npz`, old frozen fit CSVs, old `sr_geometry.py`, and `theory/continuous_geometry.py`.
3. `python independent/summarize_margins.py` — four-fold aggregate results, descriptive AUROCs, exploratory source-only threshold learning and cross-source application. Scores are rounded to 12 decimals for arithmetic ties.
4. `python independent/profiled_decision.py` — target-level nearest-family decision-rule ablation; keep separate from held-out frozen predictions.
5. `python independent/verify_continuous_witnesses.py` — direct real-matrix checks of every saved continuous margin witness.
6. `python independent/audit_observed_classes.py` — apply the same current-equivalence rule to old and new decision endpoints, preserving old candidate labels.
7. `python independent/audit_auc_precision.py` — conservative score-ranking uncertainty from continuous margin brackets; not confidence intervals.
8. `python independent/make_report_tables.py` — paired threshold cases, discount-resolution and error-projection tables.
9. `python independent/make_figures.py` — exact figures in PNG, SVG and PDF.
10. `python independent/audit_things_transfer.py` — requires completed `empirical/session_disjoint/` and original THINGS `analysis_input.npz`; independently verifies all calibration rules, splits, hashes and 324 direct-matrix/SVD raw-input conditions.

The saved `repeated_inputs.npz` allows the main branch to run without reprocessing original MEG/EEG files. Old code and frozen fit CSVs must still remain at their documented relative paths. No original file is overwritten.

## Outputs

| File | Meaning |
|---|---|
| `report.md` | Complete Chinese scientific report |
| `margin_folds.csv`, `margin_aggregate.csv` | Fixed-model margins, actual-repeat error decompositions and original outcomes |
| `continuous_geometry.csv` | Continuous lower/upper brackets, witness parameters and numerical status |
| `profiled_decision_folds.csv` | Target-level nearest-family results and matching oracle inequalities |
| `profiled_multi_threshold_cases.csv` | Full paired records for the 77 multi-class gate cases |
| `profiled_binary_threshold_cases.csv` | Full paired records for the 80 binary gate cases |
| `paired_gate_summary.csv` | Compact paired counts, split by data source |
| `observed_class_audit.csv` | Raw candidate-label and unified observation-class endpoints side by side |
| `auc_precision_bounds.csv` | Conservative numerical-ranking intervals from continuous geometry brackets |
| `discount_resolution_summary.csv` | Continuous far-discount margins and old gamma reportability |
| `predictive_auroc.csv` | Descriptive empirical score rankings; no population confidence claims |
| `cross_dataset_thresholds.csv` | Exploratory source-only threshold choice and transfer results |
| `review_theory_v1.md` | Independent final mathematical review |
| `review_things_session_disjoint_v2.md` | Independent final record-disjoint validation review |
| `review_things_numeric_audit.json` | Numerically reproduced rules/indices and scope |
| `figures/decision_rule_and_projection.*` | Projection, decision-rule and discount figure |
| `figures/figure_independent_predictive_limits.*` | Cross-source empirical prediction limits |

## Preserved revisions

`method_review_v0.md`, `analysis_contract_v1.md`, `method_freeze_v0.json` precede the new calculations. `amendment_A01_profiled_decision.md` records the authorized decision-rule ablation. `amendment_A02_numerical_ties.md` explains why score ranks merge arithmetic ties; raw exploratory summaries and code remain in `v1_raw_numeric_ties/`.

`review_things_session_disjoint_v0.md` records the original request to correct the naive noise-only primary index; v2 confirms the correction. These are transparent retrospective analyses, not preregistered new experiments.

`amendment_A03_observation_classes.md` records the final cross-review correction of a candidate-label versus observation-class endpoint mismatch. Geometry, scores and original label fields are unchanged; new visible-unique-recovery endpoints are explicitly audited. Four-fold aggregate outcomes and their empirical AUROCs are unchanged.

The continuous-domain bounds use analytic inequalities evaluated in float64. They are not outward-rounded interval arithmetic. Empirical repeat differences have no demonstrated confidence coverage. Known-transform oracle checks use known computational truth and do not identify natural cognitive mechanisms in the original participants.
