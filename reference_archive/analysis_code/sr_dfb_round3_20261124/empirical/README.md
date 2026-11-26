# THINGS real-repeat SR separation calibration

This branch connects continuous parameter-family separation with the existing frozen recovery analysis, and adds a record-disjoint 6+6 empirical transfer experiment. Read `report_zh.md` for results and limits.

## Reproduce

From the round-3 directory, run:

```bash
python empirical/run_pipeline.py
```

Required sibling inputs are `prior/sr_round2_20261124`, `prior/restored/complete/analysis/D05/analysis_input.npz`, and `theory/continuous_geometry.py`. The scripts resolve paths relative to themselves; no machine-specific absolute paths are required. Python dependencies are NumPy, SciPy, pandas, scikit-learn and Matplotlib.

The pipeline executes:

1. `analyze_margin_inputs.py`: reconstruct original frozen predictions, audit all scores, compute actual train/test error budgets and dense-grid witnesses. Does not refit original models.
2. `compute_continuous_margins.py`: continuous-domain analytic bounds through the theory branch API, with float64 arithmetic and explicit tolerances.
3. `summarize_margins.py`: join margins to recovery; compute corrected pair+q indices and oracle conditions; enforce the gamma endpoint definition.
4. `session_disjoint.py`: refit candidates separately in sessions 1–6 and 7–12. Calibration rules are written and hashed before the latter half is computed.
5. `audit_auc_precision.py`: numerical ordering bounds and gamma-stratified descriptive results.
6. `plot_results.py`: figures in PNG/PDF/SVG.

Re-running writes reproducible numeric outputs in this directory. Preserve the supplied archive and protocol files when inspecting previous versions. The historical v1 naive-radius archive is not an input to the corrected analysis.

## Definitions and scope

- `I_train_pair_q_*` is the primary training empirical index. `I_train_*` without `pair_q` in legacy joined tables is the obsolete naive comparison; do not confuse them.
- All margins profile a free current component and a nonnegative predictive gain. Rival gain zero is included.
- Oracle indices use the known imposed truth and test observations. They are retrospective deterministic checks, not confidence certificates for unknown mechanisms.
- Unit-margin intervals are analytic bounds evaluated numerically, not outward-rounded machine-certified intervals.
- `gamma_correct_0p2` uses `error < 0.2 - 1e-12`; the legacy floating-point endpoint result is retained separately.
- The 10080 grid conditions and 1680 new split conditions are computations on one participant, with repeated/nested stimuli and channels. No independence, participant-level confidence or FDR claim is made.
- THINGS was not a learned cyclic-transition task. Known transforms are imposed on real measurements, including their noise. These are method-calibration results, not evidence for neural SR in this dataset.

The two transparent amendments are `amendment_A01_error_budget.md` and `amendment_A02_gamma_endpoint.md`.
