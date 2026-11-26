# Continuous-family calibration and historical analysis sources

The `dFB_` files are verbatim compact tables or numerical summaries from `SR_dFB_complete_20261124.zip`, except `dFB_error_space_counts.csv`, which is reaggregated from its condition-level data. The original archive contains 155,150,177 bytes and has SHA-256 `8e925706bb2935709b08540e6314771ff4247b4b9cca358937e3d21471d882b8`. Exact original paths and per-file hashes are in `empirical/dFB_provenance.json`.

`code/empirical_check_dfb_tables.py` joins 912 condition-level rows, selects the 864 imposed forward/backward conditions, verifies the saved sufficient inequalities and recreates the full/common-space/current-removed counts. These are deterministic calibration conditions based on a known transform and a computable repeat discrepancy. They carry no nominal probability-coverage claim.

The `historical_` tables are verbatim outputs of the second-round analyses embedded in that archive. Their original paths and hashes appear in `empirical/historical_source_manifest.json`. They support Supplementary E10's analyses of learned reference content, actual event-history mixing, spike-history prediction and measured temporal-kernel transfer. The 16-branch catalog in `research_branch_catalog.csv` corresponds to E11. The first-round formal report and claim-to-evidence table are preserved in `reference_archive/`, with metadata in `empirical/round1_source_manifest.json`.

Current Turner event-array reference and paired-residual calibration are documented separately in `README_fig07.md`. The portable script `code/empirical_check_fig07_portable.py` recreates those reference estimates and frozen readouts. Its 17 regrouped pairs and nominal residual-test threshold are distinct from the 864 dFB deterministic conditions.

The exact finite-cycle example in Supplementary E10.5 is supported by `reference_archive/cyclic_lag_scope_proof.md`; `code/empirical_check_cyclic_lag_coordinates.py` reproduces the five-cycle numerical instance from the archived check. The example concerns finite coordinates and actual lag weights, and introduces no simulated neural participants.
