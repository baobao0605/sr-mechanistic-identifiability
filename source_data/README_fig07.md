# Figure 7 source data and reproducibility

All Figure 7 empirical values are read from the included CSV files. No manually entered replacement measurements or simulated participants are used. The compact NPZ includes measured, baseline-corrected EEG excerpts, not raw EEG.

## Default portable execution

From the manuscript package directory:

```bash
python code/empirical_check_fig07_portable.py
python code/redraw_figures.py
```

The first command requires Python, NumPy, pandas and an authorized local copy of the optional NPZ (see README_DATA.md). It refits 40 reference estimators in eight physical-recording/block/window contexts from the included measured excerpts, reconstructs two frozen J01 readouts, checks the empirical summaries, and writes `empirical/fig07_portable_checks.json`. It requires no external AB directory and no absolute workspace path.

The second command redraws all eight figures in both languages, reading the empirical CSVs for Figure 7 and invoking Inkscape to export vector PDFs and PNG previews. Source files are placed in `figures/svg`, `figures/pdf` and `figures/png`. The physical figure width is 170 mm. No plotting Figure objects or generated illustrations are used.

## Files and data levels

| File | Content and unit |
|---|---|
| `fig07_reference_blocks.csv` | Eight first-plane reference-prediction endpoints: two recordings × two physical blocks × two windows, all from one participant |
| `fig07_reference_aggregate.csv` | Ratio of mean block losses; not a mean of percentage gains or a participant effect estimate |
| `fig07_all_reference_estimators.csv` | Forty zero/mean/CV-scalar/CV-SVD/CV-Fourier endpoints used for the portable refit check |
| `fig07_temporal_splits.csv` | Earlier/later event counts and onset gaps for the actual block splits |
| `fig07_acquisition_hierarchy.csv` | Physical acquisition hierarchy, not computational repeat counts |
| `fig07_regrouping_hierarchy.csv` | State-complete occurrence maps and alternating pairs reused in conditional calibration |
| `fig07_conditional_sr.csv` | Three nonduplicate standard backward-SR conditional positives from one recording/window; nominal test quantities |
| `fig07_within_physical_block.csv` | Two block-specific restricted reanalyses, both abstaining |
| `fig07_transport_conditions.csv` | All 1,200 cross-record computed conditions underlying the stress endpoint |
| `fig07_transport_summary.csv` | Dataset-specific counts of those computed conditions |
| `fig07_decision_summary.csv` | Figure denominators with explicit unit labels |
| `fig07_empirical_arrays.npz` | Compact event-level measurements, state/time/block metadata, saved B estimates and frozen readout parameters for the portable checks |
| `empirical_ekman_reference_all.csv` | All retained Ekman reference endpoints, including zero or negative results; supplementary only |

The NPZ `R1__measurements` and `R2__measurements` arrays have dimensions event × window × channel. Windows are ordered `early`, `late`, and channels are P7, P5, P3, P1, Pz, P2, P4, P6, P8, PO7, PO3, POz, PO4, PO8, O1, Oz, O2. Arrays `labels`, `time` and `block` retain source state labels, event onset times and physical block labels. Prefix `reference__` contains saved reference estimates; prefixes ending in `B_trained_A__` contain frozen state contrasts, whitening, templates and linear readouts. No target-based readout refitting is performed by the checker.

The early window 100–150 ms, late window 200–250 ms and baseline −200–0 ms are documented in the previous statistical appendix. The original supplied NPZ stores named windows but does not contain sample-index bounds or a participant-ID field. The prior report identifies the selected participant as ID13. Consequently the package verifies the downstream derived-array analysis, not the missing upstream extraction provenance or a raw-EEG preprocessing rerun.

## Optional reconstruction from the larger supplied archive

If the full `SR_AB_complete_20261125.zip` is available, extract it and run:

```bash
python code/audit_empirical_and_export_fig07.py --ab /path/to/SR_AB_20261125
```

This command reads the AB release's source CSVs and derived NPZ arrays, re-executes its independently written reference audit, and regenerates the compact source files. It is a separate optional step; default figure drawing and core checking do not require this archive. Input-relative paths and SHA-256 values are recorded in `empirical/empirical_audit_results.json`. The two upstream EEG feature-cache hashes are:

- recording 1: `b3242dbbe34b8e4ce16e28d03c239116c69665c7b83a9a248455f1eb7fe78375`;
- recording 2: `7a2414901ece2a1490b0fa4a445b5cca4d9cc9c35b3c9ada26ce0d86b968d846`.

Original Turner paper: https://doi.org/10.1371/journal.pbio.3003189 ; data https://osf.io/sn4a7/ ; analysis https://github.com/bootstrapbill/neural-location-decoding ; archived version https://zenodo.org/records/15283636 . These original public resources do not themselves supply our custom AB derived cache. The prior AB archive or its measured-array excerpt is therefore the correct starting point for reproducing the present analysis; downloading the original raw data alone does not silently reproduce undocumented extraction steps.

## What reproduction establishes

Portable checking establishes numerical agreement for derived measured-array reference fits, frozen readouts and source-table summaries. It does not establish independent Gaussian paired sampling, prospective familywise coverage, participant-level generalization, or endogenous SR identification. In the known-transform calibration, the entire measured signal-plus-noise array is transformed. Computational conditions and regrouped pairs are not independent participants or new task acquisitions.
