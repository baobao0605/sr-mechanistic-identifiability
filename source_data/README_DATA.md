# Data availability and provenance

The `research_branch_catalog.csv` file, manuscript and Supplementary E2/E8–E11 document the primary sources and roles of all 16 empirical branches. The small tables are deterministic constructions, analysis summaries and condition-level numerical results, not original recordings. No raw recordings are included.

## Optional measured input — do not redistribute without confirmed permission

`fig07_empirical_arrays.npz` is the baseline-corrected event-level excerpt derived from Turner's EEG data. The original source repository's redistribution terms were not established in this packaging audit; the public directory therefore excludes the array. Obtain an authorized copy of the exact excerpt and either pass its path to `reproduce.py` or place it in `source_data/` locally. It has SHA-256:

`1c4cdfeb625a195530f8a08b685a03f7f654a8d033b9db37e3152308f87f60ff`

- Turner paper: https://doi.org/10.1371/journal.pbio.3003189
- Data: https://osf.io/sn4a7/
- Code: https://github.com/bootstrapbill/neural-location-decoding
- Code archive: https://zenodo.org/records/15283636

The exact upstream raw-feature extraction and sample-index provenance were not recovered in the original audit. Re-downloading the raw EEG is not equivalent to reconstructing this exact file. An authorized copy of the excerpt or the original AB archive is needed for exact downstream refits; see `README_fig07.md`. Review original license and consent conditions before releasing any measured arrays.

## Other sources

- Kern MEG: https://doi.org/10.5281/zenodo.8001755
- Ekman fMRI: https://doi.org/10.34973/bsy6-9h29
- Tarder-Stoll fMRI: https://openneuro.org/datasets/ds005125
- THINGS: https://github.com/ViCCo-Group/THINGS-data
- All remaining branches: `research_branch_catalog.csv`.

`README_dFB_and_history.md`, `../empirical/`, and `../reference_archive/analysis_code/README.md` specify historical archives, hashes and dependencies. These large archives are not included; they must be separately restored from permitted sources for historical reruns. The included portable checks should not be described as complete raw-neural-data reprocessing.
