# Bilingual manuscript figures

This directory contains eight original figures in English and Chinese. Each
figure is supplied as an editable SVG, a vector PDF and a PNG preview. The SVG
files are the editable masters. PDF text is converted to vector outlines to
preserve the Chinese typography across PDF readers; no raster image is embedded.

## Rebuild

From the manuscript root, run:

```bash
python code/redraw_figures.py
python code/write_figure_captions.py
python code/verify_redrawn_figures.py
```

The drawing script needs Python, NumPy, Inkscape, DejaVu Sans and Noto Sans CJK SC.
The verification script additionally needs PyMuPDF. No network access is needed.
Theoretical examples are deterministic. The empirical figure reads the archived
summary tables supplied in `source_data`; its counts are checked against those
tables while building the figure.

## Files and placement

Every figure has `_en` and `_zh` editions. Use the corresponding file under
`pdf` for vector insertion in a PDF typesetting workflow, or the corresponding
file under `svg` where SVG insertion is supported. PNG files are 2040 pixels
wide and are intended for previews or document workflows without SVG support.

| New figure | Stem | Scientific role | Relation to earlier figure set |
| --- | --- | --- | --- |
| 1 | `fig01_foundations` | SR intuition, Euclidean geometry equivalence, exact learned-reference alias | Earlier Fig. 1 plus new introductory panels |
| 2 | `fig02_families` | Complete continuous families, exact free-current compensation, multiscale equivalence | Earlier Fig. 2, reorganized |
| 3 | `fig03_reciprocity` | A common reciprocal task plane and the four competing families | Earlier Fig. 3, reorganized |
| 4 | `fig04_measurement` | Compatible readout and joint uncertainty | Earlier Fig. 6, reorganized |
| 5 | `fig05_quantitative` | Parameter sets, independent propagation and two-plane recovery | Earlier Fig. 4, reorganized |
| 6 | `fig06_kernels` | Finite memory, visible temporal weights and a learned reverse reference | Earlier Fig. 5, reorganized |
| 7 | `fig07_calibration` | Real-data calibration and known-transform evaluations | Earlier Fig. 7 plus archived dFB results |
| 8 | `fig08_programme` | The author's general-model/particular-model programme | Earlier Fig. 8, retaining its two-level argument |

Captions are supplied in `captions_en.md`, `captions_zh.md` and as individual
Markdown files in `captions`. The structured caption/source mapping is
`source_data/redrawn_figure_captions.json`. All additional analytic examples,
including the exact matrix entries behind the simplified diagrams, are recorded
in `source_data/redrawn_analytic_values.json`.

## Dimensions and typography

The drawing width is 170 mm. Heights vary from 110 to 146.25 mm. The smallest
label is 8.15 pt at this width, or 7.67 pt when the entire figure is placed at
160 mm. Panel letters are lowercase and aligned. The blue/orange direction
palette is consistent across all figures; green denotes additional constraints
or task-compatible measurements, and purple denotes mixtures or further work.
Detailed interpretation and assumptions are placed in the captions.

## Data and interpretation

Figures 1–6 and the analytic part of Figure 8 are deterministic constructions,
not estimates from new experiments. Figure 7 preserves the archived counts and
ratios. Its 864 conditions are analysis conditions, not independent experiments;
the 0/6/77 counts describe a deterministic sufficient-condition check using the
actual residual. Its three known-transform declarations and its physical-block
and cross-record evaluations retain their separate scopes in the caption.

The `I + J` compensation in Figures 2 and 8 belongs to the family with a freely
varying current-state coefficient. It is not an admissible positive-gain complete
SR with a discount below one. Figure 2c instead gives a legitimate complete-SR
mixture and its positively scaled effective complete SR on one task plane.

