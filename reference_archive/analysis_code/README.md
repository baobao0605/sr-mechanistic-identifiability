# Historical analysis algorithms / 历史分析算法

This directory preserves the final analysis code that generated the archived dFB and second-round research results. It is a source-code excerpt from `SR_dFB_complete_20261124.zip`, retaining the original `sr_dfb_round3_20261124/` directory tree. Every copied file is byte-for-byte identical to its ZIP member; `manifest.json` gives the original relative path, byte count and SHA-256. The excerpt contains 58 original files, including 36 Python sources, with 306,623 bytes of original content.

本目录保存产生 dFB 及第二轮研究结果的最终算法源文件，保留原始目录层级。`manifest.json` 逐项记录原归档路径、字节数及 SHA-256。复制文件与原 ZIP 中对应条目逐字节一致；Python 文件仅作语法解析检查，本次打包没有执行这些历史分析入口。

## Current self-contained checks / 当前可直接核验的入口

The submission's portable empirical checks are located in the main package's `code/` directory. Run them from the package root:

```bash
python code/empirical_check_fig07_portable.py
python code/empirical_check_dfb_tables.py
python code/empirical_check_cyclic_lag_coordinates.py
```

These checks use the submission's included measurement excerpt and source tables, or the stated algebraic example. Their executed results are recorded in `empirical/fig07_portable_checks.json`, `empirical/dFB_table_checks.json` and `empirical/finite_kernel_moment_scope.json`. They check reference estimates/readouts, deterministic sufficient-condition counts, and the cyclic lag-coordinate example, respectively.

本提交包中可直接运行、且已经核对结果的经验入口是上述 `code/` 脚本。历史脚本保留原数据依赖和运行约定，没有为新根目录重写路径。历史 README 中的“包内缓存”“运行”等说明指完整原始归档，不能理解为本轻量代码摘录已包含全部输入。

## Included algorithms / 所含算法

All paths in the following table are relative to `sr_dfb_round3_20261124/`.

| Directory | Included computational role | Input dependencies retained by the original code |
|---|---|---|
| `theory/` | Continuous point-to-family and family-to-family separation, Krylov/spectral error projection, direct implementation checks, mathematical notes | NumPy/SciPy and the documented operator/reference inputs; the algebraic checker supplies its own specified examples |
| `independent/` | Kern/Turner input reconstruction; frozen-rule margins; target-refitted nearest-family decisions; current-equivalence classes; continuous-witness checks; AUROC precision; independent THINGS transfer audit | Original measured caches under `prior/restored/complete/`, the frozen second-round fit tables, derived `repeated_inputs.npz`, and intermediate margin/decision tables |
| `empirical/` | THINGS frozen-prediction reconstruction, continuous margins, corrected error indices, and disjoint six-session calibration/evaluation | `prior/restored/complete/analysis/D05/analysis_input.npz`, second-round frozen fit results, and the sibling theory module |
| `uniform/` | Spectral compression of actual THINGS references and continuous uniform separation calculations | The same THINGS measurement input, then the branch's spectral-case outputs |
| `prior/sr_round2_20261124/math/` | Original Kern/Turner calibration; shared SR templates, affine fits and nuisance residualization; corrected reference-ball calculation | Original Kern/Turner derived arrays and trial tables under `prior/restored/complete/`; source-specific metadata and active protocol files |
| `prior/sr_round2_20261124/calibration/` | THINGS calibration, cross-repeat contrasts and reference-compatibility calculations | Original THINGS arrays and the sibling `math/sr_geometry.py` module |
| `prior/sr_round2_20261124/history/` | Allen current/history count prediction and fine time bins; FASTIMAGES response extraction, kernel competition, nuisance geometry and temporal-rank prediction | Original event/count caches or Allen NWB input; FASTIMAGES derived MEG inputs or the original FIF for extraction; branch-specific sufficient statistics |
| `prior/sr_round2_20261124/ekman/` | Localizer content and direction contrasts; actual-history mixing; held-out history prediction; independent numerical checks | Author probability arrays, participant/sequence mapping, source notebook fields and the two continuous-fMRI participants' derived inputs |

`independent/build_inputs.py` reconstructs measured analysis inputs rather than building a document. `independent/make_report_tables.py` computes paired condition counts, discount summaries and projection ratios from numerical outputs; it contains the aggregation logic behind the reported tables. Both are retained for algorithm provenance. Active amendments define the final observation classes, error budgets, arithmetic ties and discount endpoints. The shared `math/sr_geometry.py` module is included because both the independent and calibration branches import it.

## Restoring historical inputs / 恢复历史输入

The complete source archive is identified by `../../empirical/dFB_provenance.json`:

- Archive: `SR_dFB_complete_20261124.zip`
- Size: 155,150,177 bytes
- SHA-256: `8e925706bb2935709b08540e6314771ff4247b4b9cca358937e3d21471d882b8`

The provenance record supplies the archive name and cryptographic identity; obtain the archive from an authorized location. `../../empirical/historical_source_manifest.json` and `../../empirical/round1_source_manifest.json` connect the numerical tables and earlier formal reports to their sources. Original input paths are stated in the scripts and retained branch READMEs. The source reports in the parent `reference_archive/` directory document their acquisition and processing units.

To reproduce a historical branch, restore the complete archive to a separate working directory, follow its README and active amendments, and provide any raw data required for optional re-extraction. The complete archive contains the derived caches and earlier result tables that this lightweight excerpt omits; some raw NWB/FIF or source-author materials may require the original documented dataset download. Reproduction should preserve the archive snapshot and write regenerated results in a working copy.

The original `requirements_analysis.txt` records the dFB release environment. The history README additionally specifies `threadpoolctl` and `h5py`; MNE is needed for FIF extraction. Historical algorithm execution under the present submission root has not been claimed or tested.

## Selection boundaries / 选取范围

Only final top-level implementations, their small active protocols/amendments, and explanatory README/algorithm notes are included. Superseded `v0`/`v1` backup directories, initial-grid diagnostic copies, raw observations, binary caches, result-table duplicates, plotting scripts and document/archive construction scripts are omitted. The retained original READMEs therefore describe files and commands beyond this excerpt; they remain unchanged as provenance. For the current manuscript's verified portable checks, use the main `code/` entries above.

The historical computations are known-transform calibration or dataset-specific reanalyses. Transform conditions, folds and nested state/channel settings retain their original dependence structure; the source code does not convert them into independent participants or new natural-task evidence.
