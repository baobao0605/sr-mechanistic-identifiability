# 时间历史分支复现说明

最终说明见 `report_zh.md`。本分支只有真实记录和真实试次派生量，没有模拟被试。

## 仅用包内缓存复跑

在此目录运行：

```bash
python allen_event_history.py
python allen_timebins.py
python fastimages_kernel_competition.py
python fastimages_kernel_geometry.py
python fastimages_temporal_rank_nested.py
python make_history_figures.py
```

需要 Python、numpy、scipy、pandas、matplotlib、threadpoolctl、h5py。mne 仅用于从连续 FIF 重新提取输入。`allen_timebins.py` 在计数缓存已存在时不需要原始 NWB；它生成并执行 `allen_timebins_model.py`。

## 从原始数据重新提取

`extract_fastimages_input.py /absolute/path/to/sub-01_task-main_proc-clean_raw.fif` 可重新提取全部 MEG 输入。它核验实际图像事件与 FIF 注释的时间原点，并生成无 pickle 的数值缓存。

`allen_timebins.py` 在 `allen_visp_25ms_counts.npz` 缺失时读取本轮相对工作树中的原始 NWB（本目录的同级 raw/D14/ 路径），从 110 个 VISp 单位的 spike times 重新计算十个25ms窗。若移到其他工作树，须更新这一原始输入路径。提供的计数缓存和 `independent_count_epoch_audit.json` 足以复跑所有关键模型与验证计数一致性。

## 输入和结果

- `inputs/natural_scene_counts.npz`：上一轮的实际 Allen 事件、作者质量选择、失效时段、原 50–250 ms/前刺激窗计数；本轮新分析使用响应窗和真实刺激标签。
- `allen_visp_25ms_counts.npz`：本轮从原始 spike times 提取的实际细时间窗计数。
- `fastimages_kernel_input.npz`：独立定位试次和实际快速序列的传感器观测、时间索引、刺激标签。
- `fastimages_kernel_sufficient_stats.npz`：逐训练折/逐真实序列的 Gram、交叉积、传感器均值和响应能量；可用于独立复算拟合与测试误差。
- `*_outer_folds.csv`：所有测试折的预测误差。
- `*_inner_selection.csv`：训练内部正则化选择完整记录。
- `fastimages_kernel_gamma_profiles.csv`：候选折扣的训练比较；这些参数未被验证为认知折扣。
- `*_summary.csv`：报告表格的数据来源。
- `initial_grid_diagnostic/`：上界触及的初始分析，诊断用途，**不是最终结果**。

`analysis_contract.json` 是运行前分析约定。`regularization_amendment.json` 和 `temporal_rank_amendment.json` 记录随后发生的修订及其具体原因。本轮为回顾性计算研究，不宣称预注册。

四个 Allen 时间折共享同一动物；八个 FASTIMAGES 块共享同一人。折间重复不提供人群样本量。报告中的几何投影仅条件于估计参考和指定干扰空间，不能直接解释成统计功效或机制正确概率。
