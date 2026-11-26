# 连续 SR 几何与模型比较误差预算

完整数学条件见 `proofs.md`，实现为 `continuous_geometry.py`。只支持已知的实 normal 压缩传播算子，所有模型使用同一固定测量度量，仅自动消去共同的自由当前系数。若加入其他 nuisance，需要重新构造相应投影与连续界。

## 点到连续竞争族

```python
from continuous_geometry import SpectralSRGeometry
geo = SpectralSRGeometry.from_reference(P, C)
bracket = geo.point_margin(
    "forward", 0.5, "backward",
    intervals=[(0.0, 0.95)], atol=1e-4, max_evals=2049,
)
```

`lower`、`upper` 是单位真实预测残差到整个连续竞争非负射线族的距离。竞争增益允许为零；当前系数自由。它们不是前后两个单位向量的弦距离。实际距离需要乘以该参考上真实预测残差的范数：

```python
r_true = alpha_true * geo.residual_norm("forward", gamma_true)
margin_lower = r_true * bracket["lower"]
```

本轮 oracle 校准必须使用原冻结 `T_true` 在测试参考 `C` 上的真实残差长度。原 `alpha_true` 由训练设计固定，不能在测试参考上重设原真值，也不能把训练端 `rho * ||A||` 直接当作测试端残差。

`intervals` 支持不相连闭区间。例如判断真折扣 0.5 与相差至少 0.2 的同类模型，可给 `[(0,0.3),(0.7,0.95)]`。所有端点需在 `[0,1)`。

## 两个连续折扣族与强度范围

```python
bracket = geo.family_margin(
    "forward", "backward", atol=1e-3, max_evals=20001,
)
```

- `lower/upper`：真族残差幅度至少为 1，竞争族幅度可为 0。
- `chord_lower/chord_upper`：两侧残差幅度都至少为 1。
- 乘以规定的共同幅度下界得到对应模型族距离。
- `cosine_lower/cosine_upper`：连续双折扣域最大相关的上下界。
- `gamma_witness`：实际评估到的可行最优参数；它本身只给分离距离上界。

这两个强度问题不同，不能用较大的弦距离认证允许零增益的竞争族。

## 共同预测空间与 oracle 条件

```python
mu = T_true @ C
y_test = T_true @ D
epsilon = float(geo.krylov_norm(y_test - mu))
r_true_fit = float(geo.krylov_norm(y_test - prediction_from_training))
condition = margin_lower > epsilon + r_true_fit
```

`project_krylov` 和 `krylov_norm` 支持 `(q,p)` 或 `(...,q,p)` 数组。共同空间包含当前及所有候选，严格保留任意候选的平方误差差值。原相对损失分母、并列阈值仍使用原值。该条件仅在已知变换校准时直接可算；真实未知 SR 的 `T_true` 不可凭空提供。

## 明确的边界状态

| 输入情形 | 返回含义 |
|---|---|
| 真实模型没有可见预测残差 | `no_true_residual`；距离 0，不作方向认证 |
| 竞争模型为当前类 | 单位真实残差距离 1 |
| 两个族均退化成当前类 | 距离 0，`current_equivalence_merged` |
| 一个族退化成当前、另一族有预测 | 距离 1；退化类被合并为当前，不能声称该退化类满足正残差下界 |
| 参数区间为空列表 | 按数学空集合距离返回 `+inf`，并显式标记 `empty_rival_domain` / `empty_family_domain`；不代表任何非空竞争机制被认证排除 |
| 评估预算耗尽 | 保留当前上下界并标 `max_evaluations`；不得把 witness 冒充精确最优 |

实际科学数据运行不使用空竞争域。本接口只接受 F/B/S/current；`balanced_fb` 已被上一轮解析映射合并为相同对称观测类。

## 数值保证的准确表述

上下界来自覆盖连续区间的解析不等式，随后在双精度中求值。实现并未使用向外舍入区间算术，所以它们是连续域解析界的数值实现，不是对所有机器舍入均形式认证的证书。`normality_error`、`spectral_error`、模板数值零阈值、迭代容差均可检查。

`from_reference` 分组重复特征值的数值容差为 `1e-10`；本轮预定循环/Fourier 谱分离远大于此容差。真实数据的独立实现另用实 Krylov 矩阵 SVD 核对投影。若用于几乎重根的一般算子，应另外分析谱误差并调整实现，不应默认本轮数值容差适用。

## 核验

```bash
OPENBLAS_NUM_THREADS=1 python verify_continuous_geometry.py
```

输出 `geometry_verification.json`。这些条目是代数和实现核验，不是模拟被试或生物学证据。包含直接矩阵计算、独立 SVD 投影、SSE 差恒等、解析旋转、连续界与密集评价比较、未知参考歧义以及自由当前增益构造。
