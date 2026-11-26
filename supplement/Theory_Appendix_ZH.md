---
title: "理论与统计附录：从 SR 预测族的交叠到性质与参数识别"
lang: zh-CN
date: "2026-11-25"
---

# 理论与统计附录

本附录按“观测对象—不可辨边界—实验解除—可观测读出—不确定性”组织证明。所有方向声明均相对于预先规定的状态标签与转移方向。各命题的编号 T01–T28 在正文、理论清单及复现文件中一致。这里区分三种不同目标：恢复任意完整算子；辨认受限模型族共同具有的性质；在该族内估计参数。三者需要的信息不相同。

## A1. 模型、状态对比与观测层级

### A1.1 符号

设 $P\in\mathbb R^{n\times n}$ 为行随机矩阵，$P_{ij}$ 表示从状态 $i$ 到状态 $j$ 的一步转移概率。假设有严格正的平稳分布 $\pi$，作为列向量满足 $\pi^\top P=\pi^\top$、$\mathbf1^\top\pi=1$。记 $D=\operatorname{diag}(\pi)$。平稳时间反转为

$$
P^*=D^{-1}P^\top D.
$$

这一表达来自 $\Pr(S_{t-1}=j\mid S_t=i)=\pi_jP_{ji}/\pi_i$。只有平稳分布均匀时，转置才直接表示时间反转。无须假设过程非周期。

包含当前状态的完整继任表征及纯未来部分分别为

$$
M_\gamma(P)=(I-\gamma P)^{-1}
=\sum_{h\ge0}\gamma^hP^h,
\qquad
N_\gamma(P)=P(I-\gamma P)^{-1},
\qquad 0\le\gamma<1.
$$

二者满足 $M_\gamma=I+\gamma N_\gamma$。因此 $M_0=I$ 无方向信息，但 $N_0=P$ 通常有一步方向信息。本文始终保留这一区别。

令 $\widetilde P=D^{1/2}PD^{-1/2}$、$r=\sqrt\pi$，并取列标准正交的 $U\in\mathbb R^{n\times(n-1)}$，使 $U^\top r=0$。定义

$$
A=U^\top\widetilde P U=S+K,
\qquad S=(A+A^\top)/2,
\qquad K=(A-A^\top)/2.
$$

$S$ 是对称部分，$K$ 是相对于任务方向的反对称部分。令 $F_\gamma(A)=(I-\gamma A)^{-1}$，$\operatorname{skew}T=(T-T^\top)/2$，$\langle V,W\rangle_F=\operatorname{tr}(V^\top W)$。在二维旋转块中统一采用

$$
J=\begin{pmatrix}0&-1\\1&0\end{pmatrix},
\qquad J^2=-I,
\qquad J^\top=-J.
$$

### 命题 T01：平稳对比与共同测量关系

有 $\|\widetilde P\|_2\le1$、$\|A\|_2\le1$、$-I\preceq S\preceq I$。若条件 $c$ 的总体状态模式满足

$$
Y_c=g_c M_c(P)B+\mathbf1b_c^\top,
\qquad g_c>0,
$$

其中 $B$ 为各条件共用的任意线性特征映射，$b_c$ 是状态无关偏置，则其加权状态对比满足

$$
X_c=U^\top D^{1/2}Y_c=g_c M_c(A)R,
\qquad R=U^\top D^{1/2}B.
$$

**证明。** Jensen 不等式给出

$$
\sum_i\pi_i\Big(\sum_jP_{ij}x_j\Big)^2
\le\sum_{ij}\pi_iP_{ij}x_j^2
=\sum_j\pi_jx_j^2.
$$

故 $\widetilde P$ 在欧氏范数下为压缩算子。又 $\widetilde P r=r$ 与 $\widetilde P^\top r=r$，所以 $r^\perp$ 是共同约化子空间；限制后的算子为 $A$，反转为 $A^\top$。有理函数限制到该子空间与先限制再取函数相同，状态无关偏置被 $U^\top D^{1/2}\mathbf1=0$ 消去。对称部分的谱界来自 $\|A\|_2\le1$。证毕。

**适用边界。** 固定记录通道或成功解码不能独自证明生物学映射 $B$ 在条件间相同。若允许任意 $B_c$，任意可逆候选 $M_c$ 均可令 $B_c=g_c^{-1}M_c^{-1}(Y_c-\mathbf1b_c^\top)$ 来解释观测。非线性分类概率也不自动满足上述线性关系。数学中的 $P$ 是候选认知机制使用的转移；它等于任务生成器仍须明示或独立约束。

## A2. 几何与参考造成的精确等价

### 定理 T02：规范几何的正规性边界

若 $P$ 双随机且 $0<\gamma<1$，则 $M_\gamma(P)$ 与 $M_\gamma(P^\top)$ 的完整平方欧氏行距离矩阵相同，当且仅当 $P$ 正规，即 $PP^\top=P^\top P$。令 $G_f=M_\gamma(P)M_\gamma(P)^\top$、$G_b=M_\gamma(P)^\top M_\gamma(P)$，还有

$$
G_f-G_b=\gamma^2G_f(PP^\top-P^\top P)G_b,
$$

$$
\frac{\gamma^2}{(1+\gamma)^4}\|PP^\top-P^\top P\|_F
\le\|G_f-G_b\|_F
\le\frac{\gamma^2}{(1-\gamma)^4}\|PP^\top-P^\top P\|_F.
$$

**证明。** 两个 Gram 矩阵的逆之差为 $\gamma^2(PP^\top-P^\top P)$。用 $X-Y=X(Y^{-1}-X^{-1})Y$ 得第一式。两者均将 $\mathbf1$ 映为 $(1-\gamma)^{-2}\mathbf1$，故其差已中心化。平方距离确定中心化 Gram 矩阵，因此完整距离相等当且仅当 $G_f-G_b=0$；可逆性和 $\gamma>0$ 随即给出正规性。由于 $\sigma_{min}(I-\gamma P)\ge1-\gamma$、$\|I-\gamma P\|_2\le1+\gamma$，对上述乘积应用奇异值上下界即得范数界。证毕。

**非均匀扩展。** 将行模式改为 $M_\gamma(P)D^{-1/2}$ 与 $M_\gamma(P^*)D^{-1/2}$，即明确采用特征度量 $D^{-1}$，结论变为 $\widetilde P$ 正规的充要条件。令

$$
\Delta=F_\gamma(\widetilde P)F_\gamma(\widetilde P)^\top-
F_\gamma(\widetilde P)^\top F_\gamma(\widetilde P),
\qquad G=D^{-1/2}\Delta D^{-1/2}.
$$

有 $G\pi=0$。若距离相同，对称矩阵 $G$ 必为 $\mathbf1v^\top+v\mathbf1^\top$，从而 $v+(v^\top\pi)\mathbf1=0$。再左乘 $\pi^\top$ 得 $v=0$。故距离相等等价于 $\Delta=0$，之后同前。未经加权的非均匀情况不由此证明。

**纯未来扩展。** 对 $N_\gamma(P)$，同一规范几何结论在 $\gamma=0$ 也成立。行列和相同保证距离相等等价于 $N_\gamma$ 正规。正规 $P$ 的有理函数正规；反向用

$$
P=N_\gamma(I+\gamma N_\gamma)^{-1}
$$

及谱定理证明，不要求 $P$ 可逆。相反，完整 $M_0=I$ 的等距不蕴含 $P$ 正规。

这些结论针对规范特征几何、相同折扣和完整距离，不能直接推广到任意各向异性读出 $M_\gamma B$、相关距离、部分状态对或两侧不同折扣的拟合。

### 推论 T03：高度有向但几何相同的图族

对至少三个状态的循环置换 $R$，令

$$
P=pR+qR^\top+(1-p-q)I,
\qquad p,q>0,\quad p+q\le1.
$$

它是正规双随机矩阵；当 $p\ne q$ 时不可逆，平稳熵产生为 $(p-q)\log(p/q)$，在固定 $p$、$q\downarrow0$ 时无界，但所有折扣下的前后规范几何相同。

更一般地，在 $\mathbf1^\perp$ 中取正交向量对 $(u_\ell,v_\ell)$，令

$$
K_\sigma=\sum_{\ell=1}^k\sigma_\ell\epsilon_\ell
(u_\ell v_\ell^\top-v_\ell u_\ell^\top),
\quad P_\sigma=\mathcal J_n+s(I-\mathcal J_n)+K_\sigma,
\quad\mathcal J_n=\mathbf1\mathbf1^\top/n,
$$

其中 $0\le s<1$、$\sigma_\ell\in\{-1,1\}$、$\epsilon_\ell>0$ 足够小以保证逐项非负。则 $2^k$ 个已标记转移矩阵在所有折扣下等距。

**证明。** 循环中所有项可交换，故正规；熵产生直接按相反有向边配对求和。第二个构造中各平面正交，故

$$
K_\sigma^2=-\sum_\ell\epsilon_\ell^2
(u_\ell u_\ell^\top+v_\ell v_\ell^\top)
$$

与符号无关；在状态对比空间中的 Gram 矩阵为

$$
\big[(1-\gamma s)^2I-\gamma^2K_\sigma^2\big]^{-1}.
$$

常数方向也相同，结论成立。证毕。该计数不是未标记非同构图计数，也不保证随着维数增加仍有统一的正分离距离。

### 命题 T04：直接参考的方向正性与已学习参考的完整别名

若 $X_0=g_0R$ 满行秩，$X_1=g_1F_\gamma(A)R$，则 $T=X_1X_0^+=(g_1/g_0)F_\gamma(A)$。独立转移参考 $P_0=\mathbf1q^\top$ 的已学习完整 SR 在状态对比后也具有 $g_0R$ 形式。令 $K\ne0$，定义

$$
\omega(T)=\frac{\langle T,K\rangle_F}{\|K\|_F^2}.
$$

对直接完整 SR，$\gamma>0$ 时前向为正、后向为负。对同一个图的两个已学习前向参考，

$$
T=aF_\gamma(A)(I-\eta A)
=a[I+(\gamma-\eta)AF_\gamma(A)],
$$

有

$$
\omega(T)=a(\gamma-\eta)\rho_\gamma/\|K\|_F^2,
\qquad\rho_\gamma>0.
$$

**证明。** 满行秩给出 $X_0X_0^+=I$。独立参考的非恒等部分各行相同，故在对比后消失。为证明正性，设

$$
H=I-\gamma S\succ0,\quad L=H^{-1/2}KH^{-1/2},\quad R_\gamma=(I-\gamma L)^{-1}.
$$

则 $F_\gamma(A)=H^{-1/2}R_\gamma H^{-1/2}$，逆矩阵之差给出

$$
\operatorname{skew}F_\gamma=\gamma F_\gamma K F_\gamma^\top,
\qquad
\operatorname{skew}(AF_\gamma)=F_\gamma K F_\gamma^\top.
$$

第二式在 $\gamma=0$ 仍成立。迹的循环性及 $R_\gamma^\top L R_\gamma=L(I-\gamma^2L^2)^{-1}$ 给出

$$
\rho_\gamma=\langle F_\gamma K F_\gamma^\top,K\rangle_F
=-\operatorname{tr}\{L^2(I-\gamma^2L^2)^{-1}\}>0.
$$

最后一步使用 $-L^2$ 非零半正定，不要求 $H$ 与 $L$ 可交换。取内积即得所有符号公式。证毕。

**完整对比算子的反例。** 七状态 Paley 转移取 $P_{ij}=1/3$ 当 $j-i\pmod7\in\{1,2,4\}$。其六维状态对比块满足

$$
A=-I/6+K,\qquad K^2=-7I/36.
$$

由 $F_t(A)=\{(1+t/6)I+tK\}/\{(1+t/6)^2+7t^2/36\}$ 可直接验算

$$
F_{3/10}(A)F_{4/5}(A)^{-1}
=\frac{53}{48}I-\frac{25}{56}K
=\frac{634}{519}F_{75/173}(A^\top).
$$

因此两个真正前向但时域不同的已学习模式，其相对算子严格等于一个合法后向完整 SR 乘正增益。等价发生在全部六维对比算子上，超出了欧氏几何；常数状态方向明确不在该等价中。确定性均值相等也不自动意味着噪声分布相等。

## A3. 连续族为何改变可识别的科学主张

### 命题 T05：观测相容集合上的性质与方向锚定

令 $\Phi:\Theta\to\mathcal Y$ 是包含所有合法机制和干扰参数的观测映射。性质 $q(\theta)$ 能由无噪声观测确定，当且仅当它在每个非空纤维 $\Phi^{-1}(y)$ 上为常量。

**证明。** 若两个参数映到同一观测却性质不同，任何只以该观测为输入的规则都不能同时正确；若各纤维上恒定，则定义 $\widetilde q(y)=q(\theta)$，其中任取 $\theta\in\Phi^{-1}(y)$，便得到良定义的观测性质。证毕。该命题是标准可辨识性原则；本研究的数学内容在于具体计算 SR 观测族与其纤维，而非重新命名这一原则。

若允许的转移集合 $\mathcal P$ 对平稳时间反转封闭，则

$$
\{M_\gamma(P):P\in\mathcal P\}
=\{M_\gamma(P^*):P\in\mathcal P\}.
$$

这是因为 $P^{**}=P$，两边只是重新编号。若需要方向主张，学得转移必须由设计或独立资料约束。若行为参数集合的漏盖概率至多 $\beta$，神经集合的漏盖概率至多 $\alpha$，则二者同时正确的概率至少 $1-\alpha-\beta$；这是事件并集界，不要求两类误差独立。

### 定理 T06：单平面直接广义 SR 的半平面塌缩

取单位圆上的正交块

$$
A_\sigma=aI+\sigma\beta J,\quad a^2+\beta^2=1,\quad\beta>0,\quad\sigma\in\{-1,1\},
$$

并规定合格直接参考。设

$$
G_\sigma=cI+\alpha N_\gamma(A_\sigma),
\quad c\in\mathbb R,\quad\alpha\ge0,\quad0\le\gamma<1.
$$

则

$$
G_\sigma=uI+vJ,\quad
D_\gamma=1-2a\gamma+\gamma^2>0,
$$

$$
u=c+\alpha\frac{a-\gamma}{D_\gamma},
\qquad v=\sigma\alpha\frac{\beta}{D_\gamma}.
$$

该连续族的像恰为 $\sigma v\ge0$ 的整个闭半平面；若严格要求 $\alpha>0$，像为开半平面 $\sigma v>0$。给定任意非零参考 $b$，观察 $y=G_\sigma b$ 可唯一识别 $u,v$ 和净旋转符号，但无法识别 $\gamma$。

**证明。** 用 $(I-\gamma A_\sigma)^{-1}=\{(1-a\gamma)I+\sigma\beta\gamma J\}/D_\gamma$ 相乘得系数。对任意合法 $(u,v)$ 及任意 $\gamma$，令

$$
\alpha=|v|D_\gamma/\beta,
\qquad c=u-|v|(a-\gamma)/\beta,
$$

即可得到同一算子，证明整个像及折扣别名。由于 $b,Jb$ 正交等长，

$$
u=\frac{b^\top y}{\|b\|^2},\qquad
v=\frac{(Jb)^\top y}{\|b\|^2}
=\frac{\det[b,y]}{\|b\|^2}.
$$

这也说明两半平面并集对单个非零参考的任意目标均值均可零残差拟合。证毕。

**扩展与边界。** 该单平面均值关系不提供检验单一几何 SR 核的剩余限制。若允许前后纯未来项非负混合，净系数为

$$
v=\beta\left(\frac{\alpha_F}{D_{\gamma_F}}-
\frac{\alpha_B}{D_{\gamma_B}}\right).
$$

两侧权重均可严格为正而保持任意规定 $v$，再用自由 $c$ 匹配 $u$。所以正符号不证明组分纯前向。限制 $c\ge0$、增益范围或跨块共享参数时，族会缩小，必须重新求相容集合。

### 反例 T07：自由当前项使两学习条件完全换向而观测不变

在任务块 $J,-J$ 中，设 $\gamma=1/2$、$\alpha=5$。前向机制取当前项 $c_F=1$：

$$
G^F_+=-I+4J,\qquad G^F_-=-I-4J.
$$

后向机制取 $c_B=3$：

$$
G^B_+=I-4J=-G^F_+,\qquad
G^B_-=I+4J=-G^F_-.
$$

将未知共同基础映射从 $B$ 改为 $B'=-B$，就有 $G^F_\pm B=G^B_\pm B'$。完整的两个条件模式及相对算子均相同：

$$
T_{\rm rel}=(-15I-8J)/17.
$$

**证明。** $N_{1/2}(J)=(-\tfrac12I+J)/(1+\tfrac14)$，逐项代入并相乘即可。所有当前项和未来增益在例中均为正。证毕。

本反例允许共同映射未知且模型间不同；若独立合格参考已固定同一个 $B$，则不能采用 $B'=-B$，T06 的直接方向结果仍成立。它也不否定“完整 SR 加非负额外当前项”的定理：后者把当前与未来权重约束在不同的连续族中。

## A4. 用互反任务解除参照歧义

### 定理 T08：一般实压缩块的互反符号

设 $A=S+K$ 为实压缩算子，$K\ne0$。两个任务在同一对比坐标中为 $A,A^\top$。若两条件均为完整前向 SR，允许不同折扣与正增益，其精确相对算子为

$$
T=aF_\gamma(A)(I-\eta A^\top),
\qquad a>0,\quad0\le\gamma,\eta<1,\quad\gamma+\eta>0.
$$

则 $\langle T,K\rangle_F>0$。两条件均后向时为负；以同一个对称部分 $S$ 构造的对称机制为零。

**证明。** 记 $F=F_\gamma(A)$，插入 $F^{-1}$ 与 $F^{-\top}$ 得

$$
\operatorname{skew}T/a
=F\{(\gamma+\eta)K-\gamma\eta(SK+KS)\}F^\top.
$$

沿用 T04 中 $H=I-\gamma S$、$L=H^{-1/2}KH^{-1/2}$，令

$$
V=(I-\gamma^2L^2)^{-1},
\qquad Z=-L^2V\succeq0.
$$

$K\ne0$ 保证 $Z\ne0$。因 $S$ 与 $H$ 可交换，花括号可写为 $H^{1/2}L_0H^{1/2}$，其中

$$
L_0=(\gamma+\eta)L-\gamma\eta(SL+LS).
$$

迹的循环性给出

$$
\begin{aligned}
\langle T,K\rangle_F/a
&=-\operatorname{tr}(R_\gamma L_0R_\gamma^\top L)\\
&=-\operatorname{tr}(L_0LV)\\
&=\operatorname{tr}\{[(\gamma+\eta)I-2\gamma\eta S]Z\}.
\end{aligned}
$$

最后一步使用 $L$ 与 $V$ 可交换；没有假设 $S$ 与 $L$ 可交换。又

$$
(\gamma+\eta)I-2\gamma\eta S
\succeq[\gamma(1-\eta)+\eta(1-\gamma)]I\succ0.
$$

正定矩阵与非零半正定矩阵乘积的迹严格为正。后向以 $A^\top=S-K$ 替换 $A$，使标记反号；对称机制的两个因子为同一个对称矩阵的函数，故乘积对称。证毕。

**零边界。** 仅一个折扣为零仍有严格方向；两者均为零时 $T=aI$。$K=0$ 时任务本身没有方向标记。允许折扣、非可逆性或增益趋零时，严格符号不等于统一正的噪声分离距离。该定理也未包括条件间切换方向的 FB、BF 机制；这些机制在 T13 明确进入竞争族。

### 命题 T09：非负额外当前成分的已证扩展

设 $r_0,s_0\ge0$，考虑

$$
T=a[F_\gamma(A)+r_0I][F_\eta(A^\top)+s_0I]^{-1}.
$$

对二维实压缩块以及任意维正规块，T08 的方向符号保持；对满足 $A^2=0$ 的任意实块也保持。

**二维证明。** 正交对角化对称部分后，可写

$$
A=\begin{pmatrix}s_1&k\\-k&s_2\end{pmatrix},\quad
h_i=1-\gamma s_i>0,\quad j_i=1-\eta s_i>0.
$$

令 $d_\gamma=h_1h_2+\gamma^2k^2$、$d_\eta=j_1j_2+\eta^2k^2$，

$$
\Delta=(j_2+s_0d_\eta)(j_1+s_0d_\eta)+\eta^2k^2>0.
$$

显式求逆可得 $\operatorname{skew}T=cK$，其中

$$
c=\frac{a d_\eta}{2d_\gamma\Delta}
\{\eta(h_1+h_2+2r_0d_\gamma)+
\gamma(j_1+j_2+2s_0d_\eta)\}>0.
$$

所有项非负且 $\gamma+\eta>0$，证明成立。正规实矩阵正交分解为一维实块与二维旋转块，逐块加和即可。若 $A^2=0$，直接展开得

$$
\operatorname{skew}T
=a\left\{\frac{\gamma}{1+s_0}
+\frac{(1+r_0)\eta}{(1+s_0)^2}\right\}K.
$$

其余项仅为 $I$ 与 $AA^\top$ 的倍数，因而对称。证毕。此处没有声称任意高维非正规算子、任意独立特征映射或自由当前/未来比均受同一证明覆盖。

## A5. 满秩要求与任务相容测量

### 定理 T10：任意算子功能量的精确可观测条件

设 $X_1=TX_0$，暂不限制 $T$ 的参数形式，令 $\Pi=X_0X_0^+$。数据确定功能量 $\langle T,W\rangle_F$，当且仅当

$$
W=W\Pi.
$$

**证明。** 数据精确确定 $T\Pi=X_1X_0^+$。若 $W=W\Pi$，则 $\langle T\Pi,W\rangle_F=\langle T,W\rangle_F$。反之，取 $\Delta=W(I-\Pi)\ne0$。有 $\Delta X_0=0$，但

$$
\langle\Delta,W\rangle_F=\|W(I-\Pi)\|_F^2>0.
$$

因此 $T$ 与 $T+\Delta$ 观测相同而功能量不同。证毕。

完整任意算子要求 $X_0$ 满行秩；若有 $n$ 个状态并去掉常数方向，线性通道数 $p\ge n-1$ 只是满秩的必要计数条件，并不保证模式线性独立或稳定。恢复指定标记只需标记的相关方向可见。受限参数模型还可进一步减弱要求，如 T12 与 T17；不能把任意算子的必要条件当成所有机制性质的必要条件。

### 命题 T11：共同约化块与八循环的 $P/P^3$ 干预

若 $E^\top E=I$ 且

$$
AE=EA_E,\qquad A^\top E=EA_E^\top,
$$

则对在相关谱上合法的有理函数 $f$ 有 $E^\top f(A)=f(A_E)E^\top$。所以可先取共同任务块再分析 SR，而不引入其余块的未知作用。

**证明。** 两个关系表示子空间及正交补均保持不变。$A$ 在相应正交基中为块对角矩阵；求幂、求可逆矩阵的逆和取有理函数均保持相同分块，结论成立。证毕。任意状态投影不自动具有该性质。

对于八状态行循环 $P_{j,j+1\bmod8}=1$、$Q=P^3$，定义

$$
C_{1j}=\tfrac12\cos(\pi j/2),\qquad
C_{2j}=-\tfrac12\sin(\pi j/2),\qquad j=0,\ldots,7.
$$

直接相乘有

$$
C\mathbf1=0,\quad CC^\top=I,\quad
CP=JC,\quad CQ=-JC,\quad
CP^\top=-JC,\quad CQ^\top=JC.
$$

因此整图虽非互为反转，$Q\ne P^\top$，这个预定平面上却严格互反。对于完整 SR，

$$
T_E=a(I-\gamma J)^{-1}(I+\eta J)
=\frac{a}{1+\gamma^2}
\{(1-\gamma\eta)I+(\gamma+\eta)J\}.
$$

完整标记

$$
W=\tfrac12(P-P^\top-Q+Q^\top)=C^\top(2J)C
$$

的秩为二、$\|W\|_F^2=8$，且

$$
\frac{\langle a(I-\gamma P)^{-1}(I-\eta Q),W\rangle_F}{8}
=\frac{a(\gamma+\eta)}{2(1+\gamma^2)}.
$$

**恒等式验证。** 以循环幂展开 $(I-\gamma P)^{-1}=(\sum_{j=0}^7\gamma^jP^j)/(1-\gamma^8)$。各幂在 Frobenius 内积下正交，标记只提取第 1、7、3、5 次幂系数。合并分子得到 $(\gamma+\eta)(1-\gamma^2+\gamma^4-\gamma^6)$，约去 $(1-\gamma^2)(1+\gamma^4)$ 即得公式。$W=C^\top(2J)C$ 亦可逐项验证。后向转置使符号反转；两个可交换的对称化转移给零分数。证毕。

该命题是任务结构结论，不宣称任何已采集参与者的共同映射、可见性或噪声条件已经合格。

### 定理 T12：任务相容读出与弱可见性

对训练后冻结的非零模板 $R\in\mathbb R^{2\times p}$，定义

$$
z_R(X)=\frac1{\|R\|_F}
\begin{pmatrix}\langle R,X\rangle_F\\
\langle JR,X\rangle_F\end{pmatrix}.
$$

则对所有 $X,u,v$ 有

$$
z_R(JX)=Jz_R(X),\qquad
z_R((uI+vJ)X)=(uI+vJ)z_R(X).
$$

若两个总体模式满足 $X_+=T X_-$、$T=uI+vJ$，则压缩均值 $x_\pm=z_R(X_\pm)$ 满足相同关系，并且

$$
\det[x_-,x_+]=v\|x_-\|^2.
$$

所以只要真实参考投影 $x_-\ne0$，便可识别 $v$ 的符号；不需要两行参考满秩。

**证明。** $J^\top=-J$ 给出 $\langle R,JX\rangle_F=-\langle JR,X\rangle_F$、$\langle JR,JX\rangle_F=\langle R,X\rangle_F$，即交换关系。又 $R,JR$ 正交等长，故 $\|z_R(X)\|\le\|X\|_F$。行列式公式来自 $\det[x,Jx]=\|x\|^2$。证毕。

这是一项结构恒等式，不要求训练模板等于真实参考；但模板非零不保证正式均值有非零投影。落在 $\operatorname{span}\{R,JR\}$ 正交补的信号会丢失。任意线性压缩 $L$ 若要保持候选作用，须另验 $LT_\theta=\widetilde T_\theta L$；方差最大化的状态 PCA 或分类器不自动满足这一条件。

**一般旋转块的两学习条件接口。** 设 $A=\rho I+\kappa J$，$\rho^2+\kappa^2\le1$、$\kappa>0$。纯前向相对算子有

$$
u=\frac{a\{(1-\gamma\rho)(1-\eta\rho)-\gamma\eta\kappa^2\}}
{(1-\gamma\rho)^2+\gamma^2\kappa^2},
$$

$$
v=\frac{a\kappa(\gamma+\eta-2\gamma\eta\rho)}
{(1-\gamma\rho)^2+\gamma^2\kappa^2}>0.
$$

正性由 $\gamma(1-\eta\rho)+\eta(1-\gamma\rho)>0$ 得到。一般块的 $u$ 可以为负；只有进一步限制为下节 quarter-turn 块，才有特定的正 $u$ 扇区。该接口允许参考本身已学习，但不能跨越 T07 中自由当前项造成的等价。

## A6. 四种机制族、混合与部分定量识别

### 定理 T13：quarter-turn 相对连续族的精确角度像

在任务块 $J,-J$ 中，两个条件均使用纯完整 SR，各自折扣在 $[0,g]$，其中 $g<1$，两条件正增益比任意。以 F、B 分别表示相对于各自任务的前、后向机制。令 $h=\arctan g$，将相对算子写为

$$
T=r(\cos\phi\,I+\sin\phi\,J),\qquad r>0.
$$

则四族的准确角度像为

| 机制标签 | 相对角 $\phi$ 的全部范围 |
|---|---|
| FF | $[0,2h]$ |
| BB | $[-2h,0]$ |
| FB | $[-h,h]$ |
| BF | $[-h,h]$ |

每个范围中的每个角度和任意正径向尺度均可实现。

**证明。** 单个完整 SR 满足

$$
F_\gamma(J)=\frac{I+\gamma J}{1+\gamma^2}
=\frac1{\sqrt{1+\gamma^2}}
\{\cos\theta_\gamma I+\sin\theta_\gamma J\},
\qquad\theta_\gamma=\arctan\gamma.
$$

反任务上的前向角为 $-\theta_\eta$；相对算子的角是目标减参考，故四族分别为

$$
\theta_\gamma+\theta_\eta,\quad
-\theta_\gamma-\theta_\eta,\quad
\theta_\gamma-\theta_\eta,\quad
-\theta_\gamma+\theta_\eta.
$$

独立角度取值于 $[0,h]$，即得区间。FF 中取 $\theta_\gamma=\theta_\eta=\phi/2$ 实现每个角；切换族中取一个角为零、另一个取所需绝对值。增益比自由可匹配任意 $r>0$。证毕。

FF 的角与尺度具体为

$$
\phi=\arctan\gamma+\arctan\eta,\qquad
r=a\sqrt{\frac{1+\eta^2}{1+\gamma^2}}.
$$

因此 $0<\phi\le h$ 可由 FF 与切换机制同时产生，正面积不能单独识别 FF；若数据相容集合非空且全部满足 $\phi>h$，则只有 FF 相容。其可观测形式为 $u>0$ 且 $v-gu>0$。零折扣和净方向零边界必须保留。若增益有额外上下界，径向范围不能继续当作全正轴。

当折扣范围为 $[0,1)$ 时，FF、BB 角区间分别为 $[0,\pi/2)$、$(-\pi/2,0]$，切换族为 $(-\pi/4,\pi/4)$。因此精确无噪声的 $\phi=\pi/4$ 已不属于切换族，但它属于该族闭包；带稳健余量的排除仍要求置信范围严格超过 $\pi/4$。

### 命题 T14：非负完整 SR 混合保留净方向扇区

在各自任务坐标中，设一个条件允许

$$
H(J)=cI+\int_0^g F_\gamma(J)\,\mu_F(d\gamma)
+\int_0^g F_\gamma(-J)\,\mu_B(d\gamma),
$$

其中 $c\ge0$，$\mu_F,\mu_B$ 为有限非负测度，且该算子非零；所有尺度、方向共用同一个特征映射。则

$$
H(J)=uI+vJ,\quad
u=c+\int_0^g\frac{\mu_F+\mu_B}{1+\gamma^2}>0,
$$

$$
v=\int_0^g\frac{\gamma(\mu_F-\mu_B)}{1+\gamma^2},
\qquad |v|\le gu.
$$

故每个条件净相位 $\theta=\arctan(v/u)$ 均在 $[-h,h]$，按两个条件净相位的正负分类仍得到 T13 的四个扇区。若相对角大于 $h$，两条件相对于各自任务都必须净前向。

**证明。** 展开每个完整 SR 得系数；逐项用 $\gamma\le g$ 和非负权重得到不等式。两条件相对角是各自任务坐标净相位之和，区间运算同 T13。证毕。

此结论识别的是加权净方向，而非没有任何相反方向组分。它不涵盖任意自由 $cI+\alpha N_\gamma$ 的相对族；后者可有 $u<0$，T07 已给反例。

### 定理 T15：从方向到折扣的尖锐部分识别

在纯完整 SR 的 FF 模型中，设 $\gamma,\eta\in[0,g]$、$h=\arctan g$，给定精确相对角 $\phi\in[0,2h]$，且增益比自由。则

$$
\gamma\in
\left[
\tan\{\max(0,\phi-h)\},
\tan\{\min(h,\phi)\}
\right],
$$

$\eta$ 具有相同边缘范围，联合相容集为 $\arctan\gamma+\arctan\eta=\phi$。每个允许值均可实现，因而这些界是尖锐的。

**证明。** 令 $\theta_1=\arctan\gamma$、$\theta_2=\arctan\eta$。约束为 $0\le\theta_i\le h$、$\theta_1+\theta_2=\phi$，等价于

$$
\max(0,\phi-h)\le\theta_1\le\min(h,\phi).
$$

$\tan$ 严格递增，得到边缘区间。对区间内任取 $\theta_1$，令 $\theta_2=\phi-\theta_1$，再按 T13 选择增益匹配观测尺度，即达到该值。证毕。

若原范围为 $[0,1)$，则有三个精确情形：

$$
\begin{array}{ll}
0\le\phi<\pi/4:
&\gamma,\eta\text{ 的边缘识别集为 }[0,\tan\phi],\\
\phi=\pi/4:
&\gamma,\eta\text{ 的边缘识别集为 }(0,1),\\
\pi/4<\phi<\pi/2:
&\gamma,\eta\text{ 的边缘识别集为 }
(\tan(\phi-\pi/4),1).
\end{array}
$$

特别是

$$
\boxed{\phi>\pi/4\quad\Longrightarrow\quad
\gamma>\tan(\phi-\pi/4),\quad
\eta>\tan(\phi-\pi/4).}
$$

例如 $\phi=70^\circ$ 强制两个折扣均大于 $\tan25^\circ\approx0.4663$，但通常不唯一决定任一个折扣。有限上界 $g$ 下，$\phi=0$ 或 $2h$ 等边界可使两折扣唯一；“一般不可单独恢复”不能替代这些边界事实。

对 T14 的混合模型，同样的相位逻辑只约束有效净系数 $q_i=v_i/u_i$，不能逐个约束混合中的每个 $\gamma$。有限噪声下应投影整个合法联合相容集合；不能把点估计角代入上述公式后直接称作置信界。

### 命题 T16：单旋转平面不区分同方向多尺度与单一折扣

设 $H=cI+\int F_\gamma(J)\mu(d\gamma)$，$c\ge0$、$\mu\ge0$、$H\ne0$。写 $H=uI+vJ$。则

$$
\gamma_{\rm eff}=v/u,\qquad
g_{\rm eff}=u(1+\gamma_{\rm eff}^2)
$$

满足 $H=g_{\rm eff}F_{\gamma_{\rm eff}}(J)$。若所有折扣小于一，则 $0\le\gamma_{\rm eff}<1$。

**证明。** $u>0$、$v\ge0$，且 $v/u$ 为以 $1/(1+\gamma^2)$ 加权的折扣平均，并包含当前项的零折扣质量，故在合法范围内。直接代入 $F_{\gamma_{\rm eff}}(J)=(I+\gamma_{\rm eff}J)/(1+\gamma_{\rm eff}^2)$ 即得。证毕。

更多同平面模板可以检验共享线性算子，但不能拆开已在该平面完全相同的时间核。定量出口需要额外传播形状或额外任务块。

## A7. 识别折扣与时间核的建设性出口

### 定理 T17：三个传播形状是自由当前项模型的精确折扣条件

设 $P\in\mathbb R^{m\times m}$ 已知且谱半径不超过一，直接参考 $B\in\mathbb R^{m\times p}$ 已知。考虑完整参数域

$$
Y(c,\alpha,\gamma)
=cB+\alpha P(I-\gamma P)^{-1}B,
\qquad c\in\mathbb R,\quad\alpha>0,\quad0\le\gamma<1.
$$

则以下命题等价：

1. $\gamma$ 在该族内全局结构可辨识。
2. $c,\alpha,\gamma$ 三个参数在该族内共同可辨识。
3. $B,PB,P^2B$ 作为同维矩阵在实标量线性组合意义下线性独立。

若第三条不成立，任意给定合法参数点的相同观测都可由任意另一个合法折扣产生；不只是存在一对例外碰撞。

**充分性证明。** 假设两组参数产生相同 $Y$，记 $\delta$ 为第二组折扣，$c',\alpha'$ 为其余参数。所有矩阵均为 $P$ 的函数，故左乘 $(I-\gamma P)(I-\delta P)$ 得

$$
\begin{aligned}
0={}&(c-c')B\\
&+\{\alpha-\alpha'-(c-c')(\gamma+\delta)\}PB\\
&+\{(c-c')\gamma\delta-\alpha\delta+\alpha'\gamma\}P^2B.
\end{aligned}
$$

三个传播形状线性独立，系数分别为零。第一式给 $c=c'$，第二式给 $\alpha=\alpha'$，第三式变为 $\alpha(\gamma-\delta)=0$。因 $\alpha>0$，得 $\gamma=\delta$，证明共同唯一性。

**必要性证明。** 若 $B=0$，所有观测为零。若 $B\ne0$ 而 $PB=\lambda B$，则

$$
Y=\left(c+\frac{\alpha\lambda}{1-\gamma\lambda}\right)B.
$$

任意改变 $\gamma$ 后，都可保持正 $\alpha$ 并改变自由 $c$ 使观测不变。

余下情况为 $B,PB$ 线性独立，但存在实数 $s,t$ 使

$$
P^2B=sB+tPB.
$$

于是 $\operatorname{span}\{B,PB\}$ 对左乘 $P$ 不变。在有序基 $(B,PB)$ 上，左乘 $P$ 的矩阵为

$$
C=\begin{pmatrix}0&s\\1&t\end{pmatrix}.
$$

其谱是左乘 $P$ 在该子空间的谱；故 $I-\gamma C$ 对 $0\le\gamma<1$ 可逆。令

$$
D_\gamma=\det(I-\gamma C)=1-t\gamma-s\gamma^2.
$$

$D_\gamma$ 在整个连通区间中连续、不为零，且 $D_0=1$，所以始终严格为正。直接求逆得

$$
P(I-\gamma P)^{-1}B
=\frac{\gamma s B+PB}{D_\gamma},
$$

从而

$$
Y=\left(c+\frac{\alpha\gamma s}{D_\gamma}\right)B
+\frac{\alpha}{D_\gamma}PB.
$$

对任意另一个折扣 $\delta\in[0,1)$，取

$$
\alpha'=\alpha D_\delta/D_\gamma>0,
\qquad
c'=c+\alpha s(\gamma-\delta)/D_\gamma.
$$

代入即有 $Y(c',\alpha',\delta)=Y(c,\alpha,\gamma)$，证明必要性和逐点别名。证毕。

**解释与范围。** 这里的线性独立是将三个矩阵向量化后检查三列的秩，不是要求 $B$ 行满秩。比如 $J_3$ 为三阶下移位块、$B=e_1$，则

$$
B=e_1,\quad J_3B=e_2,\quad J_3^2B=e_3,\qquad
Y=(c,\alpha,\alpha\gamma)^\top.
$$

一列通道即可辨认三个参数，因为任务让三个传播形状分开。T19 将这一块嵌入合法随机任务；本例仅解释代数条件，不代表神经实验已成功。

本定理假设 $P$ 和直接 $B$ 已知、$\alpha$ 严格为正、$c$ 全实且无额外参数限制；改变这些条件需重新判断。若 $\alpha=0$ 也允许，则该边界所有 $\gamma$ 均同像。若 $P$ 或参考未知，不能借本定理绕过 T05 或 T07。它是在指定方向分支内的参数识别，不独自识别机制方向。线性独立也不是有限噪声下稳定估计的保证；近相关的传播形状会使置信区间宽。

### 定理 T18：两个不同旋转块的显式折扣恢复

采用 T17 的直接参考模型，设可见的两个任务块为

$$
A_j=a_jI+\sigma\beta_jJ,\quad
a_j^2+\beta_j^2=1,\quad\beta_j>0,\quad a_1\ne a_2,
$$

且两个块共享同一 $c,\alpha,\gamma$，$\alpha>0$。每个参考均有非零相容投影，所以可观测 $u_j,v_j$。有

$$
D_j=1-2a_j\gamma+\gamma^2,\quad
u_j=c+\alpha(a_j-\gamma)/D_j,\quad
v_j=\sigma\alpha\beta_j/D_j.
$$

令

$$
q=\frac{v_1/\beta_1}{v_2/\beta_2}=\frac{D_2}{D_1},
\qquad
t=\frac{q-1}{2(qa_1-a_2)}.
$$

则合法模型中 $t=\gamma/(1+\gamma^2)$ 且

$$
\boxed{\gamma=\frac{2t}{1+\sqrt{1-4t^2}}.}
$$

随后

$$
\alpha=\frac{\sigma v_jD_j}{\beta_j},\qquad
c=u_j-\alpha(a_j-\gamma)/D_j.
$$

**证明。** $q=(1-2a_2t)/(1-2a_1t)$，直接求解得 $t$。又

$$
\frac{dq}{d\gamma}
=\frac{2(a_1-a_2)(1-\gamma^2)}{D_1^2}
$$

在 $[0,1)$ 上符号严格固定，所以一一对应。$t$ 严格递增于 $[0,1/2)$，解其二次方程并选小于一的根，得所述逆式；写成有理化形式避免 $\gamma=0$ 时的 $0/0$。再用两块系数恢复其余参数。证毕。

**剩余可检验约束。** 两块各自算出的 $c$ 必须一致；这使四个可观测系数承担三个共享参数，增加一个均值层面的可失败限制。例如八循环的两个实际 Fourier 块可取 $(a_1,\beta_1)=(0,1)$、$(a_2,\beta_2)=(\sqrt2/2,\sqrt2/2)$。若 $v_1=4/5$、$v_2=(8+10\sqrt2)/17$，则强制 $\gamma=1/2,\alpha=1$；若又观察 $u_1=u_2=1$，两块给 $c_1=7/5,c_2=(19-6\sqrt2)/17$，模型不相容。这里明确选择有限随机任务中可实现的单位圆块；任意抽象正交块并不自动对应有限马尔可夫链的单位模特征值。

也可令 $s_j=\sigma v_j/\beta_j$，在 $s_1\ne s_2$ 时有

$$
\gamma=\frac{s_1a_1-s_2a_2-(u_1-u_2)}{s_1-s_2}.
$$

该式在 $\gamma=0$ 时退化，而比值逆式仍可用。比值对 $\gamma$ 的导数在 $\gamma\to1$ 时趋零，只说明该比值的敏感性变差，不证明完整系数信息必同样退化。两个块若可独立选择增益、当前项或折扣，以上共同参数结论不再成立。统计估计应联合反演系数和参考误差，而非对不稳定比值直接套正态近似。

### 命题 T19：合法随机任务中的可扩展幂零块

以二元字 $x=(x_1,\ldots,x_k)$ 为状态，每次删除第一位并在末尾加入独立等概率二元位。完整图有 $n=2^k$ 个状态，每状态两条概率为 $1/2$ 的入边与出边；故为双随机矩阵。定义

$$
C_{j,x}=n^{-1/2}(1-2x_j),\qquad j=1,\ldots,k.
$$

有

$$
C\mathbf1=0,\quad CC^\top=I_k,\quad
CP=J_kC,\quad CP^\top=J_k^\top C,\quad J_k^k=0,
$$

其中 $(J_k)_{i+1,i}=1$。

**证明。** 各二元位在全部状态中独立平衡，得到零均值和正交性。固定终点后，前驱第一位是独立等概率变量，其余位由终点前 $k-1$ 位确定；故 $CP$ 第一行零，第 $j>1$ 行等于 $C$ 的第 $j-1$ 行。对后继求平均得到转置关系。移位 $k$ 次为零。证毕。

幂零的是预定对比块，不是完整转移矩阵；完整任务仍为正常随机转移系统。该构造使用经典 de Bruijn 图结构，本文用其获得有限滞后上可分开的传播方向。

### 定理 T20：直接参考下的有限时间核与混合检验

在 T19 块中，若合格直接参考支持所需算子量的恢复，目标时间核为 $f(z)=\sum_h w_hz^h$，则

$$
T=a\sum_{h=0}^{k-1}w_hJ_k^h+cI.
$$

定义

$$
\psi_h=\frac1{k-h}\sum_{i=1}^{k-h}
(T_{i+h,i}-T_{i,i+h}),\qquad 1\le h<k.
$$

则 $\psi_h=aw_h$；后向核给 $-aw_h$；对称算子给零。完整单一 SR 的 $w_h=\gamma^h$，故在有非零分母时 $\psi_{h+1}/\psi_h=\gamma$，且 $k\ge4$ 时有 $\psi_2^2=\psi_1\psi_3$。

对具有共同映射的非负完整 SR 时域混合 $w_h=\int_0^1\gamma^h\mu(d\gamma)$，有

$$
w_2^2\le w_1w_3.
$$

当 $w_1>0$，令 $\nu(d\gamma)=\gamma\mu(d\gamma)/w_1$，则

$$
\frac{w_1w_3-w_2^2}{w_1^2}
=\operatorname{Var}_{\nu}(\gamma).
$$

**证明。** $J_k^h$ 仅在第 $h$ 条下副对角线上为一，逐项得到 $\psi_h$。将 Cauchy–Schwarz 应用于 $\gamma^{1/2},\gamma^{3/2}$ 得不等式；归一化后 $\mathbb E_\nu\gamma=w_2/w_1$、$\mathbb E_\nu\gamma^2=w_3/w_1$，给方差恒等式。证毕。

方差也可用 $\psi_h$ 表示，整体正增益消去。零折扣质量对正滞后不可见，等式不能排除它。有限可见滞后既不能验证无限尾部，也不能唯一确定学习算法。若各尺度具有任意不同的特征映射，混合矩不等式不再自动适用。

### 定理 T21：未知反转参考下的时间核与双折扣恢复

在同一幂零块内，令目标和反向参考核分别为

$$
f(J_k)=\sum_{h=0}^{k-1}f_hJ_k^h,\qquad
q(J_k^\top)=\sum_{h=0}^{k-1}q_h(J_k^\top)^h,
\qquad q_0\ne0.
$$

若所需相对算子量可观测，则

$$
T=a f(J_k)q(J_k^\top)^{-1}
\quad\Longrightarrow\quad
T_{h+1,1}=\frac{a}{q_0}f_h.
$$

**证明。** $J_k^\top e_1=0$，所以 $q(J_k^\top)e_1=q_0e_1$，逆矩阵将 $e_1$ 映为 $e_1/q_0$。再乘 $f(J_k)$ 得第一列公式。证毕。

这在一个共同尺度之外恢复目标的可见核，不要求参考为几何核；方向解释另需 $a/q_0>0$。当前项只改变 $f_0$ 或 $q_0$，所以正滞后比例保持。$q_0=0$ 时参考块奇异，不支持这个论证。

纯完整 SR 给出

$$
T=a(I-\gamma J_k)^{-1}(I-\eta J_k^\top),
\qquad a=T_{11},\quad
\gamma=T_{21}/T_{11},\quad
\eta=-T_{12}/T_{11}.
$$

这些式子在 $k\ge2$ 时成立。$k\ge4$ 时第一列还须满足 $T_{31}^2=T_{21}T_{41}$。参数识别以纯 SR 族为条件，其他矩阵元素和竞争核仍须比较。

不能将 T20 的直接副对角线差量误用于这个相对算子：上下副对角线均受参考核影响。精确单一 SR 在 $\gamma=0.6,\eta=0.8,a=1,k=4$ 下，误用会得到 $\psi_2^2-\psi_1\psi_3=-0.18607104$，错误地似乎支持多尺度混合。相对第一列是适用于该参考问题的量。其可观测性也应按 T10 或受限模型另证，不能自动由二维方向读出推得。

### 命题 T22：零方向不是唯一对称机制

在 $J_4$ 块上，混合两方向的完整地图

$$
\tfrac12\{F_\gamma(J_4)+F_\gamma(J_4^\top)\}
$$

具有常数对角线。而先对称化转移再取 SR 的

$$
F_\eta((J_4+J_4^\top)/2)
$$

在 $0<\eta<1$ 时满足

$$
[F_\eta]_{22}-[F_\eta]_{11}
=\frac{x^2}{1-3x^2+x^4}>0,\qquad x=\eta/2.
$$

**证明。** 幂零完整地图的所有对角线均为一，前一个结论立即成立。后一矩阵的三对角行列式为 $1-3x^2+x^4$，两个相关余子式为 $1-x^2$ 与 $1-2x^2$，相减即得公式。证毕。正增益保留差异，额外标量当前项消去；只有零反对称分数无法辨别这两种机制。

## A8. 时间历史与误差预算的边界

### 定理 T23：记录历史的后瞻别名与独立探测控制

若真实瞬时响应只有当前状态码 $b(S_t)$，但测量为

$$
y_t=\sum_{h\ge0}\lambda^h b(S_{t-h}),\qquad0\le\lambda<1,
$$

在平稳马尔可夫链下有

$$
\mathbb E[y_t\mid S_t=i]=[M_\lambda(P^*)B]_i.
$$

**证明。** 条件过去分布为 $\Pr(S_{t-h}=j\mid S_t=i)=[(P^*)^h]_{ij}$，逐项取条件期望并求和。证毕。因而图反转也会反转该测量伪影；任务互反本身不排除历史滤波机制。

若学习后的探测状态 $Z_t$ 独立同分布，并独立于学习序列，已学习响应 $r_c(z)$ 在测量期稳定，且记录滤波器绝对可和，则

$$
y_t=\sum_{h\ge0}h_h r_c(Z_{t-h})
\quad\Longrightarrow\quad
\mathbb E[y_t\mid Z_t=z]
=h_0r_c(z)+\sum_{h\ge1}h_h\mathbb E r_c(Z).
$$

后项与当前状态无关，故状态对比将其消去。这是条件均值层面的控制。若需要方向符号，当前有效增益还须与观测模型规定一致；零系数会毁掉可见性。不放回排列、禁止相邻重复、状态依赖适应、继续学习或不稳定测量均不满足这一简单证明。

### 定理 T24：完整算子恢复的确定性误差界

设

$$
\widehat X_0=X_0+E_0,\quad
\widehat X_1=TX_0+E_1,\quad
\|E_i\|_F\le e_i.
$$

若 $\widehat X_0$ 满行秩且 $s=\sigma_{\min}(\widehat X_0)>e_0$，定义 $\widehat T=\widehat X_1\widehat X_0^+$。则

$$
\|\widehat T-T\|_F\le
\delta_F=\frac{e_1+\|\widehat T\|_2e_0}{s-e_0}.
$$

**证明。** 满行秩给出

$$
\widehat T-T=(E_1-TE_0)\widehat X_0^+.
$$

令 $d=\|\widehat T-T\|_F$，则 $d\le(e_1+\|T\|_2e_0)/s$，同时 $\|T\|_2\le\|\widehat T\|_2+d$。整理得到 $d(s-e_0)\le e_1+\|\widehat T\|_2e_0$。证毕。无需假设 $\widehat X_1-\widehat T\widehat X_0=0$。

对预定标记 $W$，有

$$
|\langle\widehat T-T,W\rangle_F|\le\delta_F\|W\|_F.
$$

因此 $\omega$ 的误差半径为 $\delta_F/\|K\|_F$，八循环归一化标记半径为 $\delta_F/\sqrt8$，T20 的 $\psi_h$ 半径为 $\delta_F\sqrt{2/(k-h)}$。仅当合法误差区间排除竞争符号才声明方向；秩不足不是零方向。

该结论以真实扰动满足界为条件；只有输入误差界具有相应联合覆盖时才能继承统计解释。截断伪逆、参考漂移、任务转移不确定性和系统性模型失配不会自动包含在 $e_i$ 内。浮点容差或条件数保护也不等于经过向外舍入验证的机器证明。

## A9. 连续族分离、联合误差与拒判

### 命题 T25：全局最坏情况距离与数据局部排除

在预定观测、尺度和参数范围下，令 $\mathcal M_F,\mathcal M_B$ 为两类非空紧致预测集合，并定义

$$
d_{FB}=\inf_{f\in\mathcal M_F,\;b\in\mathcal M_B}\|f-b\|.
$$

若两类允许的扰动是完整闭范数球，半径分别为 $\varepsilon_F,\varepsilon_B$，则两误差邻域不相交当且仅当

$$
d_{FB}>\varepsilon_F+\varepsilon_B.
$$

**证明。** 若同一点同时位于两误差邻域，三角不等式给出一对预测的距离不超过半径之和。反之，紧致性保证最小距离由 $f,b$ 取得；若该距离不超过半径之和，在线段 $fb$ 上可选一点同时距两端不超过各自半径，故邻域相交。证毕。非紧集下，下确界未取得时的等号需要另外处理。

这个数是设计的最坏情况承诺。若有参数可趋向零折扣、零方向或零幅度，$d_{FB}$ 可能为零，但某一真实观测邻域仍可能完全排除一个竞争族。局部推断应检查同一个有效不确定性区域 $\mathcal C$ 与各族的交，而不要求两个整个族处处正分离。若参考 $B$ 有误差，将 $\widehat B$ 直接代入条件距离也不会自动产生覆盖保证。

**算法相关的充分条件。** 若真实点到错误族距离为 $g$，实际观测误差为 $e$，最近整个族的正确判决由 $g>2\|e\|$ 保证。若正确类只使用训练后冻结的预测，且其距真值为 $h$，则充分条件变为 $g>2\|e\|+h$。

证明是：错误族到观测至少 $g-\|e\|$，正确真点到观测最多 $\|e\|$，冻结正确点则最多 $\|e\|+h$。所以不能用整族距离保证为一个有迁移误差的冻结点判决背书。以上均是确定性几何，不是额外创造的抽样概率。

### 定理 T26：候选均值残差反演与单一属性声明的错误控制

给定独立训练资料后，固定状态标签、任务块、通道处理、模板、ROI、时间窗和完整候选空间。设正式观测由 $n$ 个完整配对向量

$$
Z_i=(B_i^\top,Y_i^\top)^\top
$$

构成，并在条件训练分布下独立同分布高斯。对于任意预定候选 $\theta$，定义 $q$ 维残差

$$
E_i(\theta)=Y_i-T_\theta B_i.
$$

设真实 $\theta_0$ 属于候选空间，满足 $\mathbb E Y_i=T_{\theta_0}\mathbb E B_i$，其残差协方差正定，且 $n>q$。令 $S_{E,\theta}$ 为无偏样本协方差，

$$
R(\theta)=\bar E_\theta^\top S_{E,\theta}^{-1}\bar E_\theta,
\qquad
K_q(\alpha)=\frac{q(n-1)}{n(n-q)}F_{q,n-q}(1-\alpha).
$$

将全部连续候选反演为

$$
\widehat\Theta_\alpha
=\{\theta:R(\theta)\le K_q(\alpha)\}.
$$

则真实参数的保留概率为 $1-\alpha$。若只在该集合非空且所有候选的预定属性 $\ell(\theta)$ 均相同时声明该属性，则

$$
\Pr_{\theta_0}\{\text{作出声明且属性错误}\mid D_{\rm train}\}
\le\alpha.
$$

**证明。** 真实残差独立同分布高斯，均值为零。经典 Hotelling 均值分布给出

$$
\frac{n(n-q)}{q(n-1)}R(\theta_0)\sim F_{q,n-q}.
$$

故 $\Pr(\theta_0\in\widehat\Theta_\alpha\mid D_{\rm train})=1-\alpha$。若真实参数保留，所有相容候选共有的属性必与真实属性相同。因此错误声明事件包含于真实参数被排除事件，其概率不超过 $\alpha$。证毕。

这一统计分布是经典结果；本文使用它把参考和目标的共同误差接到 SR 连续族与机制标签上。对 FF、BB、FB、BF 应以带标签的不交并参数空间反演；不同标签即使映为同一 $T$，也必须全部保留。证明不要求按候选族数或无穷参数点数再做 Bonferroni，因为控制的是一个固定真参数被误排除的事件。事后选择窗口、ROI、模板、状态重排或参考来源属于另一类选择，仍须独立冻结或适当调整。

该概率不是 $\Pr(\text{错误}\mid\text{已经声明})$，不是候选的后验概率，也不是同时覆盖全部观测等价参数的保证。它以真实机制和共同参考关系属于候选空间为条件；模型遗漏或系统漂移不能由名义 $p$ 值修复。

**参考—目标共同协方差。** 单模板时 $B_i,Y_i\in\mathbb R^2$。若

$$
S_Z=\begin{pmatrix}S_B&S_{BY}\\S_{YB}&S_Y\end{pmatrix},
$$

则

$$
\bar E_T=\bar Y-T\bar B,\qquad
V_T=S_Y+TS_BT^\top-TS_{BY}-S_{YB}T^\top.
$$

这正是残差的样本协方差；均值残差协方差另除以 $n$。配对内相关完全允许，配对之间相关不由该定理覆盖。

**多个模板。** 对 $r$ 个固定二维模板，将所有参考、目标保留在同一个实际重复单位中。共同算子残差通常维数为 $q=2r$，使用 $I_r\otimes T_\theta$ 和完整交叉协方差，阈值为 $K_{2r}$。模板不是额外独立被试。已知确定冗余可预先消去；不能任意截断评价协方差再沿用未经推导的普通 $F$ 阈值。

### 命题 T27：二维残差与四维联合均值椭球不能混同

单模板且 $S_Z\succ0$ 时，令 $L_T=(-T,I)$。则

$$
\inf_{\mu:\,\mu_Y=T\mu_B}
(\bar Z-\mu)^\top S_Z^{-1}(\bar Z-\mu)
=
(\bar Y-T\bar B)^\top
(L_TS_ZL_T^\top)^{-1}
(\bar Y-T\bar B)
=R(T).
$$

**证明。** 令 $h=\bar Z-\mu$，约束为 $L_Th=L_T\bar Z$。对严格凸二次型作拉格朗日优化，最小化向量为

$$
h_*=S_ZL_T^\top(L_TS_ZL_T^\top)^{-1}L_T\bar Z.
$$

代回得到等式。证毕。

因此二维残差阈值

$$
K_2(\alpha)=\frac{2(n-1)}{n(n-2)}F_{2,n-2}(1-\alpha)
$$

可以通过一个四维**辅助**椭球的相交计算实现，但该辅助椭球不是四维均值的 $1-\alpha$ 置信域。真正的联合均值置信域是

$$
\mathcal C_{1-\alpha}
=\{\mu:(\bar Z-\mu)^\top S_Z^{-1}(\bar Z-\mu)\le K_4(\alpha)\},
$$

$$
K_4(\alpha)=\frac{4(n-1)}{n(n-4)}F_{4,n-4}(1-\alpha),
\qquad n>4.
$$

其覆盖率由四维 Hotelling 分布保证。将同一个 $\mathcal C$ 与所有候选均值关系集合相交，可同时保护由真实均值关系推出的性质；任意额外均值性质或参考漂移范围若要沿用这条区域路线，必须使用 $K_4$ 或另一个真实有效的联合区域。

二维候选残差理论只需 $n>2$ 及残差协方差可逆；采用四维正定辅助几何的具体实现需 $n>4$，不可混写。作为纯阈值算例，$n=17,\alpha=0.05$ 时 $K_2\approx0.4621$、$K_4\approx0.9206$；该算例不认证任何实际 17 个配对的独立性。

### 命题 T28：有向面积偏差与有界参考漂移敏感性

设单次配对 $B_i=b+\varepsilon_{B,i}$、$Y_i=y+\varepsilon_{Y,i}$，噪声均值为零，交叉协方差为 $\Sigma_{BY}$。则

$$
\mathbb E\det(B_i,Y_i)
=\det(b,y)+(\Sigma_{BY})_{12}-(\Sigma_{BY})_{21}.
$$

若 $n$ 次配对 iid，

$$
\mathbb E\det(\bar B,\bar Y)
=\det(b,y)+\frac{(\Sigma_{BY})_{12}-(\Sigma_{BY})_{21}}n.
$$

**证明。** 展开行列式，线性噪声项期望为零，交叉二次项分别为两个交叉协方差元素；样本均值的交叉协方差再除以 $n$。证毕。因此逐次面积平均或其 plug-in 均值不自动是无偏机制指标；T26 直接检验总体均值残差，避免把面积乘积当作无偏的简单单变量量。

若实际目标参考为 $b+\delta$，外部依据保证 $\|\delta\|\le\rho$，则

$$
\inf_{\|\delta\|\le\rho}\det[b+\delta,y]
=\det[b,y]-\rho\|y\|.
$$

对 $y\ne0$，将方向推至零边界的最小参考变化为

$$
\rho_*=\frac{|\det[b,y]|}{\|y\|}.
$$

**证明。** $\det[\delta,y]=\delta^\top(y_2,-y_1)^\top$，其系数范数为 $\|y\|$。Cauchy–Schwarz 给极值，取 $\delta$ 与系数向量反向即可达到；解零行列式的最小范数线性约束得到第二式。证毕。

若采样误差也存在，可在真正联合均值置信域 $\mathcal C$ 上检验

$$
\inf_{(b,y)\in\mathcal C}
\{\det[b,y]-\rho\|y\|\}>0.
$$

当区域覆盖和外部漂移界同时有效时，该条件保证允许漂移下的正方向。若 $\rho$ 只是人为扫描的一组数值，输出只是敏感性分析；它不能证明实际漂移已被涵盖，也不能用 T27 的 $K_2$ 辅助椭球冒充联合域。

## A10. 论文中每种结论的正确用途

| 要回答的问题 | 核心结果 | 实施中不可省略的条件 |
|---|---|---|
| 规范距离能否决定前后向？ | T02–T03 | 明确距离、度量、状态集合与读出；同折扣比较 |
| 完整相对算子拟合能否自动解除参考歧义？ | T04、T07 | 已学习参考和自由当前项必须进入竞争族 |
| 一个方向属性是否由数据确定？ | T05–T06、T13–T14 | 所有合法相容解释共享该属性；保留切换和混合标签 |
| 什么互反实验支持方向判断？ | T08–T09、T11 | 共同映射、已知任务结构及相应模型范围 |
| 有多少观测秩才足够？ | T10、T12 | 区分完整算子、指定功能量与受限性质；检查真实投影可见 |
| 定性方向后还可量化什么？ | T15、T17–T21 | 尖锐参数区间，或三个传播形状、双平面、有限核约束 |
| 单平面能否证明单一 SR 核？ | T06、T16 | 不能；额外模板只检查共同作用，核识别需额外任务信息 |
| 时间持续性会不会制造后瞻模式？ | T23 | 记录历史替代机制与实际探测依赖需独立检查 |
| 有误差时何时可以声明？ | T24–T28 | 有效误差界或合法抽样覆盖；连续族完整排除及数值可靠性 |

最终推断区分四种输出：非空相容集合具有唯一预定属性；多个属性仍相容；所有规定模型均不相容；数据或推断条件不足。最后一类包括重复单位不足、协方差不可用、相关信号不可见或覆盖假设无依据。在 T06 的单平面自由均值并集中，均值饱和性通常使“所有模型均不相容”无法出现；额外任务块和共享参数限制才提供真正可失败的核适配检查。

实际独立重复应按完整配对图谱的真实抽样单位定义，而不是按重组 occurrence 数、读出模板数或施加的计算变换数定义。不重复事件编号不等于没有共享时间支持；共同预处理、慢漂移和长时程相关均需审计。跨参与者也不能在个体参考与机制不同的情况下直接当作同一机制的 iid 重复。

## A11. 代数复现与结果索引

`code/math_checks.py` 与 `code/math_additions.py` 复现精确构造、补偿公式、逆解和扰动实例。计算输出与部分折扣界位于 `source_data/`；登记表 `math/theorem_registry.json` 索引 T01–T28。一般结论由以上证明给出，真实测量方法及来源见 B 部分。
