### Tasks, time reversal and the common observation model

A task is defined by a finite row-stochastic matrix $P$ and fixed state labels, with $P_{ij}$ the one-step probability of moving from $i$ to $j$. Let $\pi$ be a strictly positive stationary distribution satisfying $\pi^\top P=\pi^\top$, and set $D=\operatorname{diag}(\pi)$. Backward conditional probabilities are $P^*=D^{-1}P^\top D$. Complete SR and its pure-future part are $M_\gamma(P)=(I-\gamma P)^{-1}$ and $N_\gamma(P)=P(I-\gamma P)^{-1}$, with $0\leq\gamma<1$. They satisfy $M_\gamma=I+\gamma N_\gamma$; at zero discount, complete SR equals $I$ and the pure-future part equals $P$.

Population patterns in condition $c$ are $Y_c=g_cf_c(P)B+\mathbf1b_c^\top$. Here $B$ is a shared linear mapping, $g_c>0$ is overall gain and $b_c$ is a state-independent offset. Set $\widetilde P=D^{1/2}PD^{-1/2}$ and let $U$ have orthonormal columns spanning the orthogonal complement of $\sqrt\pi$. Define

$$
A=U^\top\widetilde P U,\qquad
X_c=U^\top D^{1/2}Y_c,\qquad
R=U^\top D^{1/2}B.
$$

Then $X_c=g_cf_c(A)R$, with the state-independent component removed. The real matrix $A$ is a contraction, and time reversal is represented by $A^\top$ in these coordinates. Common-mapping and visibility results are defined at this level (Supplementary Proposition T01). Experimental and behavioural constraints can specify the allowed internal transitions, which are then retained jointly with the neural constraints.

### Observation levels, candidate families and target properties

The parameter space contains mechanism directions, discounts, overall gains, current and predictive weights, and the references and internal transitions allowed by each analysis. Different direction labels form a labelled disjoint union: two labels remain distinct explanations even when they yield the same operator. For a specified task and observation, $\Phi(\theta)$ is the candidate prediction and $\mathcal C(y)=\{\theta:\Phi(\theta)=y\}$ is the exact compatibility set. Target properties include pure direction combinations, net direction of mixtures, discount ranges and kernel relations. Global property identification requires $\Phi(\theta_1)=\Phi(\theta_2)$ to imply agreement of the property (Supplementary Proposition T05).

Geometric results use squared Euclidean row distances for every state pair in canonical state features. The same-discount comparison for doubly stochastic tasks is treated in T02. Its nonuniform extension compares row distances of $M_\gamma(P)D^{-1/2}$ and $M_\gamma(P^*)D^{-1/2}$, specifying the feature metric $D^{-1}$. General neural observations $f(P)B$ retain their actual mapping. Direct references, learned same-task references and learned reciprocal references use their corresponding observation relations.

The directly referenced generalized direction family has $c\in\mathbb R$, $\alpha\geq0$ and $0\leq\gamma<1$; point-parameter identification uses $\alpha>0$. Reciprocal four-family angular results use pure complete SR, positive overall gain and $0\leq\gamma,\eta\leq g<1$. The full discount domain is obtained as $g\uparrow1$. Nonnegative mixture extensions require a common feature mapping across scales and a nonnegative additional current term. Directional properties and component purity are specified separately.

### Deterministic task constructions and figure calculations

In the seven-state Paley task, each state moves with probability $1/3$ to relative indices $1,2,4$. After removing the constant direction, rational arithmetic verifies that forward target discount $3/10$ and reference discount $4/5$ produce a relative operator equal to complete backward SR with discount $75/173$ and gain $634/519$. The identity holds on the complete six-dimensional state-contrast space (T04).

The eight-state cycle has $P_{ij}=1$ when $j=i+1\pmod8$. The prespecified plane

$$
C_{1j}=\tfrac12\cos(\pi j/2),\qquad
C_{2j}=-\tfrac12\sin(\pi j/2)
$$

satisfies $CC^\top=I$, $CP=JC$ and $CP^3=-JC$, where $J=\begin{psmallmatrix}0&-1\\1&0\end{psmallmatrix}$. The two-plane example uses the $90^\circ$ and $45^\circ$ Fourier blocks of the same cycle, with shared current weight, predictive gain and discount. These planes are specified by the task structure.

The finite-memory task contains all binary words of length $k$. It moves from $x_1\cdots x_k$ to $x_2\cdots x_kb$, with $b\in\{0,1\}$ equally likely. Contrasts $C_{j,x}=2^{-k/2}(1-2x_j)$ realize the lower-shift block $N$, with $CP=NC$, $CP^\top=N^\top C$ and $N^k=0$. Main-text examples use $k=4$. Paths, sectors, compensation curves, kernel matrices and confidence geometry in the figures follow analytic relations and deterministic calculations; the parameters and source data accompany each figure.

### Task-compatible readout and functional visibility

For a channel-pattern matrix $X$ on the two-dimensional task block, let $R_0$ be a template determined from independent training data. Define

$$
z_{R_0}(X)=
\begin{pmatrix}
\langle R_0,X\rangle_F\\
\langle JR_0,X\rangle_F
\end{pmatrix}.
$$

It satisfies $z_{R_0}(JX)=Jz_{R_0}(X)$. Common positive normalization and a frozen right-side channel transformation preserve the relation. If $X_1=(uI+vJ)X_0$, the readouts obey the same relation and $\det[z(X_0),z(X_1)]=v\|z(X_0)\|^2$. Directional information requires a nonzero population reference projection; a nonzero training template alone does not specify population visibility (T12).

More generally, for $Y=TX_0$, set $\Pi=X_0X_0^+$. With unrestricted $T$, the functional $\langle T,W\rangle_F$ is observable if and only if $W=W\Pi$ (T10). A full operator requires a full-row-rank reference; a particular directional or kernel quantity requires visibility of the relevant action. Multiple templates are organized within the same actual repeated unit, preserving their cross-covariances.

### Recovery of parameters and finite temporal kernels

The three-propagation-shape condition follows by comparing two generalized SR observations, clearing the inverse matrices and examining the resulting linear combination of $B,PB,P^2B$. For the two rotation blocks, let $q=(v_1/\beta_1)/(v_2/\beta_2)$ and

$$
t=\frac{q-1}{2(qa_1-a_2)},\qquad
\gamma=\frac{2t}{1+\sqrt{1-4t^2}}.
$$

Within the domain of T18, this gives the unique discount. Each block then recovers the common $\alpha,c$. Cross-block consistency is retained as a model check. Parameter ranges under noise are obtained by projecting the compatibility set.

The directly referenced finite-kernel relation is $T=a\sum_{h=0}^{k-1}w_hN^h+cI$. Define

$$
\psi_h=\frac{1}{k-h}\sum_{i=1}^{k-h}
\bigl(T_{i+h,i}-T_{i,i+h}\bigr),\qquad 1\leq h<k.
$$

A forward kernel gives $\psi_h=aw_h$ and a backward kernel gives $-aw_h$. Single-geometric ratios and mixture variance use these quantities with their common scale removed. For a learned reverse reference, recovery uses the first-column relation of T21. For pure complete SRs with $k\geq2$, $a=T_{11}$, $\gamma=T_{21}/T_{11}$ and $\eta=-T_{12}/T_{11}$. The two reference settings have different recovery formulas; remaining matrix relations test the specified kernel model.

### Joint uncertainty and continuous-candidate inversion

Conditional on independent training data, state labels, task blocks, regions of interest, windows, channel processing, templates and the full candidate family are fixed. Formal repeated units are $Z_i=(B_i^\top,Y_i^\top)^\top$. Candidate $T_\theta$ defines a $q$-dimensional residual $E_i(\theta)=Y_i-T_\theta B_i$. Assume independent and identically distributed Gaussian pairs, a true candidate satisfying the population mean relation, positive-definite residual covariance and $n>q$. Following Hotelling’s multivariate mean test [26], use unbiased sample covariance $S_{E,\theta}$ to calculate

$$
R(\theta)=\bar E_\theta^\top S_{E,\theta}^{-1}\bar E_\theta,
\qquad
K_q(\alpha)=\frac{q(n-1)}{n(n-q)}F_{q,n-q}(1-\alpha).
$$

The retained set is $\widehat\Theta_\alpha=\{\theta:R(\theta)\leq K_q(\alpha)\}$. Under these conditions, the true candidate is retained with probability $1-\alpha$. A property is declared when the set is nonempty and all candidates agree; the probability of making an incorrect declaration is at most $\alpha$ (T26). The guarantee follows from retention of the true candidate and applies to continuous inversion. Selection of windows, templates and tasks is handled through frozen training or the corresponding multiplicity control.

For a single template, write reference–target sample covariance blocks as $S_B,S_Y,S_{BY},S_{YB}$. The residual covariance is

$$
V_T=S_Y+TS_BT^\top-TS_{BY}-S_{YB}T^\top.
$$

Cross terms retain within-pair dependence. A two-dimensional residual uses $K_2$; a four-dimensional joint-mean region uses $K_4$ and $n>4$. Computing a two-dimensional candidate residual through four-dimensional auxiliary geometry additionally requires positive-definite $S_Z$, generally with $n>4$; the residual test itself needs only $n>2$ and positive-definite residual covariance. Additional population-mean properties use an actual joint-coverage region (T27). Cross-covariance bias in directional area and worst-case margins for externally bounded reference drift are given in T28.

### Full operators and deterministic error budgets

For complete-operator recovery, set $\widehat T=\widehat Y\widehat B^+$. If the population relation is $Y=TB$, $\widehat B$ has full row rank, $\|\widehat B-B\|_F\leq e_B$, $\|\widehat Y-Y\|_F\leq e_Y$ and $\sigma_{\min}(\widehat B)>e_B$, then

$$
\|\widehat T-T\|_F\leq
\frac{e_Y+\|\widehat T\|_2e_B}
{\sigma_{\min}(\widehat B)-e_B}.
$$

Errors for linear functionals follow by their norms (T24). For a nearest-family decision, distance $d$ from the true prediction to the incorrect family and actual observation error $e$ give the sufficient condition $d>2\|e\|$. If the correct class instead uses a frozen training prediction with transfer error $h$ relative to the truth, the condition becomes $d>2\|e\|+h$ (T25). Each calibration matches its error space, decision rule and reference relation.
