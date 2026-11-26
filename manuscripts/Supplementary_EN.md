---
title: "Supplementary Information: Identifying direction and temporal weighting in predictive cognitive maps"
lang: en
author: "Yongzhang Ji"
affiliation: "Faculty of Psychology, Southwest University, Chongqing, China"
orcid: "0009-0006-4111-2000"
date: "2026-11-26"
---

This file provides the general proofs supporting the main text, complete methods and results for real repeated measurements and known-transformation calibration, and related-work context. Results T01–T28 are cross-referenced in the article; empirical sections E1–E11 collect all analysis branches and their sources. Editable vector figures, source data and reproduction code accompany the paper.

# A. Theory and statistical proofs

This supplement proceeds from the observation model to non-identifiability, experimental remedies, observable readouts and uncertainty. Every directional claim is relative to prespecified state labels and a prespecified transition direction. Results T01–T28 have the same identifiers in the manuscript, theorem registry and Chinese supplement. We distinguish recovery of an unrestricted operator, identification of a property shared by a restricted compatible family, and estimation of parameters within that family. These objectives require different information.

## A1. Models, state contrasts and observation levels

### A1.1 Notation

Let $P\in\mathbb R^{n\times n}$ be row stochastic, with $P_{ij}$ the probability of transitioning from state $i$ to state $j$. Assume a strictly positive stationary distribution, written as a column vector $\pi$, with $\pi^\top P=\pi^\top$ and $\mathbf1^\top\pi=1$. Write $D=\operatorname{diag}(\pi)$. Stationary time reversal is

$$
P^*=D^{-1}P^\top D.
$$

This follows from $\Pr(S_{t-1}=j\mid S_t=i)=\pi_jP_{ji}/\pi_i$. Transposition alone represents time reversal only when the stationary distribution is uniform. Aperiodicity is unnecessary.

The complete successor representation, including the current state, and its future-only counterpart are

$$
M_\gamma(P)=(I-\gamma P)^{-1}
=\sum_{h\ge0}\gamma^hP^h,\qquad
N_\gamma(P)=P(I-\gamma P)^{-1},\qquad0\le\gamma<1.
$$

They satisfy $M_\gamma=I+\gamma N_\gamma$. Thus $M_0=I$ contains no directional information, whereas $N_0=P$ generally contains one-step directional information. We retain this distinction throughout.

Let $\widetilde P=D^{1/2}PD^{-1/2}$ and $r=\sqrt\pi$, and choose $U\in\mathbb R^{n\times(n-1)}$ with orthonormal columns and $U^\top r=0$. Define

$$
A=U^\top\widetilde P U=S+K,\qquad
S=(A+A^\top)/2,\qquad K=(A-A^\top)/2.
$$

$S$ is the symmetric component and $K$ the antisymmetric component relative to task direction. Write $F_\gamma(A)=(I-\gamma A)^{-1}$, $\operatorname{skew}T=(T-T^\top)/2$ and $\langle V,W\rangle_F=\operatorname{tr}(V^\top W)$. For two-dimensional rotation blocks we consistently use

$$
J=\begin{pmatrix}0&-1\\1&0\end{pmatrix},
\qquad J^2=-I,\qquad J^\top=-J.
$$

### Proposition T01: stationary contrasts and a common observation relation

We have $\|\widetilde P\|_2\le1$, $\|A\|_2\le1$ and $-I\preceq S\preceq I$. Suppose the population state patterns in condition $c$ satisfy

$$
Y_c=g_cM_c(P)B+\mathbf1b_c^\top,\qquad g_c>0,
$$

where $B$ is an arbitrary linear feature mapping shared across conditions and $b_c$ is a state-independent offset. Their stationary-weighted contrasts satisfy

$$
X_c=U^\top D^{1/2}Y_c=g_cM_c(A)R,
\qquad R=U^\top D^{1/2}B.
$$

**Proof.** Jensen's inequality gives

$$
\sum_i\pi_i\Big(\sum_jP_{ij}x_j\Big)^2
\le\sum_{ij}\pi_iP_{ij}x_j^2
=\sum_j\pi_jx_j^2.
$$

Hence $\widetilde P$ is a Euclidean contraction. Moreover, $\widetilde P r=r$ and $\widetilde P^\top r=r$, so $r^\perp$ is a common reducing subspace. The restricted operator is $A$, and its reversal is $A^\top$. Restricting a rational function to this subspace is equivalent to applying the function after restriction. The offset vanishes because $U^\top D^{1/2}\mathbf1=0$. The spectral bounds on $S$ follow from $\|A\|_2\le1$. ∎

**Scope.** Fixed recording channels or successful decoding do not by themselves establish a biologically common mapping $B$. If arbitrary $B_c$ are allowed, any invertible candidate $M_c$ can explain the data by choosing $B_c=g_c^{-1}M_c^{-1}(Y_c-\mathbf1b_c^\top)$. Nonlinear class probabilities do not automatically satisfy this linear observation relation. The mathematical $P$ is the transition operator used by the cognitive candidate; equating it with the task generator must be stated or independently constrained.

## A2. Exact equivalence caused by geometry and reference measurements

### Theorem T02: normality characterizes canonical geometric equivalence

If $P$ is doubly stochastic and $0<\gamma<1$, the complete squared Euclidean row-distance matrices of $M_\gamma(P)$ and $M_\gamma(P^\top)$ are equal if and only if $P$ is normal, $PP^\top=P^\top P$. Let

$$
G_f=M_\gamma(P)M_\gamma(P)^\top,\qquad
G_b=M_\gamma(P)^\top M_\gamma(P).
$$

Then

$$
G_f-G_b=\gamma^2G_f(PP^\top-P^\top P)G_b
$$

and

$$
\frac{\gamma^2}{(1+\gamma)^4}\|PP^\top-P^\top P\|_F
\le\|G_f-G_b\|_F
\le\frac{\gamma^2}{(1-\gamma)^4}\|PP^\top-P^\top P\|_F.
$$

**Proof.** The difference between the inverse Gram matrices is $\gamma^2(PP^\top-P^\top P)$. Applying $X-Y=X(Y^{-1}-X^{-1})Y$ gives the first identity. Both Gram matrices map $\mathbf1$ to $(1-\gamma)^{-2}\mathbf1$, so their difference is already centered. Squared distances determine the centered Gram matrix; complete distance equality is therefore equivalent to $G_f-G_b=0$. Invertibility and $\gamma>0$ then give normality. Finally, $\sigma_{\min}(I-\gamma P)\ge1-\gamma$ and $\|I-\gamma P\|_2\le1+\gamma$ provide the singular-value bounds for the product. ∎

**Nonuniform extension.** Use row patterns $M_\gamma(P)D^{-1/2}$ and $M_\gamma(P^*)D^{-1/2}$, explicitly choosing the feature metric $D^{-1}$. The corresponding condition is normality of $\widetilde P$. To prove this, let

$$
\Delta=F_\gamma(\widetilde P)F_\gamma(\widetilde P)^\top
-F_\gamma(\widetilde P)^\top F_\gamma(\widetilde P),
\qquad G=D^{-1/2}\Delta D^{-1/2}.
$$

We have $G\pi=0$. If all distances agree, symmetry implies $G=\mathbf1v^\top+v\mathbf1^\top$ for some $v$. Thus $v+(v^\top\pi)\mathbf1=0$; premultiplying by $\pi^\top$ gives $v=0$. Distance equality is therefore equivalent to $\Delta=0$, after which the preceding argument applies. This does not establish the same result for unweighted distances with nonuniform stationarity.

**Future-only extension.** For $N_\gamma(P)$, the same canonical result includes $\gamma=0$. Equal row and column sums again reduce distance equality to normality of $N_\gamma$. Rational functions of a normal $P$ are normal; conversely,

$$
P=N_\gamma(I+\gamma N_\gamma)^{-1}
$$

and the spectral theorem give normality of $P$, without requiring $P$ itself to be invertible. In contrast, equality of the complete maps $M_0=I$ does not imply normality.

These statements concern canonical feature geometry, the same discount, and all distances. They do not automatically cover arbitrary anisotropic readouts $M_\gamma B$, correlation distances, a subset of state pairs or fits with different discounts on the two sides.

### Corollary T03: strongly directed and exponentially large labelled isometric families

For a cyclic permutation $R$ on at least three states, let

$$
P=pR+qR^\top+(1-p-q)I,\qquad
p,q>0,\quad p+q\le1.
$$

This matrix is normal and doubly stochastic. When $p\ne q$ the process is irreversible, with stationary entropy production $(p-q)\log(p/q)$, which diverges as $q\downarrow0$ at fixed $p$. Nevertheless, forward and backward canonical geometries agree at every discount.

More generally, choose orthonormal pairs $(u_\ell,v_\ell)$ in $\mathbf1^\perp$ and define

$$
K_\sigma=\sum_{\ell=1}^k\sigma_\ell\epsilon_\ell
(u_\ell v_\ell^\top-v_\ell u_\ell^\top),\qquad
P_\sigma=\mathcal J_n+s(I-\mathcal J_n)+K_\sigma,
\quad \mathcal J_n=\mathbf1\mathbf1^\top/n,
$$

where $0\le s<1$, $\sigma_\ell\in\{-1,1\}$ and the positive $\epsilon_\ell$ are sufficiently small to keep all transition entries nonnegative. These $2^k$ labelled transition matrices are isometric at every discount.

**Proof.** All terms in the cycle construction commute, establishing normality; pairing opposite directed edges gives the entropy-production formula. In the second construction, orthogonality of the planes gives

$$
K_\sigma^2=-\sum_\ell\epsilon_\ell^2
(u_\ell u_\ell^\top+v_\ell v_\ell^\top),
$$

independent of the signs. The contrast-space Gram matrix is

$$
\big[(1-\gamma s)^2I-\gamma^2K_\sigma^2\big]^{-1}.
$$

The constant direction is also identical across the family. ∎

This counts labelled matrices, not non-isomorphic unlabelled graphs, and does not establish a uniform positive separation as dimension grows.

### Proposition T04: direct-reference positivity and complete aliases from a learned reference

If $X_0=g_0R$ has full row rank and $X_1=g_1F_\gamma(A)R$, then $T=X_1X_0^+=(g_1/g_0)F_\gamma(A)$. A learned complete SR for an independent-transition reference $P_0=\mathbf1q^\top$ also reduces to $g_0R$ after state contrasts. For $K\ne0$, define

$$
\omega(T)=\frac{\langle T,K\rangle_F}{\|K\|_F^2}.
$$

For a directly referenced complete SR, $\gamma>0$ yields a positive forward score and a negative backward score. For two learned forward maps of the same graph,

$$
T=aF_\gamma(A)(I-\eta A)
=a[I+(\gamma-\eta)AF_\gamma(A)]
$$

and

$$
\omega(T)=a(\gamma-\eta)\rho_\gamma/\|K\|_F^2,
\qquad \rho_\gamma>0.
$$

**Proof.** Full row rank implies $X_0X_0^+=I$. The nonidentity component of the independent-transition reference has identical rows and disappears under contrasts. To establish positivity, set

$$
H=I-\gamma S\succ0,\quad
L=H^{-1/2}KH^{-1/2},\quad R_\gamma=(I-\gamma L)^{-1}.
$$

Then $F_\gamma(A)=H^{-1/2}R_\gamma H^{-1/2}$, and the difference-of-inverses identity gives

$$
\operatorname{skew}F_\gamma=\gamma F_\gamma K F_\gamma^\top,
\qquad
\operatorname{skew}(AF_\gamma)=F_\gamma K F_\gamma^\top.
$$

The second identity also holds at $\gamma=0$. Cyclicity of the trace and $R_\gamma^\top L R_\gamma=L(I-\gamma^2L^2)^{-1}$ give

$$
\rho_\gamma
=\langle F_\gamma K F_\gamma^\top,K\rangle_F
=-\operatorname{tr}\{L^2(I-\gamma^2L^2)^{-1}\}>0.
$$

Here $-L^2$ is nonzero positive semidefinite; no commutation of $H$ and $L$ is assumed. Taking inner products establishes the stated scores. ∎

**A complete contrast-operator counterexample.** Define the seven-state Paley transition by $P_{ij}=1/3$ when $j-i\pmod7\in\{1,2,4\}$, and zero otherwise. Its six-dimensional contrast block satisfies

$$
A=-I/6+K,\qquad K^2=-7I/36.
$$

Using $F_t(A)=\{(1+t/6)I+tK\}/\{(1+t/6)^2+7t^2/36\}$ gives the exact identity

$$
F_{3/10}(A)F_{4/5}(A)^{-1}
=\frac{53}{48}I-\frac{25}{56}K
=\frac{634}{519}F_{75/173}(A^\top).
$$

Thus the relative operator between two genuinely forward maps with different horizons is exactly a lawful backward complete SR multiplied by a positive gain. Equality holds over the entire six-dimensional contrast operator, not merely its Euclidean geometry. The constant state direction is excluded. Equality of deterministic means does not by itself imply equality of arbitrary noise distributions.

## A3. How continuous families change the identifiable claim

### Proposition T05: properties on observation fibres and transition anchoring

Let $\Phi:\Theta\to\mathcal Y$ be an observation map including all lawful mechanisms and nuisance parameters. A property $q(\theta)$ is determined by a noiseless observation if and only if it is constant on each nonempty fibre $\Phi^{-1}(y)$.

**Proof.** If two parameters have the same observation but different properties, no rule using that observation alone can be correct for both. Conversely, if each fibre has a constant property, define $\widetilde q(y)=q(\theta)$ for any $\theta\in\Phi^{-1}(y)$; this is well defined. ∎

This is a standard identification principle. The mathematical work here is to calculate concrete SR observation families and their fibres, not to rename that principle.

If the admissible transition set $\mathcal P$ is closed under stationary reversal, then

$$
\{M_\gamma(P):P\in\mathcal P\}
=\{M_\gamma(P^*):P\in\mathcal P\}.
$$

Because $P^{**}=P$, the two sets are merely re-indexed. Directional claims therefore require learned dynamics to be constrained by design or independent evidence. If a behavioural parameter set has noncoverage probability at most $\beta$ and a neural set at most $\alpha$, both are valid with probability at least $1-\alpha-\beta$. This union bound does not require independence.

### Theorem T06: collapse of the directly referenced single-plane family to half-planes

Let

$$
A_\sigma=aI+\sigma\beta J,\quad
a^2+\beta^2=1,\quad\beta>0,\quad\sigma\in\{-1,1\},
$$

be an orthogonal rotation block, and assume a qualified direct reference. Consider

$$
G_\sigma=cI+\alpha N_\gamma(A_\sigma),
\quad c\in\mathbb R,\quad\alpha\ge0,\quad0\le\gamma<1.
$$

Then

$$
G_\sigma=uI+vJ,\qquad D_\gamma=1-2a\gamma+\gamma^2>0,
$$

$$
u=c+\alpha\frac{a-\gamma}{D_\gamma},
\qquad v=\sigma\alpha\frac{\beta}{D_\gamma}.
$$

The image of the entire family is exactly the closed half-plane $\sigma v\ge0$; if $\alpha>0$ is required, it is the open half-plane $\sigma v>0$. For any nonzero reference $b$, observing $y=G_\sigma b$ identifies $u,v$ and the sign of net rotation, but does not identify $\gamma$.

**Proof.** Multiplication using $(I-\gamma A_\sigma)^{-1}=\{(1-a\gamma)I+\sigma\beta\gamma J\}/D_\gamma$ gives the coefficients. For any admissible $(u,v)$ and any $\gamma$, choose

$$
\alpha=|v|D_\gamma/\beta,\qquad
c=u-|v|(a-\gamma)/\beta.
$$

This produces the same operator, proving both the exact image and the discount alias. Since $b,Jb$ are orthogonal with equal length,

$$
u=\frac{b^\top y}{\|b\|^2},\qquad
v=\frac{(Jb)^\top y}{\|b\|^2}
=\frac{\det[b,y]}{\|b\|^2}.
$$

The union of the two half-planes can therefore fit any target mean exactly for one nonzero reference. ∎

**Extensions and limits.** This single-plane mean relation leaves no residual restriction with which to test a single geometric SR kernel. If nonnegative forward and backward future-only components are mixed, the net coefficient becomes

$$
v=\beta\left(\frac{\alpha_F}{D_{\gamma_F}}
-\frac{\alpha_B}{D_{\gamma_B}}\right).
$$

Both weights may be strictly positive while producing any specified $v$; free $c$ then matches $u$. A positive sign therefore does not establish pure forward composition. Restricting $c\ge0$, gains or cross-block parameter sharing changes the family and requires recomputing the compatible set.

### Counterexample T07: reversing both learned mechanisms without changing the observed patterns

In task blocks $J,-J$, set $\gamma=1/2$ and $\alpha=5$. A forward mechanism with current coefficient $c_F=1$ gives

$$
G^F_+=-I+4J,\qquad G^F_-=-I-4J.
$$

A backward mechanism with $c_B=3$ gives

$$
G^B_+=I-4J=-G^F_+,\qquad
G^B_-=I+4J=-G^F_-.
$$

Changing the unknown common baseline mapping from $B$ to $B'=-B$ yields $G^F_\pm B=G^B_\pm B'$. Both complete condition patterns and the relative operator agree:

$$
T_{\rm rel}=(-15I-8J)/17.
$$

**Proof.** Use $N_{1/2}(J)=(-\tfrac12 I+J)/(1+\tfrac14)$ and multiply directly. All current coefficients and future gains in the construction are positive. ∎

This counterexample permits an unknown common mapping that differs between candidate explanations. An independently qualified direct reference fixing the same $B$ rules out the change to $-B$, so the direct-direction result T06 remains intact. Nor does the example refute the theorem for a complete SR plus a nonnegative additional current term: that theorem concerns a different family that ties current and future weights.

## A4. Resolving reference ambiguity with reciprocal tasks

### Theorem T08: reciprocal signs for any real contraction block

Let $A=S+K$ be a real contraction with $K\ne0$. The two tasks are $A,A^\top$ in the same contrast coordinates. If both conditions use a complete forward SR, allowing separate discounts and positive gains, their exact relative operator is

$$
T=aF_\gamma(A)(I-\eta A^\top),
\qquad a>0,\quad0\le\gamma,\eta<1,\quad\gamma+\eta>0.
$$

Then $\langle T,K\rangle_F>0$. If both conditions are backward, the sign is negative. A symmetric mechanism constructed from the same symmetric part $S$ gives zero.

**Proof.** Write $F=F_\gamma(A)$. Inserting $F^{-1}$ and $F^{-\top}$ gives

$$
\operatorname{skew}T/a
=F\{(\gamma+\eta)K-\gamma\eta(SK+KS)\}F^\top.
$$

Use $H=I-\gamma S$ and $L=H^{-1/2}KH^{-1/2}$ as in T04, and define

$$
V=(I-\gamma^2L^2)^{-1},\qquad Z=-L^2V\succeq0.
$$

$K\ne0$ implies $Z\ne0$. Since $S$ commutes with $H$, the expression in braces equals $H^{1/2}L_0H^{1/2}$, where

$$
L_0=(\gamma+\eta)L-\gamma\eta(SL+LS).
$$

Cyclicity of the trace yields

$$
\begin{aligned}
\langle T,K\rangle_F/a
&=-\operatorname{tr}(R_\gamma L_0R_\gamma^\top L)\\
&=-\operatorname{tr}(L_0LV)\\
&=\operatorname{tr}\{[(\gamma+\eta)I-2\gamma\eta S]Z\}.
\end{aligned}
$$

The last equality uses commutation of $L$ and $V$, not of $S$ and $L$. Moreover,

$$
(\gamma+\eta)I-2\gamma\eta S
\succeq[\gamma(1-\eta)+\eta(1-\gamma)]I\succ0.
$$

The trace of a positive-definite matrix times a nonzero positive-semidefinite matrix is strictly positive. Replacing $A$ with $A^\top=S-K$ changes the directional marker's sign. In the symmetric mechanism, both factors are functions of one symmetric matrix and their product is symmetric. ∎

**Zero boundaries.** A single zero discount still permits strict directionality; when both vanish, $T=aI$. If $K=0$, the task provides no directional marker. Allowing discounts, irreversibility or gains to approach zero destroys any implication of a uniform positive noise margin. This theorem does not include condition-switching FB or BF mechanisms; those enter explicitly in T13.

### Proposition T09: established extensions with nonnegative additional current components

Let $r_0,s_0\ge0$ and consider

$$
T=a[F_\gamma(A)+r_0I][F_\eta(A^\top)+s_0I]^{-1}.
$$

The sign in T08 is retained for a two-dimensional real contraction block, a normal block of arbitrary dimension, and any real block satisfying $A^2=0$.

**Two-dimensional proof.** Orthogonally diagonalizing the symmetric part gives

$$
A=\begin{pmatrix}s_1&k\\-k&s_2\end{pmatrix},\quad
h_i=1-\gamma s_i>0,\quad j_i=1-\eta s_i>0.
$$

Set $d_\gamma=h_1h_2+\gamma^2k^2$, $d_\eta=j_1j_2+\eta^2k^2$, and

$$
\Delta=(j_2+s_0d_\eta)(j_1+s_0d_\eta)+\eta^2k^2>0.
$$

Explicit inversion yields $\operatorname{skew}T=cK$, with

$$
c=\frac{a d_\eta}{2d_\gamma\Delta}
\{\eta(h_1+h_2+2r_0d_\gamma)
+\gamma(j_1+j_2+2s_0d_\eta)\}>0.
$$

All terms are nonnegative and $\gamma+\eta>0$, proving the result. A real normal matrix decomposes orthogonally into real one-dimensional and two-dimensional rotation blocks; summing their contributions establishes that case. If $A^2=0$, direct expansion gives

$$
\operatorname{skew}T
=a\left\{\frac{\gamma}{1+s_0}
+\frac{(1+r_0)\eta}{(1+s_0)^2}\right\}K.
$$

The remaining terms are multiples of $I$ or $AA^\top$ and hence symmetric. ∎

No general claim is made here for arbitrary high-dimensional nonnormal operators, independent extra feature mappings, or freely varying current-to-future ratios.

## A5. Rank requirements and task-compatible measurement

### Theorem T10: the exact observability condition for an unrestricted operator functional

Suppose $X_1=TX_0$ without restricting $T$ to a parameterized family, and let $\Pi=X_0X_0^+$. The functional $\langle T,W\rangle_F$ is determined by the data if and only if

$$
W=W\Pi.
$$

**Proof.** The data determine $T\Pi=X_1X_0^+$. If $W=W\Pi$, then $\langle T\Pi,W\rangle_F=\langle T,W\rangle_F$. Otherwise, take $\Delta=W(I-\Pi)\ne0$. We have $\Delta X_0=0$ but

$$
\langle\Delta,W\rangle_F=\|W(I-\Pi)\|_F^2>0.
$$

Thus $T$ and $T+\Delta$ produce the same observations and different functional values. ∎

An unrestricted full operator requires $X_0$ to have full row rank. With $n$ states and the constant direction removed, having $p\ge n-1$ linear channels is only a necessary counting condition; it does not ensure independence or stability of the patterns. A specified marker requires only its relevant directions to be observed. A restricted parameter family can require still less, as in T12 and T17. The necessary condition for an unrestricted operator must not be treated as necessary for every mechanistic property.

### Proposition T11: common reducing blocks and the eight-cycle $P/P^3$ intervention

If $E^\top E=I$ and

$$
AE=EA_E,\qquad A^\top E=EA_E^\top,
$$

then $E^\top f(A)=f(A_E)E^\top$ for every rational function defined on the relevant spectra. SR analysis can therefore be performed after selecting a common task block without admitting unknown contributions from the other blocks.

**Proof.** The two relations imply that both the subspace and its orthogonal complement are invariant. In the corresponding orthonormal basis, $A$ is block diagonal. Powers, inverses and rational functions preserve this structure. ∎ An arbitrary state projection need not have this property.

For the eight-state row cycle $P_{j,j+1\bmod8}=1$ and $Q=P^3$, define

$$
C_{1j}=\tfrac12\cos(\pi j/2),\qquad
C_{2j}=-\tfrac12\sin(\pi j/2),\qquad j=0,\ldots,7.
$$

Direct multiplication gives

$$
C\mathbf1=0,\quad CC^\top=I,\quad
CP=JC,\quad CQ=-JC,\quad
CP^\top=-JC,\quad CQ^\top=JC.
$$

Although the complete graphs are not reversals, $Q\ne P^\top$, this prespecified plane is exactly reciprocal. For complete SRs,

$$
T_E=a(I-\gamma J)^{-1}(I+\eta J)
=\frac{a}{1+\gamma^2}
\{(1-\gamma\eta)I+(\gamma+\eta)J\}.
$$

The full marker

$$
W=\tfrac12(P-P^\top-Q+Q^\top)=C^\top(2J)C
$$

has rank two and $\|W\|_F^2=8$, with

$$
\frac{\langle a(I-\gamma P)^{-1}(I-\eta Q),W\rangle_F}{8}
=\frac{a(\gamma+\eta)}{2(1+\gamma^2)}.
$$

**Verification.** Expand $(I-\gamma P)^{-1}=(\sum_{j=0}^7\gamma^jP^j)/(1-\gamma^8)$. Distinct cyclic powers are Frobenius-orthogonal, and the marker extracts only the coefficients of powers 1, 7, 3 and 5. Combining them gives numerator $(\gamma+\eta)(1-\gamma^2+\gamma^4-\gamma^6)$; cancellation of $(1-\gamma^2)(1+\gamma^4)$ yields the formula. The relation $W=C^\top(2J)C$ can also be checked entrywise. Reversing both mechanisms transposes the relevant operator and reverses the score. The two commuting symmetrized transitions yield a zero score. ∎

This is a statement about task structure, not a claim that any recorded participant already satisfies the common-mapping, visibility or noise assumptions.

### Theorem T12: task-compatible readout and weak visibility

For a nonzero template $R\in\mathbb R^{2\times p}$ frozen after training, define

$$
z_R(X)=\frac1{\|R\|_F}
\begin{pmatrix}\langle R,X\rangle_F\\
\langle JR,X\rangle_F\end{pmatrix}.
$$

For all $X,u,v$,

$$
z_R(JX)=Jz_R(X),\qquad
z_R((uI+vJ)X)=(uI+vJ)z_R(X).
$$

If population patterns satisfy $X_+=TX_-$ with $T=uI+vJ$, the compressed means $x_\pm=z_R(X_\pm)$ retain that relation, and

$$
\det[x_-,x_+]=v\|x_-\|^2.
$$

Thus a nonzero true reference projection, $x_-\ne0$, suffices to identify the sign of $v$; full row rank of the two-dimensional reference is unnecessary.

**Proof.** Since $J^\top=-J$, $\langle R,JX\rangle_F=-\langle JR,X\rangle_F$ and $\langle JR,JX\rangle_F=\langle R,X\rangle_F$, proving the intertwining relation. Also, $R,JR$ are orthogonal and equally long, implying $\|z_R(X)\|\le\|X\|_F$. The determinant identity follows from $\det[x,Jx]=\|x\|^2$. ∎

This structural identity does not require the trained template to equal the true reference. It also does not guarantee that the evaluation mean projects nontrivially onto that template. Signal orthogonal to $\operatorname{span}\{R,JR\}$ is lost. A different linear compression $L$ must independently satisfy $LT_\theta=\widetilde T_\theta L$ to preserve the candidate action. Variance-based state PCA or a classifier is not automatically compatible.

**Connecting two learned conditions in a general rotation block.** Let $A=\rho I+\kappa J$, $\rho^2+\kappa^2\le1$, $\kappa>0$. The pure forward relative operator has

$$
u=\frac{a\{(1-\gamma\rho)(1-\eta\rho)-\gamma\eta\kappa^2\}}
{(1-\gamma\rho)^2+\gamma^2\kappa^2},
$$

$$
v=\frac{a\kappa(\gamma+\eta-2\gamma\eta\rho)}
{(1-\gamma\rho)^2+\gamma^2\kappa^2}>0.
$$

Positivity follows from $\gamma(1-\eta\rho)+\eta(1-\gamma\rho)>0$. In a general block, $u$ can be negative. The more restrictive positive-$u$ sectors below require the quarter-turn block. This connection allows a learned reference but does not remove the alias caused by free current components in T07.

## A6. Four mechanism families, mixtures and partial quantitative identification

### Theorem T13: exact angular images of the relative quarter-turn families

In task blocks $J,-J$, let both conditions use pure complete SRs, each with discount in $[0,g]$, $g<1$, and allow any positive gain ratio. F and B denote forward and backward mechanisms relative to each condition's own task. Put $h=\arctan g$ and write

$$
T=r(\cos\phi\,I+\sin\phi\,J),\qquad r>0.
$$

The exact angular images are

| Mechanism label | Complete range of relative angle $\phi$ |
|---|---|
| FF | $[0,2h]$ |
| BB | $[-2h,0]$ |
| FB | $[-h,h]$ |
| BF | $[-h,h]$ |

Every listed angle and every positive radial scale are attainable.

**Proof.** A single complete SR satisfies

$$
F_\gamma(J)=\frac{I+\gamma J}{1+\gamma^2}
=\frac1{\sqrt{1+\gamma^2}}
\{\cos\theta_\gamma I+\sin\theta_\gamma J\},
\qquad\theta_\gamma=\arctan\gamma.
$$

A forward mechanism on the reversed task has angle $-\theta_\eta$. The relative angle is target minus reference; hence the four angles are

$$
\theta_\gamma+\theta_\eta,\quad
-\theta_\gamma-\theta_\eta,\quad
\theta_\gamma-\theta_\eta,\quad
-\theta_\gamma+\theta_\eta.
$$

The independent component angles lie in $[0,h]$, giving the ranges. For FF, choose both angles equal to $\phi/2$. For the switching families, set one angle to zero and the other to the required magnitude. The free gain ratio matches any $r>0$. ∎

For FF specifically,

$$
\phi=\arctan\gamma+\arctan\eta,\qquad
r=a\sqrt{\frac{1+\eta^2}{1+\gamma^2}}.
$$

Thus $0<\phi\le h$ is compatible with FF and a switching mechanism: positive area alone does not identify FF. If a nonempty compatible set lies entirely at $\phi>h$, only FF remains. In observable coefficients, the condition is $u>0$ and $v-gu>0$. Zero discounts and zero net-direction boundaries must remain in the candidate set. Additional bounds on gain also restrict the radial image.

For the open discount domain $[0,1)$, the FF and BB ranges are $[0,\pi/2)$ and $(-\pi/2,0]$, whereas switching families cover $(-\pi/4,\pi/4)$. Exactly noiseless $\phi=\pi/4$ is therefore outside the switching families, but belongs to their closure. Exclusion with a robust margin still requires a confidence range strictly beyond $\pi/4$.

### Proposition T14: nonnegative complete-SR mixtures preserve net-direction sectors

In the task coordinates of one condition, let

$$
H(J)=cI+\int_0^gF_\gamma(J)\,\mu_F(d\gamma)
+\int_0^gF_\gamma(-J)\,\mu_B(d\gamma),
$$

where $c\ge0$, the measures are finite and nonnegative, and the operator is nonzero. All scales and directions share one feature mapping. Then

$$
H(J)=uI+vJ,\qquad
u=c+\int_0^g\frac{\mu_F+\mu_B}{1+\gamma^2}>0,
$$

$$
v=\int_0^g\frac{\gamma(\mu_F-\mu_B)}{1+\gamma^2},
\qquad |v|\le gu.
$$

Consequently, each condition's net phase $\theta=\arctan(v/u)$ lies in $[-h,h]$. Classifying the signs of the two net phases again gives the four sectors in T13. A relative angle greater than $h$ requires both conditions to be net forward relative to their respective tasks.

**Proof.** Expanding each complete SR gives the coefficients. Nonnegative weights and $\gamma\le g$ yield the inequality termwise. The relative angle is the sum of the two net phases expressed in their own task coordinates; interval addition follows T13. ∎

This identifies weighted net direction, not the absence of opposite-direction components. It does not cover the arbitrary relative family $cI+\alpha N_\gamma$, whose scalar coefficient can be negative, as T07 demonstrates.

### Theorem T15: sharp partial identification of the two discounts

In the pure complete-SR FF model, let $\gamma,\eta\in[0,g]$, $h=\arctan g$, and suppose the exact relative angle is $\phi\in[0,2h]$, with free positive gain ratio. Then

$$
\gamma\in
\left[
\tan\{\max(0,\phi-h)\},
\tan\{\min(h,\phi)\}
\right].
$$

The marginal set for $\eta$ is identical, and the joint compatible set is the curve $\arctan\gamma+\arctan\eta=\phi$. Every allowed value is attainable, so the bounds are sharp.

**Proof.** Put $\theta_1=\arctan\gamma$ and $\theta_2=\arctan\eta$. The constraints $0\le\theta_i\le h$ and $\theta_1+\theta_2=\phi$ are equivalent to

$$
\max(0,\phi-h)\le\theta_1\le\min(h,\phi).
$$

Monotonicity of tangent gives the interval. For any allowed $\theta_1$, choose $\theta_2=\phi-\theta_1$ and use the gain freedom in T13 to match the observed scale. ∎

For the original open domain $[0,1)$, the exact cases are

$$
\begin{array}{ll}
0\le\phi<\pi/4:
&\text{each marginal set is }[0,\tan\phi],\\
\phi=\pi/4:
&\text{each marginal set is }(0,1),\\
\pi/4<\phi<\pi/2:
&\text{each marginal set is }(\tan(\phi-\pi/4),1).
\end{array}
$$

In particular,

$$
\boxed{\phi>\pi/4\quad\Longrightarrow\quad
\gamma>\tan(\phi-\pi/4),\quad
\eta>\tan(\phi-\pi/4).}
$$

For example, $\phi=70^\circ$ requires both discounts to exceed $\tan25^\circ\approx0.4663$, but generally identifies neither uniquely. With a finite upper bound $g$, boundary angles such as $0$ or $2h$ can uniquely identify both discounts. General non-uniqueness must not be used to erase these special cases.

For the mixtures in T14, the same phase reasoning constrains only effective net coefficients $q_i=v_i/u_i$, not every component discount. Under noise, the entire lawful joint compatible set must be projected. Substitution of a point-estimated angle does not itself produce a confidence bound.

### Proposition T16: one rotation plane aliases same-direction multiscale mixtures with one effective discount

Let $H=cI+\int F_\gamma(J)\mu(d\gamma)$ with $c\ge0$, $\mu\ge0$ and $H\ne0$, and write $H=uI+vJ$. Define

$$
\gamma_{\rm eff}=v/u,\qquad
g_{\rm eff}=u(1+\gamma_{\rm eff}^2).
$$

Then $H=g_{\rm eff}F_{\gamma_{\rm eff}}(J)$. If all component discounts are below one, $0\le\gamma_{\rm eff}<1$.

**Proof.** We have $u>0$ and $v\ge0$. The ratio $v/u$ is a discount average weighted by $1/(1+\gamma^2)$, with the current component adding mass at zero, so it is lawful. Substitution into $F_{\gamma_{\rm eff}}(J)=(I+\gamma_{\rm eff}J)/(1+\gamma_{\rm eff}^2)$ proves equality. ∎

Additional templates within the same plane can test a shared linear action, but cannot separate temporal kernels that already give exactly the same planar operator. Quantitative identification requires additional propagation shapes or task blocks.

## A7. Constructive routes to discount and temporal-kernel identification

### Theorem T17: three propagation shapes exactly characterize discount identification with a free current component

Let $P\in\mathbb R^{m\times m}$ be known with spectral radius at most one, and let the direct reference $B\in\mathbb R^{m\times p}$ be known. Consider the full parameter domain

$$
Y(c,\alpha,\gamma)
=cB+\alpha P(I-\gamma P)^{-1}B,
\qquad c\in\mathbb R,\quad\alpha>0,\quad0\le\gamma<1.
$$

The following are equivalent:

1. $\gamma$ is globally structurally identifiable within this family.
2. All three parameters $c,\alpha,\gamma$ are jointly identifiable.
3. The matrices $B,PB,P^2B$ are linearly independent under real scalar combinations.

If the third condition fails, the observation at every lawful parameter point can be reproduced at every other lawful discount. Non-identifiability is not confined to one exceptional collision.

**Sufficiency.** Suppose two parameter triples give the same $Y$. Denote the second discount by $\delta$ and its other parameters by $c',\alpha'$. All matrix factors are functions of $P$. Premultiplying by $(I-\gamma P)(I-\delta P)$ gives

$$
\begin{aligned}
0={}&(c-c')B\\
&+\{\alpha-\alpha'-(c-c')(\gamma+\delta)\}PB\\
&+\{(c-c')\gamma\delta-\alpha\delta+\alpha'\gamma\}P^2B.
\end{aligned}
$$

Linear independence requires all three coefficients to vanish. The first gives $c=c'$, the second gives $\alpha=\alpha'$, and the third becomes $\alpha(\gamma-\delta)=0$. Since $\alpha>0$, $\gamma=\delta$, proving joint uniqueness.

**Necessity.** If $B=0$, all observations are zero. If $B\ne0$ but $PB=\lambda B$, then

$$
Y=\left(c+\frac{\alpha\lambda}{1-\gamma\lambda}\right)B.
$$

Any change of discount can be absorbed by free $c$, while keeping $\alpha$ positive.

The remaining case has independent $B,PB$ but real numbers $s,t$ such that

$$
P^2B=sB+tPB.
$$

Thus $\operatorname{span}\{B,PB\}$ is invariant under left multiplication by $P$. In the ordered basis $(B,PB)$, this action is

$$
C=\begin{pmatrix}0&s\\1&t\end{pmatrix}.
$$

Its spectrum is that of the restriction of left multiplication by $P$ to this invariant subspace. Hence $I-\gamma C$ is invertible for $0\le\gamma<1$. Define

$$
D_\gamma=\det(I-\gamma C)=1-t\gamma-s\gamma^2.
$$

$D_\gamma$ is continuous and nonzero on this connected interval, and $D_0=1$, so it is strictly positive throughout. Explicit inversion gives

$$
P(I-\gamma P)^{-1}B=\frac{\gamma sB+PB}{D_\gamma},
$$

and therefore

$$
Y=\left(c+\frac{\alpha\gamma s}{D_\gamma}\right)B
+\frac{\alpha}{D_\gamma}PB.
$$

For any alternative discount $\delta\in[0,1)$, set

$$
\alpha'=\alpha D_\delta/D_\gamma>0,\qquad
c'=c+\alpha s(\gamma-\delta)/D_\gamma.
$$

Substitution yields $Y(c',\alpha',\delta)=Y(c,\alpha,\gamma)$, proving necessity and the pointwise alias. ∎

**Interpretation and scope.** Linear independence here means rank three after vectorizing the three matrices, not full row rank of $B$. For example, if $J_3$ is the three-step lower shift and $B=e_1$,

$$
B=e_1,\quad J_3B=e_2,\quad J_3^2B=e_3,\qquad
Y=(c,\alpha,\alpha\gamma)^\top.
$$

A single channel can identify all three parameters because the task separates three propagation shapes. T19 embeds this block in a lawful stochastic task. The example explains the algebraic condition and is not a successful neural experiment.

The theorem assumes known $P$, a known direct $B$, strictly positive $\alpha$, unrestricted real $c$ and the specified discount interval. Changing these assumptions requires a new analysis. If $\alpha=0$ is allowed, all discounts coincide on that boundary. Unknown dynamics or reference measurements cannot be bypassed using this theorem; T05 and T07 still apply. Identification is within a specified directional branch and does not itself identify that branch. Linear independence also does not guarantee stable estimation under finite noise: nearly dependent propagation shapes can produce broad confidence sets.

### Theorem T18: explicit discount recovery from two distinct rotation blocks

Use the directly referenced model of T17. Suppose two visible task blocks are

$$
A_j=a_jI+\sigma\beta_jJ,\quad
a_j^2+\beta_j^2=1,\quad\beta_j>0,\quad a_1\ne a_2,
$$

and both share the same $c,\alpha,\gamma$ with $\alpha>0$. Each reference has a nonzero compatible projection, making $u_j,v_j$ observable. Then

$$
D_j=1-2a_j\gamma+\gamma^2,\qquad
u_j=c+\alpha(a_j-\gamma)/D_j,\qquad
v_j=\sigma\alpha\beta_j/D_j.
$$

Define

$$
q=\frac{v_1/\beta_1}{v_2/\beta_2}=\frac{D_2}{D_1},
\qquad
t=\frac{q-1}{2(qa_1-a_2)}.
$$

Within the lawful model, $t=\gamma/(1+\gamma^2)$ and

$$
\boxed{\gamma=\frac{2t}{1+\sqrt{1-4t^2}}.}
$$

The other parameters follow as

$$
\alpha=\frac{\sigma v_jD_j}{\beta_j},\qquad
c=u_j-\alpha(a_j-\gamma)/D_j.
$$

**Proof.** We have $q=(1-2a_2t)/(1-2a_1t)$, which gives the formula for $t$. Also,

$$
\frac{dq}{d\gamma}
=\frac{2(a_1-a_2)(1-\gamma^2)}{D_1^2},
$$

whose sign is strictly fixed over $[0,1)$, establishing injectivity. The function $t$ increases strictly over $[0,1/2)$. Solving its quadratic equation and selecting the root below one gives the inverse. The rationalized expression avoids a $0/0$ form at $\gamma=0$. Substitution into the block coefficients recovers the other parameters. ∎

**A remaining falsifiable constraint.** The two block-specific values of $c$ must agree. Four observable coefficients therefore satisfy three shared parameters, providing a remaining mean-level restriction. For two actual Fourier blocks of an eight-state cycle, choose $(a_1,\beta_1)=(0,1)$ and $(a_2,\beta_2)=(\sqrt2/2,\sqrt2/2)$. Coefficients $v_1=4/5$ and $v_2=(8+10\sqrt2)/17$ force $\gamma=1/2,\alpha=1$. If $u_1=u_2=1$, however, the blocks yield $c_1=7/5$ and $c_2=(19-6\sqrt2)/17$, rejecting the shared model. This example deliberately uses unit-circle blocks realizable in a finite stochastic task; an arbitrary abstract orthogonal block need not correspond to a unit-modulus eigenvalue of a finite Markov chain.

Alternatively, put $s_j=\sigma v_j/\beta_j$. When $s_1\ne s_2$,

$$
\gamma=\frac{s_1a_1-s_2a_2-(u_1-u_2)}{s_1-s_2}.
$$

This expression degenerates at $\gamma=0$, whereas the ratio inverse remains available. The derivative of the ratio tends to zero as $\gamma\to1$, showing reduced sensitivity of that ratio; it does not prove that the full coefficient information degenerates in the same way. If the blocks may choose independent gains, current coefficients or discounts, the shared-parameter result no longer applies. Statistical estimation should invert the joint model with coefficient and reference errors rather than apply a normal approximation directly to unstable ratios.

### Proposition T19: a scalable nilpotent block in a lawful stochastic task

Represent each state by a binary word $x=(x_1,\ldots,x_k)$. Each step deletes the first bit and appends an independent fair bit. The complete graph has $n=2^k$ states, each with two incoming and two outgoing edges of probability $1/2$, and is doubly stochastic. Define

$$
C_{j,x}=n^{-1/2}(1-2x_j),\qquad j=1,\ldots,k.
$$

Then

$$
C\mathbf1=0,\quad CC^\top=I_k,\quad
CP=J_kC,\quad CP^\top=J_k^\top C,\quad J_k^k=0,
$$

where $(J_k)_{i+1,i}=1$.

**Proof.** Across all states, the bits are independent and balanced, giving zero mean and orthogonality. Given a destination, the predecessor's first bit is independent and fair, while its remaining bits equal the destination's first $k-1$ bits. Therefore the first row of $CP$ is zero and its row $j>1$ is row $j-1$ of $C$. Averaging over successors gives the transpose relation. Applying the shift $k$ times gives zero. ∎

It is the prespecified contrast block, not the complete transition matrix, that is nilpotent. The full task remains a stochastic transition system. The graph is the classical de Bruijn construction, used here to separate propagation directions at finitely many lags.

### Theorem T20: finite temporal kernels and mixture tests with a direct reference

In the block of T19, suppose a qualified direct reference supports recovery of the required operator quantities and the target kernel is $f(z)=\sum_h w_hz^h$. Then

$$
T=a\sum_{h=0}^{k-1}w_hJ_k^h+cI.
$$

Define

$$
\psi_h=\frac1{k-h}\sum_{i=1}^{k-h}
(T_{i+h,i}-T_{i,i+h}),\qquad1\le h<k.
$$

We obtain $\psi_h=aw_h$. A backward kernel gives $-aw_h$, and a symmetric operator gives zero. A single complete SR has $w_h=\gamma^h$, giving $\psi_{h+1}/\psi_h=\gamma$ wherever the denominator is nonzero, and $\psi_2^2=\psi_1\psi_3$ when $k\ge4$.

For a nonnegative complete-SR scale mixture with a common feature mapping, $w_h=\int_0^1\gamma^h\mu(d\gamma)$,

$$
w_2^2\le w_1w_3.
$$

If $w_1>0$, define $\nu(d\gamma)=\gamma\mu(d\gamma)/w_1$. Then

$$
\frac{w_1w_3-w_2^2}{w_1^2}
=\operatorname{Var}_\nu(\gamma).
$$

**Proof.** $J_k^h$ is one exactly on the $h$th lower subdiagonal and zero elsewhere, giving the formula for $\psi_h$. Cauchy–Schwarz applied to $\gamma^{1/2}$ and $\gamma^{3/2}$ gives the inequality. Under the normalized measure, $\mathbb E_\nu\gamma=w_2/w_1$ and $\mathbb E_\nu\gamma^2=w_3/w_1$, yielding the variance identity. ∎

The variance can also be expressed using $\psi_h$, eliminating the common positive gain. Mass at zero discount is invisible at positive lags and cannot be excluded by equality. Finitely many visible lags establish neither the infinite tail nor a unique learning algorithm. If each scale uses an arbitrary distinct feature mapping, the moment inequality need not hold.

### Theorem T21: recovering a kernel from an unknown learned reverse reference

In the same nilpotent block, let the target and reverse-reference kernels be

$$
f(J_k)=\sum_{h=0}^{k-1}f_hJ_k^h,\qquad
q(J_k^\top)=\sum_{h=0}^{k-1}q_h(J_k^\top)^h,\qquad q_0\ne0.
$$

If the required relative-operator quantities are observable, then

$$
T=a f(J_k)q(J_k^\top)^{-1}
\quad\Longrightarrow\quad
T_{h+1,1}=\frac{a}{q_0}f_h.
$$

**Proof.** Since $J_k^\top e_1=0$, we have $q(J_k^\top)e_1=q_0e_1$, so its inverse maps $e_1$ to $e_1/q_0$. Multiplication by $f(J_k)$ gives the first-column identity. ∎

Thus the visible target kernel is recovered up to a common scale without requiring a geometric reference kernel. Directional interpretation additionally requires $a/q_0>0$. Current components change only $f_0$ or $q_0$, leaving positive-lag ratios intact. When $q_0=0$, the reference block is singular and this argument fails.

For pure complete SRs,

$$
T=a(I-\gamma J_k)^{-1}(I-\eta J_k^\top),
\qquad
a=T_{11},\quad
\gamma=T_{21}/T_{11},\quad
\eta=-T_{12}/T_{11}.
$$

These formulas hold for $k\ge2$. With $k\ge4$, the first column must also satisfy $T_{31}^2=T_{21}T_{41}$. Parameter recovery is conditional on the pure-SR family; the remaining entries and competing kernels must still be evaluated.

The direct subdiagonal-difference statistic in T20 must not be applied unchanged to this relative operator: the reference affects both upper and lower subdiagonals. At $\gamma=0.6,\eta=0.8,a=1,k=4$, an exact single SR gives $\psi_2^2-\psi_1\psi_3=-0.18607104$ under that mistaken application, spuriously resembling a scale mixture. The relative first column is the appropriate quantity for this reference problem. Its observability must follow from T10 or a separate restricted-model argument; a two-dimensional direction readout does not automatically recover it.

### Proposition T22: zero direction does not identify a unique symmetric mechanism

On the $J_4$ block, the average of the two directional complete maps,

$$
\tfrac12\{F_\gamma(J_4)+F_\gamma(J_4^\top)\},
$$

has a constant diagonal. Symmetrizing the transition before constructing an SR gives

$$
F_\eta((J_4+J_4^\top)/2),
$$

which, for $0<\eta<1$, satisfies

$$
[F_\eta]_{22}-[F_\eta]_{11}
=\frac{x^2}{1-3x^2+x^4}>0,\qquad x=\eta/2.
$$

**Proof.** The nilpotent complete maps have diagonal entries one, proving the first assertion. The tridiagonal determinant for the second inverse is $1-3x^2+x^4$, and the relevant cofactors are $1-x^2$ and $1-2x^2$. Their difference gives the stated expression. ∎ Positive gain preserves the difference, while an extra scalar current term cancels. A zero antisymmetric score alone cannot distinguish these mechanisms.

## A8. Limits imposed by temporal history and error budgets

### Theorem T23: a backward alias from measurement history and its control by independent probes

Suppose the true instantaneous response is only the current-state code $b(S_t)$, but the measurement is

$$
y_t=\sum_{h\ge0}\lambda^h b(S_{t-h}),\qquad0\le\lambda<1.
$$

Under a stationary Markov chain,

$$
\mathbb E[y_t\mid S_t=i]=[M_\lambda(P^*)B]_i.
$$

**Proof.** The conditional past distribution is $\Pr(S_{t-h}=j\mid S_t=i)=[(P^*)^h]_{ij}$. Take conditional expectations termwise and sum. ∎ Reversing the task also reverses this measurement artifact; task reciprocity alone therefore does not exclude a history-filter mechanism.

Suppose instead that post-learning probes $Z_t$ are iid and independent of the learning sequence, the learned response $r_c(z)$ is stable during measurement, and the filter is absolutely summable. Then

$$
y_t=\sum_{h\ge0}h_h r_c(Z_{t-h})
\quad\Longrightarrow\quad
\mathbb E[y_t\mid Z_t=z]
=h_0r_c(z)+\sum_{h\ge1}h_h\mathbb E r_c(Z).
$$

The second term is state independent and disappears under state contrasts. This is a conditional-mean control. A directional interpretation also requires the effective current gain to agree with the observation model; a zero coefficient destroys visibility. Sampling without replacement, prohibiting adjacent repeats, state-dependent adaptation, continuing learning and unstable measurements do not satisfy this simple proof.

### Theorem T24: a deterministic perturbation bound for complete operator recovery

Suppose

$$
\widehat X_0=X_0+E_0,\quad
\widehat X_1=TX_0+E_1,\quad
\|E_i\|_F\le e_i.
$$

If $\widehat X_0$ has full row rank and $s=\sigma_{\min}(\widehat X_0)>e_0$, define $\widehat T=\widehat X_1\widehat X_0^+$. Then

$$
\|\widehat T-T\|_F\le
\delta_F=\frac{e_1+\|\widehat T\|_2e_0}{s-e_0}.
$$

**Proof.** Full row rank gives

$$
\widehat T-T=(E_1-TE_0)\widehat X_0^+.
$$

Writing $d=\|\widehat T-T\|_F$, we have $d\le(e_1+\|T\|_2e_0)/s$ and $\|T\|_2\le\|\widehat T\|_2+d$. Rearranging yields $d(s-e_0)\le e_1+\|\widehat T\|_2e_0$. ∎ This proof does not require $\widehat X_1-\widehat T\widehat X_0=0$.

For a prespecified marker $W$,

$$
|\langle\widehat T-T,W\rangle_F|\le\delta_F\|W\|_F.
$$

The uncertainty radius is therefore $\delta_F/\|K\|_F$ for $\omega$, $\delta_F/\sqrt8$ for the normalized eight-cycle marker and $\delta_F\sqrt{2/(k-h)}$ for $\psi_h$ in T20. A direction is declared only if a valid interval excludes competing signs. Rank failure is not evidence for zero direction.

The result is conditional on the actual perturbations satisfying the bounds. Statistical interpretation requires valid joint coverage of those input bounds. Truncated pseudoinverses, reference drift, transition uncertainty and systematic model discrepancy are not automatically included in $e_i$. Floating-point tolerances or conditioning safeguards are not a proof using validated outward rounding.

## A9. Continuous-family separation, joint uncertainty and abstention

### Proposition T25: worst-case family separation and local exclusion

Under a prespecified observation, scale and parameter domain, let $\mathcal M_F,\mathcal M_B$ be nonempty compact prediction sets and define

$$
d_{FB}=\inf_{f\in\mathcal M_F,\;b\in\mathcal M_B}\|f-b\|.
$$

If the allowed perturbations are complete closed norm balls of radii $\varepsilon_F,\varepsilon_B$, their error neighbourhoods are disjoint if and only if

$$
d_{FB}>\varepsilon_F+\varepsilon_B.
$$

**Proof.** If a point belongs to both neighbourhoods, the triangle inequality places the corresponding predictions at distance at most the sum of radii. Conversely, compactness ensures that a pair $f,b$ attains the minimum distance. If that distance does not exceed the sum of radii, a point on the segment joining them lies within the respective radius of each endpoint. The neighbourhoods therefore intersect. ∎ Without compactness, equality when an infimum is not attained requires separate treatment.

This distance is a worst-case design guarantee. Parameters approaching zero discount, zero direction or zero amplitude may force $d_{FB}=0$, while a neighbourhood of one actual observation can still exclude an entire competitor. Local inference should intersect the same valid uncertainty region $\mathcal C$ with each family; it need not demand positive separation everywhere between the full families. Nor does substituting an estimated reference $\widehat B$ into a conditional distance automatically provide statistical coverage.

**Procedure-specific sufficient conditions.** Let $g$ be the distance from the true point to the wrong family and $e$ the observation error. Nearest-family classification is guaranteed correct when $g>2\|e\|$. If the correct class instead uses a training-frozen prediction at distance $h$ from truth, a sufficient condition becomes $g>2\|e\|+h$.

Indeed, the wrong family is at least $g-\|e\|$ from the observation, while the true point is at most $\|e\|$ away and the frozen correct prediction at most $\|e\|+h$. Whole-family geometry cannot therefore justify a frozen-point procedure while ignoring its transfer error. These statements are deterministic geometry, not additional sampling probabilities.

### Theorem T26: candidate-mean residual inversion and error control for one property declaration

Conditional on independent training data, freeze state labels, task blocks, channel processing, templates, ROI, time window and the full candidate space. Suppose the formal observations are $n$ complete paired vectors

$$
Z_i=(B_i^\top,Y_i^\top)^\top
$$

that are iid Gaussian under the conditional training distribution. For each prespecified candidate $\theta$, define the $q$-dimensional residual

$$
E_i(\theta)=Y_i-T_\theta B_i.
$$

Assume that the true $\theta_0$ belongs to the candidate space, satisfies $\mathbb E Y_i=T_{\theta_0}\mathbb E B_i$, has positive-definite residual covariance, and $n>q$. Let $S_{E,\theta}$ be the unbiased sample covariance and set

$$
R(\theta)=\bar E_\theta^\top S_{E,\theta}^{-1}\bar E_\theta,
\qquad
K_q(\alpha)=\frac{q(n-1)}{n(n-q)}F_{q,n-q}(1-\alpha).
$$

Invert all continuous candidates to form

$$
\widehat\Theta_\alpha
=\{\theta:R(\theta)\le K_q(\alpha)\}.
$$

The true parameter is retained with probability $1-\alpha$. If a prespecified property $\ell(\theta)$ is declared only when the set is nonempty and every compatible candidate has that property, then

$$
\Pr_{\theta_0}\{\text{declare and err}
\mid D_{\rm train}\}\le\alpha.
$$

**Proof.** The true residuals are iid Gaussian with zero mean. The classical Hotelling mean distribution gives

$$
\frac{n(n-q)}{q(n-1)}R(\theta_0)\sim F_{q,n-q}.
$$

Hence $\Pr(\theta_0\in\widehat\Theta_\alpha\mid D_{\rm train})=1-\alpha$. If truth is retained, a property common to all compatible candidates must agree with the true property. The event of an erroneous declaration is contained in the event that the true parameter is excluded, and therefore has probability at most $\alpha$. ∎

The sampling distribution is classical. Its role here is to connect uncertainty in both reference and target to continuous SR families and mechanistic labels. FF, BB, FB and BF should be represented as a labelled disjoint union of parameter spaces. If different labels map to the same numerical $T$, all must be retained. No additional Bonferroni factor over the number of families or infinitely many parameter points is required by this argument, which controls exclusion of a fixed true parameter. Post hoc choice of windows, ROIs, templates, state permutations or references is a separate selection problem requiring independent freezing or appropriate adjustment.

The controlled probability is not $\Pr(\text{wrong}\mid\text{a declaration was made})$, a candidate's posterior probability, or simultaneous coverage of every observationally equivalent parameter. It assumes the true mechanism and common-reference relation are in the candidate space. Model omission and systematic drift cannot be repaired by a nominal $p$ value.

**Joint reference–target covariance.** With one template, $B_i,Y_i\in\mathbb R^2$. If

$$
S_Z=\begin{pmatrix}S_B&S_{BY}\\S_{YB}&S_Y\end{pmatrix},
$$

then

$$
\bar E_T=\bar Y-T\bar B,\qquad
V_T=S_Y+TS_BT^\top-TS_{BY}-S_{YB}T^\top.
$$

This is exactly the residual sample covariance; the mean-residual covariance requires an additional factor $1/n$. Within-pair dependence is fully allowed. Dependence across pairs is not covered by the theorem.

**Multiple templates.** For $r$ fixed two-dimensional templates, retain every reference and target within the same real repeated unit. The residual for a common operator generally has dimension $q=2r$, uses $I_r\otimes T_\theta$ and the complete cross-covariance, and requires threshold $K_{2r}$. Templates are not extra independent participants. Known deterministic redundancy may be removed in advance. Arbitrary truncation of the evaluation covariance does not license an otherwise unproved ordinary $F$ threshold.

### Proposition T27: a two-dimensional residual region is not a four-dimensional joint-mean region

For one template with $S_Z\succ0$, let $L_T=(-T,I)$. Then

$$
\inf_{\mu:\,\mu_Y=T\mu_B}
(\bar Z-\mu)^\top S_Z^{-1}(\bar Z-\mu)
=
(\bar Y-T\bar B)^\top
(L_TS_ZL_T^\top)^{-1}
(\bar Y-T\bar B)
=R(T).
$$

**Proof.** Set $h=\bar Z-\mu$, so the constraint becomes $L_Th=L_T\bar Z$. Lagrange optimization of the strictly convex quadratic gives

$$
h_*=S_ZL_T^\top(L_TS_ZL_T^\top)^{-1}L_T\bar Z.
$$

Substituting gives the identity. ∎

The two-dimensional residual threshold

$$
K_2(\alpha)=\frac{2(n-1)}{n(n-2)}F_{2,n-2}(1-\alpha)
$$

can therefore be implemented through intersection with a four-dimensional **auxiliary** ellipsoid. That auxiliary ellipsoid is not a $1-\alpha$ confidence region for the four-dimensional mean. The genuine joint-mean region is

$$
\mathcal C_{1-\alpha}
=\{\mu:(\bar Z-\mu)^\top S_Z^{-1}(\bar Z-\mu)\le K_4(\alpha)\},
$$

$$
K_4(\alpha)=\frac{4(n-1)}{n(n-4)}F_{4,n-4}(1-\alpha),
\qquad n>4.
$$

Its coverage follows from the four-dimensional Hotelling distribution. Intersecting one such $\mathcal C$ with all candidate mean-relation sets protects properties implied by the true mean relation. Any additional arbitrary mean property or reference-drift claim using this region-based route must use $K_4$ or another genuinely valid joint region.

The two-dimensional residual theory only needs $n>2$ and invertible residual covariance. An implementation using positive-definite four-dimensional auxiliary geometry needs $n>4$; these conditions must not be conflated. As an arithmetic illustration, at $n=17,\alpha=0.05$, $K_2\approx0.4621$ and $K_4\approx0.9206$. These numbers do not certify independence of any actual set of 17 pairs.

### Proposition T28: directional-area bias and sensitivity to bounded reference drift

Write a paired observation as $B_i=b+\varepsilon_{B,i}$ and $Y_i=y+\varepsilon_{Y,i}$, with zero-mean errors and cross-covariance $\Sigma_{BY}$. Then

$$
\mathbb E\det(B_i,Y_i)
=\det(b,y)+(\Sigma_{BY})_{12}-(\Sigma_{BY})_{21}.
$$

For $n$ iid pairs,

$$
\mathbb E\det(\bar B,\bar Y)
=\det(b,y)+
\frac{(\Sigma_{BY})_{12}-(\Sigma_{BY})_{21}}n.
$$

**Proof.** Expand the determinant. Linear error terms have zero expectation, and the two cross-products contribute the respective cross-covariances. Cross-covariance of sample means is divided by $n$. ∎ Averaging per-repeat areas or simply substituting sample means does not automatically produce an unbiased mechanistic statistic. T26 instead tests the population mean residual without treating the bilinear area as an unbiased univariate outcome.

If the actual target reference is $b+\delta$ and external evidence guarantees $\|\delta\|\le\rho$, then

$$
\inf_{\|\delta\|\le\rho}\det[b+\delta,y]
=\det[b,y]-\rho\|y\|.
$$

For $y\ne0$, the smallest reference change reaching the zero-direction boundary is

$$
\rho_*=\frac{|\det[b,y]|}{\|y\|}.
$$

**Proof.** We have $\det[\delta,y]=\delta^\top(y_2,-y_1)^\top$, whose coefficient vector has norm $\|y\|$. Cauchy–Schwarz gives the extremum, attained by choosing $\delta$ in the opposite direction. Solving the zero-determinant linear constraint with minimum norm yields the second formula. ∎

With sampling error as well, a genuine joint-mean confidence region $\mathcal C$ can be used to check

$$
\inf_{(b,y)\in\mathcal C}
\{\det[b,y]-\rho\|y\|\}>0.
$$

If both coverage and the external drift bound are valid, this guarantees positive direction under the permitted drift. If $\rho$ is merely a scanned sensitivity parameter, the output remains a sensitivity analysis. It neither establishes that real drift has been covered nor permits the $K_2$ auxiliary ellipsoid of T27 to be treated as a joint confidence region.

## A10. The proper use of each conclusion in an empirical paper

| Question | Main results | Conditions that cannot be omitted |
|---|---|---|
| Can canonical distances identify temporal direction? | T02–T03 | Specify distance, metric, states and readout; the comparison uses a common discount |
| Does full relative-operator fitting automatically resolve reference ambiguity? | T04, T07 | Learned reference structure and free current components must enter the competitor set |
| Is a direction property determined by the data? | T05–T06, T13–T14 | Every lawful compatible explanation must share it; retain switching and mixture labels |
| Which reciprocal experiment supports a direction claim? | T08–T09, T11 | Common mapping, anchored dynamics and the stated mechanism family |
| How much observation rank is sufficient? | T10, T12 | Separate full operators, linear functionals and restricted properties; assess true projection visibility |
| What can be quantified after direction is identified? | T15, T17–T21 | Sharp parameter sets or additional propagation shapes, blocks and finite-kernel restrictions |
| Can one plane establish a single SR kernel? | T06, T16 | No; extra templates test common action, whereas kernel identification needs additional task information |
| Can temporal persistence create a backward pattern? | T23 | Independently examine the history-filter alternative and actual probe dependence |
| When can a statement be made under error? | T24–T28 | Valid error bounds or sampling coverage; exclusion of entire continuous competitors and reliable numerics |

The final procedure distinguishes four outcomes: a nonempty compatible set with one prespecified property; several properties remaining compatible; rejection of every specified model; and insufficient data or inferential qualification. The last includes too few repeated units, unusable covariance, invisible relevant signal, or unsupported coverage assumptions. Under the freely varying single-plane mean union in T06, saturation generally prevents rejection of every model at the mean level. Additional task blocks and parameter sharing provide a genuinely falsifiable kernel-fit check.

Independent repeats must be defined by the actual sampling unit of a complete paired map, not by recombined occurrence counts, numbers of templates or imposed computational transformations. Distinct event identifiers do not imply disjoint temporal support. Shared preprocessing, slow drift and long-range dependence require examination. Participants with different references and mechanisms cannot simply be pooled as iid repetitions of one individual's mechanism.

## A11. Algebraic reproduction and result index

Programs `code/math_checks.py` and `code/math_additions.py` reproduce the exact constructions, compensation formulas, inverse relations and perturbation examples. Outputs and partial-discount bounds are in `source_data/`. The registry `math/theorem_registry.json` indexes T01–T28. General statements are established by the proofs above; real-measurement methods and provenance are in Part B.

# B. Real repeated measurements and calibration

## E1. Purpose and levels of evidence

The empirical analyses connect the proposed measurement procedure to existing neural recordings. We ask whether a reference pattern predicts later measurements, whether a task-compatible readout preserves directional information under a known transformation, and how acquisition structure and reference transfer affect the resulting decision. Known-transform calibration applies a specified operator to a measured response, providing a ground truth for the transformation while retaining actual measurement variation.

We distinguish three endpoints: reference prediction on actual earlier/later measurements; known-transform calibration under explicit paired-error assumptions; and limitations revealed by the physical acquisition hierarchy or transfer across recordings. Successful prediction, conditional exclusion of a direction, and proof of an SR-specific biological mechanism are separate claims.

## E2. Data, states and acquisition units

The main empirical figure uses a subset of the circular-position EEG localizer from Turner et al. (2025) [E1]. The original study included 18 participants; the present analysis uses the two previously downloaded localizer recordings of one participant, labelled ID13 in the prior analysis report. Each supplied recording contains 1,980 eligible events, 40 position labels and two physical blocks separated by the localizer break. A block is not a participant. We retain the original paper's 17 posterior electrodes: P7, P5, P3, P1, Pz, P2, P4, P6, P8, PO7, PO3, POz, PO4, PO8, O1, Oz and O2.

The early and late windows are 100–150 ms and 200–250 ms after stimulus onset, respectively, with subtraction of a −200–0-ms baseline. The analysis inputs are the derived event arrays in the archived calibration package. Their named windows are `early` and `late`; exact window bounds and the ID13 mapping are recorded in the source statistical appendix rather than embedded as extraction indices or participant fields in the NPZ. The reproducibility unit for the present calibration is therefore the supplied event-level measurement array.

Cross-record stress additionally uses existing repeated localizer measurements from the Kern MEG archive, DOI 10.5281/zenodo.8001755 [E2]. The Ekman et al. localizer provides a supplementary reference-estimation example comprising two selected participants and two regions, V1 and hippocampus [E3]. None of these subsets is presented as the original study's full participant sample.

## E3. Reference estimation with separated training and evaluation

### E3.1 Temporal splits and estimator candidates

Within each physical Turner block, the midpoint of observed event-onset times separates earlier controls from later evaluation measurements. Events within 0.5 s of either side of that midpoint are excluded. Across eight recording-by-block-by-window contexts, there are 491 or 492 control events and 492 or 495 evaluation events. State counts in the control pools may differ, ranging from five to 22; extra eligible controls are not discarded solely to balance states. Training and evaluation event identities do not overlap. The closest onset gap at the split is approximately 1.20 s. This avoids direct boundary reuse, but does not establish temporal independence of EEG errors.

Internal validation holds out one control event per state on each fold and forms training means from the other controls. Candidates include the direct mean, zero pattern, scalar shrinkage, truncated singular-value estimates and Fourier state smoothing that retains complete sine/cosine pairs. Shrinkage coefficients range from zero to one in increments of 0.1. Fourier candidates retain one, two, four, eight or 19 harmonic pairs, or the entire centered state space. Retaining complete pairs does not favor one rotation direction. Selection uses control-only full-matrix mean squared error (MSE), without reading the later evaluation target or transformation labels.

The main empirical figure presents the outer-evaluation results for the internally selected Fourier estimator family. Each point is a block-specific earlier/later prediction endpoint. Estimator development and evaluation are retrospective. Within each fit, the later segment is reserved for evaluation and all parameter selection uses the earlier controls.

### E3.2 Task-plane loss and the zero baseline

Let the state-centered reference estimate be \(\widehat B\in\mathbb R^{40\times17}\), and let \(X\) contain later state-wise means. The fixed first circular harmonic plane is

\[
C=\sqrt{\frac{2}{40}}
\begin{bmatrix}
\cos(2\pi j/40)\\
\sin(2\pi j/40)
\end{bmatrix}_{j=0}^{39}.
\]

The two task-plane losses are

\[
L_C(\widehat B,X)=\frac{\|C(\widehat B-X)\|_F^2}{2\times17},
\qquad L_C(0,X)=\frac{\|CX\|_F^2}{2\times17}.
\]

Block-specific gain is \(G_C=1-L_C(\widehat B,X)/L_C(0,X)\). Zero is defined after removing the component shared by all states: it predicts no state differences, rather than no brain activity. This baseline is necessary because heavy shrinkage can outperform a noisy empirical mean while retaining no useful state-dependent signal.

If the later target satisfies \(X=B_*+\varepsilon\) with \(E[\varepsilon\mid\text{training data}]=0\), the expected difference between two fixed training estimates' prediction losses equals the difference between their squared errors relative to \(B_*\); the common test-noise term cancels. A realized prediction loss nevertheless is neither a deterministic error relative to \(B_*\) nor an upper bound on that error. Drift or changes in conditional reference means alter this interpretation.

### E3.3 Reference-prediction results

| Physical recording and block | Early-window gain over zero | Late-window gain over zero |
|---|---:|---:|
| Recording 1, block 1 | −8.3292% | 16.9238% |
| Recording 1, block 2 | 0.0000% | 6.8747% |
| Recording 2, block 1 | 25.1302% | 32.3970% |
| Recording 2, block 2 | 29.3631% | 57.9619% |

Averaging losses over the four physical blocks before taking their ratio gives an early-window plane MSE of 10.313989 against 12.094698 for zero, a **14.7231%** gain. The corresponding late-window values are 4.386882 and 5.998847, a **26.8713%** gain. The two windows are reported separately and are not pooled as independent samples.

All four late-window blocks improve, whereas early-window improvement occurs only in the second recording. These data support predictability in the specified plane of this subset, not a general claim of within-person stability across recordings. Moreover, late-window full-matrix average MSE remains higher than its full-matrix zero baseline. Predicting a useful low-dimensional reference does not recover the entire state operator. The window with the best reference prediction also need not provide the strongest directional exclusion, because these objectives weight different aspects of error.

The block comparisons describe realized prediction losses. Stability of latent reference precision or equality of expected prediction risk across blocks is a separate statistical question.

## E4. Known-SR transformations and a task-compatible readout

### E4.1 Regrouped state maps and training restrictions

We form a complete state matrix by combining the kth occurrence of each state. The minimum state count in Turner recording 2 is 38, yielding 38 matrices of size 40×17. The first four are used only to fit the reference, sensor whitening and readout. The remaining 34 alternate between reference and target-base measurements, giving 17 pairs. Training uses 160 original events and inference uses 1,360, with disjoint event identities.

This grouping does not create new independent acquisitions. Eighteen of the 38 maps span both physical blocks, and the median temporal span of a map is 133.318 s. Occurrence-map time intervals overlap even though individual events are not reused.

The reference estimator is selected among scalar, SVD and Fourier families by training-only validation. The task may retain the native state order or use a training-derived order. For the latter, the leading two training-reference components, deterministic sign choices and a linear assignment place states on a circular order so that the quarter-cycle plane \(k=10\) retains stronger training contrasts. This is a concrete candidate order for a future experiment; the participant did not relearn that order during the present calibration.

### E4.2 Preserving the operator action

Let \(Z=C_{\rm eff}XW\) be a measurement projected to the fixed state plane and whitened with fixed sensor weights, and let \(R=C_{\rm eff}\widehat BW\) be the training template. With \(J=\left[\begin{smallmatrix}0&-1\\1&0\end{smallmatrix}\right]\), define

\[
\mathcal L_R(Z)=\frac{1}{\|R\|_F^2}
\begin{bmatrix}\langle R,Z\rangle_F\\\langle JR,Z\rangle_F\end{bmatrix}.
\]

**Implementation and theoretical normalization.** The theoretical readout in Supplementary T12 divides by \(\|R\|_F\), whereas the empirical implementation above divides by \(\|R\|_F^2\). With \(s=1/\|R\|_F>0\), the implemented readout is the theoretical readout multiplied by the same training-frozen positive scale \(s\) for both reference and target. Their relative angle is unchanged, their oriented area is multiplied by \(s^2\) without changing sign, and the task-equivariance identity is preserved. For a fixed candidate operator, the residual and its sample covariance become \(r_T'=s r_T\) and \(V_T'=s^2V_T\). Thus, with nondegenerate covariance and in exact arithmetic,

\[
(r_T')^\top(V_T')^{-1}r_T'=r_T^\top V_T^{-1}r_T.
\]

At fixed sample size and nominal threshold, the two normalization conventions therefore yield identical candidate-residual inversion. This equivalence concerns only a **common scale after the other training choices and candidate family have been frozen**. Reselecting whitening, templates, regularization strengths or other processing may change the readout; the identity does not imply that all regularization parameters or fitted results are automatically invariant. Degenerate norms or covariance still require abstention.

The identity \(\mathcal L_R(JZ)=J\mathcal L_R(Z)\) preserves the action of every plane operator \(uI+vJ\). If the template is zero or its specified plane vanishes numerically, the procedure abstains rather than fitting a favorable replacement on evaluation data. Whitening, state order, template and readout are fixed from training, and the same map is applied to reference and target measurements.

### E4.3 What the known transformation changes

A known \(T\) acts on each two-dimensional base measurement \(v_i\) to create \(y_i=Tv_i\). Canonical SR calibration uses \(T=(I-\gamma A)^{-1}\) or \((I-\gamma A^\top)^{-1}\), where \(A\) is the cycle operator on the selected plane. The calibration discounts are .25, .5, .8 and .95. Additional archived conditions allow separate current and predictive gains.

The transformation acts on the entire measured signal-plus-noise vector. No virtual participants or added simulated noise are generated, but the transformed target is also not a newly acquired response to a natural learning task. These conditions test the readout's algebraic behavior and conditional discriminability in measured repeat structure; they do not estimate conventional experimental power with a fixed biological-noise distribution or establish SR-specific biology.

This empirical analysis calibrates a directional rule under a specified reference relationship. Participants did not perform the proposed reciprocal tasks, and the exclusive FF/BB/FB/BF regions of the pure complete-SR families were not empirically validated. A known canonical backward SR specifies the imposed calibration operator. Its label identifies the calibration ground truth; the measured endpoint is directional exclusion under that transformation.

## E5. Joint error, continuous exclusion and three conditional positives

The four-dimensional paired samples \((b_i,y_i)\) provide a joint mean and covariance. The primary rule inverts **two-dimensional paired-residual tests of candidate operators**; its auxiliary ellipsoid is not a simultaneous four-dimensional mean-confidence ellipsoid. Directional boundaries are optimized continuously rather than checked only at a few discount grid points. Conditional coverage requires pairs that are independent and identically Gaussian distributed given training, a correct common untransformed-reference mean relationship, nondegenerate residual covariance for the true candidate, frozen training choices and the required independence between training and evaluation. Exchangeability alone, or distinct event identities, is insufficient.

The main empirical figure reports the minimum quadratic statistic at the zero-direction boundary divided by the critical value. With 17 pairs and nominal \(\alpha=.05/54\), the critical value is 1.447336. The denominator 54 refers to 27 measurement/task contexts multiplied by two template versions. Each context contains further known-transformation conditions, so the full result table is not limited to 54 rows. A ratio exceeding one excludes the relevant direction boundary under the specified error model. The plotted segment is a margin over that threshold, not an uncertainty interval.

| Frozen training readout | Known backward-SR discount | Boundary/critical ratio | Nominal 54-context adjusted p | Phase |
|---|---:|---:|---:|---:|
| Native order, k=1 | .95 | 1.206961 | .020613 | −67.4544° |
| Designed order, k=10 | .80 | 1.156972 | .025294 | −38.6598° |
| Designed order, k=10 | .95 | 1.341880 | .012182 | −43.5312° |

All three reuse the same participant, second recording, early window and 17 pairs; they are not three independent replications. Original-A and B-trained-A versions of the native k=1 case are statistically equivalent through common scaling and appear only once. The two designed-plane positives use the updated training template. Positive scalar shrinkage by itself does not improve a scale-invariant decision.

The corresponding uniform covariance increases needed to reach the nominal boundary are approximately 20.70%, 15.70% and 34.19%. These values describe geometric sensitivity; they are not estimated corrections for actual EEG dependence. All forward-labelled conditions abstain at the same joint threshold, so no symmetric forward/backward accuracy is reported.

The nominal adjustment controls the stated 54-context comparison set under the paired-error model. The analysis belongs to a retrospective development sequence, so this adjustment does not include the full sequence of methodological choices. Temporal dependence, shared upstream processing and reference drift determine whether the conditional sampling model describes an acquisition. The three results are consequently reported as conditional calibration outcomes.

## E6. Physical-block restrictions and cross-record stress

Restricting the designed k=10, canonical backward \(\gamma=.95\) condition to each physical block leaves six and seven complete pairs. Both analyses abstain at the same nominal joint threshold. This restriction changes available repeats, means and covariance estimates, so its outcome cannot be uniquely attributed to temporal dependence. Treating only the two original blocks as independent units instead leaves insufficient units for the required four-dimensional joint covariance. The resulting constraint belongs in the design of future acquisitions.

Cross-record stress preserves the source-trained \(C_{\rm eff}\), \(W\), \(R\) and \(H\), without refitting favorable coordinates in the target recording. Within-record and cross-record comparisons use the same source reference and matched target sample budgets. At the nominal joint threshold, **1,200/1,200** cross-record computed conditions abstain: 570 Turner and 630 Kern conditions. Discounts, transformation families and readouts reuse the same recordings; these counts are neither participants nor independent experiments.

The true reference mean may change between recordings. Consequently, these conditions do not estimate a type-I error rate under a known null. Universal abstention establishes neither reference equivalence nor the absence of SR; it means that the current procedure does not support a unique directional declaration in these tested transport settings.

## E7. Supplementary fMRI reference example

The Ekman example uses previously processed localizer time series and extracted pattern matrices. Each participant/region contains eight positions with eight events each; the first four events per state train the estimator and the last four evaluate it. Actual event times and the known extraction windows support that split, while shared whole-record processing and haemodynamic dependence remain.

S01 V1 selects scalar shrinkage λ=.3. Full-pattern outer MSE is .05297993, versus .06589265 for zero and .06411860 for the direct mean: gains of 19.5966% and approximately 17.37%, respectively. Both participants' hippocampal estimates select zero. S02 V1 scalar and SVD estimates also select zero, while its Fourier estimate gains approximately .5454% over zero. All contexts remain in `source_data/empirical_ekman_reference_all.csv`. This example illustrates why reference prediction needs a zero baseline. It does not diagnose a direction error in the original study or estimate a general magnitude of localizer contamination.

## E8. Numerical verification and reproducible analysis units

The empirical package preserves event-level inputs, source tables and the checks connecting them to the reported results. The portable script `code/empirical_check_fig07_portable.py` refits the 40 Turner reference estimates from the included measured-array excerpt, reconstructs the two frozen readouts used in the positive conditions, and checks the aggregation of prediction gains, boundary ratios and cross-record outcomes. The maximum discrepancy in the reconstructed reference matrices is zero; the maximum discrepancy in loss values is $3.55\times10^{-15}$. Both frozen readouts are reproduced with zero discrepancy. Training-event identities and disjoint training/evaluation event sets also agree with the archived inputs. These computations start from derived measurements; raw-EEG preprocessing and the acquisition-level independence model lie upstream of this reproducibility unit.

For the earlier continuous-family calibration, `code/empirical_check_dfb_tables.py` joins the archived fold-level decisions and observation-class audit, verifies the sufficient inequalities, and recreates the counts reported in E9. Source-file hashes and original archive paths are recorded in `empirical/dFB_provenance.json`, `empirical/historical_source_manifest.json` and `empirical/round1_source_manifest.json`. Verification results appear in `empirical/fig07_portable_checks.json` and `empirical/dFB_table_checks.json`. Figure values are read from the corresponding CSV files.

The empirical endpoints are defined at their actual scale: event-level prediction for reference estimates; computational conditions for known-transform calibration; physical blocks and recordings for acquisition and transfer checks; and participants only for analyses whose sampling unit is a participant. Repetition of a computational condition does not increase the number of independent acquisitions.

## E9. Continuous-family separation in measured Kern and Turner repeats

### E9.1 Data structure and candidate families

This calibration uses one participant from the Kern MEG archive, DSMR101 [E2, E4], and one participant from the Turner EEG archive, ID13 [E1]. Each dataset contributes two response windows, three state spaces and four reused reference splits. The state spaces are a single Fourier plane, two Fourier planes and the full centered state space: dimensions 2, 4 and 11 for Kern, and 2, 4 and 39 for Turner. Thus there are $2\times2\times3\times4=48$ reference-fold settings.

Each setting contains 18 imposed directional transformations, formed by two directions, three discounts $\gamma_0\in\{.2,.5,.8\}$ and three predictive strengths $\rho\in\{.1,.3,1\}$, together with one current-only control. There are 912 fold-level conditions in total, including 864 directional conditions. Four-fold aggregation gives 228 conditions, of which 216 are directional. These counts describe transformations and splits of two participants' recordings.

The fitted candidates are

$$
T_j(c,\alpha,\gamma)S=cS+\alpha N_\gamma(A_j)S,
\qquad N_\gamma(A)=A(I-\gamma A)^{-1},
$$

with a free current coefficient $c\in\mathbb R$, predictive gain $\alpha\geq0$, and $\gamma\in[0,.95]$. Forward and backward candidates use $A$ and $A^\top$; the symmetric-transition candidate uses $(A+A^\top)/2$. A current-only model is also included. Two symmetric formulas that produce the same predictions on these cycle tasks are merged into one observation class. A fit with $\alpha\leq10^{-12}$, or with a prediction observationally equivalent to current-only, is assigned to the current class regardless of the candidate's nominal directional label.

Let $R,V$ denote the training reference and training target-base measurement, and $S,Z$ their test counterparts. The known operator $T_0$ produces targets $T_0V$ and $T_0Z$. Both the signal and the measurement fluctuation are transformed. The imposed strength is defined on the training target-base input:

$$
\rho=\frac{\alpha_0\|Q_VN_{\gamma_0}(A_j)V\|_F}{\|V\|_F},
$$

where $Q_V$ removes the vectorized current-pattern direction. On a different reference $S$, the predictive residual is computed using the same frozen $\alpha_0$; its relative magnitude need not remain $\rho$. Sensor whitening is estimated from the training reference. Kern's upstream signal-space projection was estimated from the complete first localizer, so the four subsequent splits share that preprocessing.

### E9.2 A fixed true point and a continuous competing family

Allowing $\alpha=0$ makes every directional family intersect the current family. The minimum distance between the unrestricted forward and backward families is therefore zero. The useful quantity here is the distance from a specified true prediction with a nonzero predictive component to the entire competing family.

For a fixed reference $S$, write

$$
r_S=\alpha_0\|Q_SN_{\gamma_0}(A_j)S\|_F,
\qquad
u_0=\frac{Q_SN_{\gamma_0}(A_j)S}{\|Q_SN_{\gamma_0}(A_j)S\|_F}.
$$

For a competing family $k$, let $\nu_k(\gamma)$ be its similarly normalized, current-residualized predictive direction. Optimizing the current coefficient over $\mathbb R$ and the predictive gain over $[0,\infty)$ gives

$$
m_{j\rightarrow k,S}
=r_S\sqrt{1-\left[\max\left\{0,
\sup_{\gamma\in[0,.95]}\langle\nu_0,\nu_k(\gamma)\rangle_F
\right\}\right]^2}.
$$

When a competing predictive residual vanishes, that candidate belongs to the current observation class and its distance is $r_S$; normalization is applied only to nonzero residuals. The binary margin compares forward with backward. The multiclass margin also includes the symmetric and current-only alternatives and takes the minimum distance. Continuous optimization supplies a lower bound $\underline m_S$ on this point-to-family separation. The calculations use numerical analytic interval bounds with independent witness checks; they are floating-point bounds, without outward-rounded interval arithmetic. The finite discount grid used for fitting and the continuous domain used for exclusion are recorded separately.

### E9.3 The error that can change a model comparison

All candidate predictions in a setting lie in the common space

$$
U_S=\operatorname{span}\{\operatorname{vec}(S),
\operatorname{vec}(AS),\ldots,\operatorname{vec}(A^{q-1}S)\},
$$

where $q$ is the number of state contrasts. For these cycle operators, polynomial reduction also includes the backward and symmetric candidates. The space has dimension at most $q$, although a full measurement has $q$ times the number of channels coordinates. For every candidate $f\in U_S$,

$$
\|Y-f\|^2=\|\Pi_{U_S}Y-f\|^2+
\|(I-\Pi_{U_S})Y\|^2.
$$

The last term is common to all candidate scores. Its removal preserves their pairwise squared-error differences and the nearest-family decision exactly.

The calibration center is $\mu_S=T_0S$, so the realized discrepancy of the transformed test measurement is

$$
E=T_0(Z-S),\qquad T_0Z=\mu_S+E.
$$

Across the 864 directional conditions, the median ratios are

$$
\frac{\|\Pi_{U_S}E\|}{\|E\|}=.596376,
\qquad
\frac{\|Q_S\Pi_{U_S}E\|}{\|E\|}=.191525.
$$

These are ratios of error norms, not percentages of explained variance or reductions in prediction MSE. The second projection is appropriate when each candidate can refit its current coefficient freely over $\mathbb R$ on the target being classified. A current coefficient frozen from training defines a different decision rule and requires its own transfer-error term.

Independent SVD and spectral constructions of $\Pi_{U_S}$ agree within $6.01\times10^{-13}$ across the 48 reference settings. The corresponding score-difference identity agrees within $1.25\times10^{-12}$ across the 4,560 archived candidate scores. The archived input reconstruction reproduces training and test SSE within $3.98\times10^{-13}$ and $2.28\times10^{-13}$, respectively.

### E9.4 Sufficient conditions for the matching nearest-family rule

In target-wise nearest-family classification, each candidate refits $c$, $\alpha$ and $\gamma$ against the same reference $S$ and target $T_0Z$. The fitting grid has 96 discounts from 0 to .95 in steps of .01 and contains every imposed true discount. Consequently, the correct fitting family contains the true prediction. Competitor exclusion uses the larger continuous family.

If $\varepsilon$ bounds the error in the score-relevant space, the correct prediction lies within $\varepsilon$ of the target and every competing prediction is at least $\underline m_S-\varepsilon$ away. Thus $\underline m_S>2\varepsilon$ is sufficient for a correct nearest-family decision. Using the realized calibration error gives the following deterministic check:

| Error norm used in the sufficient condition | Binary passes / 864 | Multiclass passes / 864 | Incorrect matching-rule decisions among passes |
|---|---:|---:|---:|
| Full measurement, $\|E\|$ | 0 | 0 | 0 |
| Common prediction space, $\|\Pi_{U_S}E\|$ | 6 | 6 | 0 |
| Common space after removal of current, $\|Q_S\Pi_{U_S}E\|$ | 80 | 77 | 0 |

The smallest positive inequality margin is approximately $6.74\times10^{-4}$, well above the observed projection discrepancies. These counts are a post hoc check of a deterministic sufficient condition with a known transformation and a computable error. They have no nominal confidence level and are distinct from the 17-pair residual inference in E5.

The same conditions also reveal why a guarantee must name its decision rule. Under the earlier rule, parameters are estimated on $R,T_0V$ and then frozen for prediction on $S$. Using the same observation classes for both rules gives:

| Conditions satisfying the target-wise bound | Frozen rule: correct visible direction | Frozen rule: current-only | Frozen rule: symmetric | Target-wise rule: correct visible direction |
|---|---:|---:|---:|---:|
| Binary, 80 conditions | 42 | 38 | 0 | 80 |
| Multiclass, 77 conditions | 39 | 37 | 1 | 77 |

The 37 current-only outcomes and one symmetric outcome in the second row are failures to recover a unique direction; they are not 38 forward/backward reversals. The two procedures answer different questions: prediction with transferred parameters, and classification after fitting the target. Their difference therefore illustrates a missing error term, rather than an improvement in independent prediction accuracy.

Specifically, with $\widehat T_{\rm true}$ denoting the correctly labelled family's training fit,

$$
q_S=\|\widehat T_{\rm true}S-T_0S\|
$$

measures parameter-transfer error. A sufficient bound for the frozen rule is

$$
\underline m_S>2\|\Pi_{U_S}E\|+q_S.
$$

No directional condition passes this conservative bound. The tighter realized-residual bound

$$
\underline m_S>\|\Pi_{U_S}E\|
+\|\Pi_{U_S}(T_0Z-\widehat T_{\rm true}S)\|
$$

passes 10 binary and 10 multiclass fold-level conditions, all correctly classified. No four-fold aggregated condition passes either frozen-rule bound. A sufficient condition can abstain while a particular decision happens to be correct; its purpose is to identify a region in which that decision is assured by the stated error comparison.

### E9.5 Predicting unknown errors and resolving discounts

A future target does not reveal $T_0(Z-S)$. The archived analysis therefore also examined empirical scales derived from training-pair differences and reference repeats. These scales describe observed variation; they are not confidence bounds on the next target's error. Their ability to rank successful frozen-rule conditions differed between the two datasets. For the multiclass comparison, the pair-plus-transfer-error index had AUROC .537 for Kern and .901 for Turner. These values describe rankings of reused computational conditions. Participant and dataset are confounded because each dataset contributes one participant.

Cross-dataset threshold exploration selected a source threshold with at least ten passing conditions and at least 90% observed correct decisions, then applied it unchanged to the other dataset. Only two of 20 direction-by-index-by-transfer combinations supplied such a source threshold. They changed from 12/12 to 9/14 and from 9/10 to 27/36 correct conditions on transfer. Neither retained the source criterion. The corresponding tables preserve arithmetic ties and numerical bounds on rankings; those bounds describe optimization precision, not sampling confidence. A separate THINGS-fMRI calibration used disjoint groups of six sessions for threshold construction and evaluation, and likewise did not yield the specified transferable low-error threshold. This analysis is catalogued under D05 in E11.

Discount resolution poses an additional question: can the true point be separated from same-direction candidates with $|\gamma-\gamma_0|\geq.2$? In a single rotation plane, removing the free current term leaves only one predictive direction. The free predictive gain absorbs the effect of $\gamma$, so this margin is exactly zero. Two planes can create a nonzero shape difference:

| Dataset and state space | Median unit discount-margin lower bound | Frozen fits with reportable $\gamma$ / 144 | Reportable and $\lvert\widehat\gamma-\gamma_0\rvert<.2$ / 144 |
|---|---:|---:|---:|
| Kern, one plane | 0 | 0 | 0 |
| Kern, two planes | .234863 | 124 | 52 |
| Kern, full 11-dimensional space | .228635 | 120 | 43 |
| Turner, one plane | 0 | 0 | 0 |
| Turner, two planes | .081264 | 112 | 27 |
| Turner, full 39-dimensional space | .234912 | 144 | 71 |

Single-plane numerical upper bounds of approximately $2.11\times10^{-8}$ reflect rounding around the exact zero result. Positive multi-plane margins establish geometric resolution, but no discount condition passes the realized-error sufficient bound in these measurements. A reportable fitted value and a guaranteed discount resolution are therefore separate endpoints. Source files for E9 have the prefix `source_data/dFB_`.

## E10. Reference content, physical history and the interpretation of temporal weights

### E10.1 A learned reference can contain temporal organization

The Ekman analysis separates the content of a post-learning reference from its temporal direction [E3]. For each of 35 participants, the author's $8\times8$ localizer probability matrix is reduced to its 56 off-diagonal cells. The regression controls cue-position preference, output-position preference, the four physical distances on the circle, and mean level for pairs of sequence positions. Two additional contrasts represent temporal proximity, coded by negative absolute distance within the learned sequence, and direction, coded $+1$ for future and $-1$ for past sequence positions. The contrasts are defined among the four learned positions. Adding them increases design rank from 18 to 20.

The primary hippocampal input is the HPa array actually called by the source Figure 6 notebook. HPa and the HP array used for Figure 5 are separate derived outputs; the filename alone does not identify an anatomical subdivision. Participant-level bootstrap intervals use 20,000 resamples, and sign-flip tests use 100,000 permutations. The four primary comparisons are the two contrasts in hippocampus and V1, with Holm adjustment.

| Primary source and contrast | Mean coefficient | Participant bootstrap 95% interval | Holm-adjusted sign-flip $p$ |
|---|---:|---:|---:|
| Hippocampal Figure 6 array: temporal proximity | .0126685 | [.0049010, .0206462] | .0164798 |
| Hippocampal Figure 6 array: direction | −.0004862 | [−.0027173, .0016285] | 1 |
| V1: temporal proximity | −.0017576 | [−.0052556, .0018532] | 1 |
| V1: direction | −.0002642 | [−.0016733, .0011401] | 1 |

Temporal-proximity coefficients have probability units per sequence step; directional coefficients have probability units for the specified signed contrast. The result supports temporal organization in the hippocampal reference while leaving its directional component unresolved. The participant is the resampling unit; the 56 cells are within-participant observations. The source table is `historical_localizer_unique_group.csv`.

### E10.2 Physical timing and response overlap

For the two available continuous-fMRI participants, actual event times define a linear mixing matrix under a specified haemodynamic response function. With a response mean of 6 s, standard deviation of 3 s, an unfiltered linear time-invariant calculation and the 3–13.5-s extraction window, the largest off-diagonal contribution relative to the mean diagonal contribution is 20.75% for S01 and 20.98% for S02. The corresponding centered seven-dimensional mixing matrices retain full contrast rank, with condition numbers 1.298 and 1.405. These are conditional properties of the timing schedule and assumed response kernel, rather than estimates of the proportion of neural activity caused by the preceding stimulus.

Predictive evaluation asks a different question. In held-out rounds of actual V1 measurements, adding previous-stimulus identity to current identity changes pooled SSE by −3.283% for S01 and +.616% for S02, where a positive value means improvement. The eight folds within a participant overlap in training data. Thus the timing-based potential for overlap does not by itself determine the observed benefit of a history model. Source tables retain alternative response kernels, filtering choices and extraction windows: `historical_actual_history_mixing_summary.csv` and `historical_localizer_history_prediction_summary.csv`.

### E10.3 Fine-time-scale history in spiking measurements

The Allen analysis uses one mouse and one Visual Coding Neuropixels session, including 118 natural images and a blank condition. The fixed VISp subset contains 110 quality-selected units. From 5,950 presentations, the analysis excludes invalid intervals, blank targets, and a three-event guard around temporal-group boundaries and interruptions, leaving 4,559 target events. The four test groups contain 667, 1,463, 1,466 and 963 events. Blank events remain in the history predictors. Each outer fit trains on three temporal groups, with ridge penalties and history-feature scaling selected within those training groups.

Models predict spike counts from current image identity (C), current and preceding image identity (H1), the current and three preceding identities (H3), or current and next identity (F1). Historical gains may be zero, so every augmented model contains C. The reported quantity is the reduction in standardized test SSE relative to C. Counts are extracted in nonoverlapping 25-ms bins. The minimum observed interval to the next stimulus is 250.188498 ms, placing all 0–250-ms bins before its onset.

| Post-stimulus window | H1 reduction in SSE relative to C | F1 reduction in SSE relative to C |
|---|---:|---:|
| 0–25 ms | 1.1156% | −.0024% |
| 25–50 ms | 1.3945% | −.0031% |
| 50–75 ms | 1.0083% | .0000% |
| 75–100 ms | .1644% | −.0029% |

H1 improves prediction in each of the four outer folds for the first three windows. These small increments establish predictable history information in this recording after conditioning on current identity. Residual sensory response, adaptation and other history-dependent processes remain possible contributors. The four folds characterize one animal's record. Full window and fold results are in `historical_allen_timebins_summary.csv` and `historical_allen_timebins_fold_improvements.csv`.

### E10.4 Transferring a measured response kernel to rapid sequences

The FASTIMAGES continuous-MEG example uses sub-01, 102 magnetometers, 100-Hz sampling and eight alternating localizer/sequence blocks. The localizer contributes 160 ordinary stimulus presentations; 40 distractor presentations are excluded from kernel estimation. The rapid-sequence data contain 64 five-image sequences, with nominal stimulus-onset intervals of 132, 164, 228 and 612 ms. Actual onset times enter each prediction.

Leave-one-block-out fitting estimates a separate 0–1-s response for every stimulus class and sensor from the other seven blocks. Training localizer residuals determine channel scaling. Each sequence uses its own prestimulus baseline. One candidate truncates a stimulus response when the next physical stimulus begins; the other sums the complete estimated responses to all stimuli already presented. Both share position-dependent nonnegative gains and a training-estimated common sensor offset.

| Stimulus-onset interval | Full response superposition: held-out $R^2$ | Response truncated at next onset: held-out $R^2$ |
|---|---:|---:|
| 132 ms | .031142 | .014922 |
| 164 ms | .013198 | .020630 |
| 228 ms | .026549 | .041725 |
| 612 ms | .003114 | .009622 |
| All sequences | .015091 | .019519 |

The richer response-overlap model improves prediction at the shortest interval and loses at the others. Adding nonnegative forward, backward or symmetric activation templates to the full sensory-history baseline changes pooled SSE by −.00295%, −.03332% and −.01346%, respectively. These templates use the physical order of each five-image sequence and a training-selected discount grid from .1 to .9. They calibrate a specified model comparison in these measurements; the physical sequence does not supply an independently established learned cognitive graph.

For fixed discounts and unit predictive strength, projecting out the stated history nuisance retains 88.94–97.13% of the forward/backward template distance, averaged across the eight folds. This geometric survival is conditional on that design and fixed parameter comparison. It does not establish exclusion of the complete continuous families or power at the measured noise level. The response kernel itself also requires validation when transferred from a localizer. Source tables are `historical_fastimages_kernel_summary.csv` and `historical_fastimages_kernel_discrimination_geometry.csv`.

### E10.5 Exact cyclic example: finite propagation coordinates are not physical lag weights

A separate algebraic example explains a condition needed before applying temporal moment tests. For true positive mixtures of discounts, the physical lag weights obey a moment inequality such as $w_1w_3-w_2^2\geq0$. Coordinates estimated in a finite propagation basis need not equal these $w_h$, because longer lags can return to the same propagation shapes.

Let $P$ be an $n$-state cycle with $n\geq5$, and let $A=CPC^\top$ act on its complete centered $(n-1)$-dimensional contrast space. Then

$$
A^n=I,
\qquad A^{n-1}=-\sum_{k=0}^{n-2}A^k.
$$

Grouping the infinite geometric series by its remainder modulo $n$ and then eliminating $A^{n-1}$ gives

$$
(I-\gamma A)^{-1}
=\sum_{k=0}^{n-2}b_kA^k,
\qquad
b_k=\frac{\gamma^k-\gamma^{n-1}}{1-\gamma^n}.
$$

The finite basis is linearly independent, so with reference $B=I$ all its coefficients are uniquely identified. Nevertheless,

$$
b_1b_3-b_2^2
=-\frac{\gamma^n(1-\gamma)^2}{(1-\gamma^n)^2}<0,
\qquad 0<\gamma<1.
$$

For $n=5$ and $\gamma=1/2$, the value is exactly $-8/961$. The residualized design formed by $AB,A^2B,A^3B$ after removing the current direction has rank three. Thus unique estimation of all three finite coordinates coexists with a negative determinant. The true infinite-kernel coefficients are $w_1=\gamma$, $w_2=\gamma^2$, $w_3=\gamma^3$ and satisfy $w_1w_3-w_2^2=0$. The negative finite-coordinate determinant is produced entirely by folding the tail into the basis.

To interpret fitted coefficients as actual lag weights, the measurement design must separate those lags from the tail. A nilpotent transient operator with $N^4=0$, together with an identifiable finite design, makes all higher propagation powers vanish and permits direct interpretation of its three positive-lag coefficients. Independent task manipulations or an explicit identifiable tail model provide other routes. The general cyclic example applies equally to the 12- and 40-state cycles used in calibration; the resolvent fits in E9 model the full series and do not apply the invalid finite-coordinate moment test. The complete archived derivation is `reference_archive/cyclic_lag_scope_proof.md`.

## E11. Research coverage and source provenance

The research programme examined complementary measurement situations. The table records what was available and the role each branch plays in the present argument. A branch with insufficient usable material remains an acquisition-status entry. Sample sizes across branches are not pooled: the two continuous-fMRI Ekman participants belong to its 35-person cohort, FASTIMAGES sub-01 belongs to the 30-person MEG cohort, and Huang's EEG and fMRI records concern the same 33 participants.

| Branch | Available material and analysis unit | Contribution to the argument | Primary source |
|---|---|---|---|
| D01, Tarder-Stoll | Four localizer scans and two anticipation scans from one participant; task and archive-field audit | Qualification of task-specific state matrices; the available event/path mapping was insufficient for path-specific operator reconstruction | [Article](https://doi.org/10.1038/s41467-024-53293-3); OpenNeuro ds005125 |
| D02, Ekman | Author-derived arrays for 35 participants; continuous fMRI from two of them | Post-learning reference content, cue-held-out comparison, reference prediction and physical-history calibration; E7 and E10.1–E10.2 | [Article](https://doi.org/10.7554/eLife.78904); [data](https://doi.org/10.34973/bsy6-9h29) |
| D03, Kern | DSMR101 MEG localizer and target measurements | Reference-target alias construction, continuous-family calibration and frozen-readout transfer; E6 and E9 | [Article](https://doi.org/10.7554/eLife.93357); [data](https://doi.org/10.5281/zenodo.8001755) |
| D04, He | EEG decoding/RDM/RSA outputs for 26 analyzable participants; 30-person supplied set for sensitivity | Learning-related content and direct directional contrasts at the available probability/geometry level | [Article](https://doi.org/10.1371/journal.pbio.3003740); [OSF](https://osf.io/sn9kv/) |
| D05, THINGS-fMRI | One participant, 100 images, 12 sessions and 72 preselected voxels | Full-space versus low-dimensional prediction; known-transform calibration, disjoint-session threshold evaluation | [THINGS-data](https://github.com/ViCCo-Group/THINGS-data); OpenNeuro ds004192 |
| D06, THINGS-MEG | One participant, one run in each of two sessions; five common images with one retained presentation each | Cross-measurement correspondence and qualification of available repeats | [THINGS-data](https://github.com/ViCCo-Group/THINGS-data) |
| D07, CNeuroMod THINGS | One participant, 36 sessions, 11,646 retained trials and seven ROIs | Across-image and across-session channel prediction under the supplied voxel-selection procedure | [Official code](https://github.com/courtois-neuromod/cneuromod-things); [dataset documentation](https://docs.cneuromod.ca/datasets/things.html) |
| D08, Turner | ID13, two EEG localizer recordings, 1,980 eligible events each | Temporal reference prediction, task-compatible readout, block restrictions and transfer; E2–E6 and E9 | [Article](https://doi.org/10.1371/journal.pbio.3003189) |
| D09, FASTIMAGES-MEG | Probability outputs for 30 participants; continuous MEG for sub-01 | Physical sequence timing, measured response-kernel transfer and history-template comparison; E10.4 | [Official data and code](https://github.com/skjerns/FASTIMAGES) |
| D10, FASTIMAGES-fMRI | Probability outputs for 40 participants, 3,000 trials and 13 TRs per trial | Temporal calibration at the probability-output level | [Official data and code](https://github.com/skjerns/FASTIMAGES) |
| D11, unavailable trial-level MEG branch | No qualifying public trial-level MEG/localizer files in the audited material | Data-availability boundary; no empirical test | Source status preserved in the research report |
| D12, Kahn | 322 participants and 535 complete graph-learning sequences | Behavioral readout competition, recency and local-trace alternatives, cross-person and cross-graph prediction | [Article](https://doi.org/10.1371/journal.pcbi.1013696); [OSF](https://osf.io/jupbz/) |
| D13, Garvert | ROI pair-response data for 23 participants and processed reaction times for 26 | Stationary time reversal on an undirected nonuniform graph, candidate geometry and run-held-out comparisons | [Article](https://doi.org/10.7554/eLife.17086) |
| D14, Allen Neuropixels | One mouse and one session; 5,950 natural-scene presentations; 110 VISp units for fine time bins | Observed sensory history, temporal measurement windows and rank/noise calibration; E10.3 | Allen Visual Coding Neuropixels; file `sub-707296975_ses-721123822.nwb` |
| D15, inaccessible source-archive branch | Original Dryad archive unavailable in the audit and no equivalent input obtained | Data-availability boundary; no empirical test | Source status preserved in the research report |
| D16, Huang | EEG outputs from 33 participants and fMRI outputs from the same 33; 2,561 fMRI trials | Probability-normalization constraints, trial-boundary verification and sequence statistics | [Article](https://doi.org/10.1038/s41467-024-51582-5); [data](https://zenodo.org/records/12547774) |

The formal first-round report, `SR_real_data_full_report_20261123.pdf`, and its claim-to-evidence table document all 16 branches. The second-round report and appendix preserve the reference, temporal-history, design-rank and known-transform analyses. The continuous-family release, `SR_dFB_complete_20261124.zip`, supplies the fold-level dFB records, independent reconstructions and continuous optimization checks used in E9. These reports and the compact source tables are preserved under `reference_archive/` and `source_data/`; their archive paths and hashes appear in the empirical provenance manifests. The current event-array calibration is separately documented in `source_data/README_fig07.md`.

The catalog has two functions. It retains the complete research coverage, including analyses that stopped at measurement qualification, and makes each empirical claim traceable to its own data and statistical unit. The original publications establish the identity and design of the datasets. Numerical results introduced here belong to the archived reanalyses described above.

### Empirical source references

[E1] Turner W, Sexton C, Johnson PA, Wilson EM, Hogendoorn H. Predictable motion is progressively extrapolated across temporally distinct processing stages in the human visual cortex. *PLOS Biology* 23, e3003189 (2025). https://doi.org/10.1371/journal.pbio.3003189

[E2] Kern S. Reactivation strength during cued recall is modulated by graph distance within cognitive maps. Public data archive, version 1 (2023). https://doi.org/10.5281/zenodo.8001755

[E3] Ekman M, Kusch S, de Lange FP. Successor-like representation guides the prediction of future events in human visual cortex and hippocampus. *eLife* 12, e78904 (2023). https://doi.org/10.7554/eLife.78904

[E4] Kern S, Nagel J, Gerchen MF, Gürsoy Ç, Meyer-Lindenberg A, Kirsch P, Dolan RJ, Gais S, Feld GB. Reactivation strength during cued recall is modulated by graph distance within cognitive maps. *eLife* 12, RP93357 (2024). https://doi.org/10.7554/eLife.93357

# C. Relation to prior work

## C1. Predictive maps, temporal context and direction

Dayan introduced the successor representation as a discounted occupancy prediction that separates environmental transition structure from reward [1]. Stachenfeld and colleagues developed its implications for hippocampal population codes, place fields and representational similarity [2]. The present analysis begins with that mechanistic interpretation and asks which parts survive an observation procedure.

Temporal direction has also been studied through alternative learning rules and behaviour. Gershman and colleagues related SR learning to temporal-context models and analysed conditions producing common and distinct predictions [3]. Keck and colleagues studied how symmetry in local learning rules changes predictive representations and generalization [10]. Momennejad and colleagues used revaluation to distinguish planning strategies [4], while Sharp and Eldar manipulated conditions under which forward and backward predictions support different decisions [9]. These studies establish both the relevance of competing mechanisms and the value of experimental intervention. The reciprocal-task results here add sufficient conditions for a direction property when a neural feature mapping is unknown, including the competing possibilities that the mechanism switches between conditions. The distinction between mixing maps and mixing transitions is made explicit by T22.

## C2. Neural observations and reference measurements

Relational geometry, future-related content and temporal scale are related empirical targets. Garvert and colleagues investigated abstract relational knowledge in the hippocampal–entorhinal system [5]. Brunec and Momennejad examined predictive representations across hippocampal and prefrontal hierarchies [6]. Ekman and colleagues compared future and past content using sequence learning, omission probes and localizer-based measurements [7]. Tarder-Stoll and colleagues studied hierarchical representations of past and future during multistep anticipation [8], and He and colleagues related sleep to successor-category evidence following learning [19].

The present results connect those observation levels to specific mathematical targets: a direction shared by every compatible explanation, a net direction in a nonnegative mixture, a discount range, or a finite sequence of temporal weights. A learned reference contributes its own predictive structure to the operator relating measurements. Temporal persistence in measurement supplies a further mechanism for apparent backward structure; T23 states the corresponding mean relation and the role of independent probes.

The empirical sources retain their original scientific aims. The Turner EEG study examined motion extrapolation across processing stages [18], and the Kern MEG study examined graph distance during cued recall [25]. Repeated measurements from these sources provide the reference estimates and known-transformation calibrations analysed here.

## C3. Representational models and statistical identification

Diedrichsen and Kriegeskorte placed encoding, pattern-component and representational-similarity methods in a common framework [12]. Gardner and Liu, and Sprague and colleagues, explained how encoding bases and model assumptions shape the interpretation of reconstructed responses [13,14]. Kriegeskorte and colleagues analysed circular selection [15], and Cai and colleagues examined task-related bias in representational similarity and a Bayesian treatment of that bias [16].

Here, task-compatible compression preserves a specified action of the task rather than requiring recovery of every neural feature. The compatible set is defined over complete labelled model families, so identical numerical predictions can retain different mechanistic labels. Property identification then asks whether all surviving explanations agree. Hotelling’s multivariate mean distribution [26] supplies the finite-sample residual threshold under the sampling conditions in T26–T27. This statistical step converts a population relation into a conditional error-controlled declaration; the actual sampling unit and covariance remain part of the empirical design.

## C4. Temporal recovery and the broader modelling programme

Momennejad and Howard used known multiscale SRs to recover a coarse future timeline through an approximate inverse Laplace transform [11]. The quantitative results here address unknown discounts and weights. Independent propagation shapes or distinct visible rotation blocks constrain shared parameters, while a finite-memory task based on a binary de Bruijn construction [27] separates the observable lags into different matrix positions. Direct and learned reverse references then permit different forms of kernel recovery.

The Discussion relates these results to the general question of what a successful fit establishes [17,20]. Its timing example refers to the models of Youguo Chen and colleagues [21,22]. Continuous parameter compensation is demonstrated analytically in the SR family studied here. The proposed use of Lie symmetries connects the broader programme to established symmetry-based work on identifiability and parameterization [23,24]. The mathematical results in this paper concern the stated SR families; the symmetry programme supplies a direction for constructing and comparing further model classes.

# D. References

1. Dayan, P. Improving Generalization for Temporal Difference Learning: The Successor Representation. *Neural Computation* **5**, 613–624 (1993). https://doi.org/10.1162/neco.1993.5.4.613

2. Stachenfeld, K. L., Botvinick, M. M. & Gershman, S. J. The hippocampus as a predictive map. *Nature Neuroscience* **20**, 1643–1653 (2017). Corrected version; correction: 10.1038/s41593-018-0133-1. https://doi.org/10.1038/nn.4650

3. Gershman, S. J., Moore, C. D., Todd, M. T., Norman, K. A. & Sederberg, P. B. The Successor Representation and Temporal Context. *Neural Computation* **24**, 1553–1568 (2012). https://doi.org/10.1162/NECO_a_00282

4. Momennejad, I., Russek, E. M., Cheong, J. H., Botvinick, M. M., Daw, N. D. & Gershman, S. J. The successor representation in human reinforcement learning. *Nature Human Behaviour* **1**, 680–692 (2017). https://doi.org/10.1038/s41562-017-0180-8

5. Garvert, M. M., Dolan, R. J. & Behrens, T. E. J. A map of abstract relational knowledge in the human hippocampal–entorhinal cortex. *eLife* **6**, e17086 (2017). https://doi.org/10.7554/eLife.17086

6. Brunec, I. K. & Momennejad, I. Predictive representations in hippocampal and prefrontal hierarchies. *Journal of Neuroscience* **42**, 299–312 (2022). https://doi.org/10.1523/JNEUROSCI.1327-21.2021

7. Ekman, M., Kusch, S. & de Lange, F. P. Successor-like representation guides the prediction of future events in human visual cortex and hippocampus. *eLife* **12**, e78904 (2023). https://doi.org/10.7554/eLife.78904

8. Tarder-Stoll, H., Baldassano, C. & Aly, M. The brain hierarchically represents the past and future during multistep anticipation. *Nature Communications* **15**, 9094 (2024). https://doi.org/10.1038/s41467-024-53293-3

9. Sharp, P. B. & Eldar, E. Humans adaptively deploy forward and backward prediction. *Nature Human Behaviour* **8**, 1726–1737 (2024). Publisher correction: 10.1038/s41562-024-01978-6. https://doi.org/10.1038/s41562-024-01930-8

10. Keck, J., Barry, C., Doeller, C. F. & Jost, J. Impact of symmetry in local learning rules on predictive neural representations and generalization in spatial navigation. *PLOS Computational Biology* **21**, e1013056 (2025). https://doi.org/10.1371/journal.pcbi.1013056

11. Momennejad, I. & Howard, M. W. Predicting the Future with Multi-scale Successor Representations. *bioRxiv* (2018). Preprint, version 1. https://doi.org/10.1101/449470

12. Diedrichsen, J. & Kriegeskorte, N. Representational models: A common framework for understanding encoding, pattern-component, and representational-similarity analysis. *PLOS Computational Biology* **13**, e1005508 (2017). https://doi.org/10.1371/journal.pcbi.1005508

13. Gardner, J. L. & Liu, T. Inverted Encoding Models Reconstruct an Arbitrary Model Response, Not the Stimulus. *eNeuro* **6**, ENEURO.0363-18.2019 (2019). https://doi.org/10.1523/ENEURO.0363-18.2019

14. Sprague, T. C., Boynton, G. M. & Serences, J. T. The Importance of Considering Model Choices When Interpreting Results in Computational Neuroimaging. *eNeuro* **6**, ENEURO.0196-19.2019 (2019). https://doi.org/10.1523/ENEURO.0196-19.2019

15. Kriegeskorte, N., Simmons, W. K., Bellgowan, P. S. F. & Baker, C. I. Circular analysis in systems neuroscience—the dangers of double dipping. *Nature Neuroscience* **12**, 535–540 (2009). https://doi.org/10.1038/nn.2303

16. Cai, M. B., Schuck, N. W., Pillow, J. W. & Niv, Y. Representational structure or task structure? Bias in neural representational similarity analysis and a Bayesian method for reducing bias. *PLOS Computational Biology* **15**, e1006299 (2019). https://doi.org/10.1371/journal.pcbi.1006299

17. Wilson, R. C. & Collins, A. G. E. Ten simple rules for the computational modeling of behavioral data. *eLife* **8**, e49547 (2019). https://doi.org/10.7554/eLife.49547

18. Turner, W., Sexton, C., Johnson, P. A., Wilson, E. M. & Hogendoorn, H. Predictable motion is progressively extrapolated across temporally distinct processing stages in the human visual cortex. *PLOS Biology* **23**, e3003189 (2025). https://doi.org/10.1371/journal.pbio.3003189

19. He, X., Büchel, P. K., Faghel-Soubeyrand, S., Klingspohr, J., Kehl, M. S. & Staresina, B. P. Sleep strengthens successor representations of learned sequences in humans. *PLOS Biology* **24**, e3003740 (2026). https://doi.org/10.1371/journal.pbio.3003740

20. Roberts, S. & Pashler, H. How persuasive is a good fit? A comment on theory testing. *Psychological Review* **107**, 358–367 (2000). https://doi.org/10.1037/0033-295X.107.2.358

21. Chen, Y., Peng, C. & Avitt, A. A unifying Bayesian framework accounting for spatiotemporal interferences with a deceleration tendency. *Vision Research* **187**, 66–74 (2021). https://doi.org/10.1016/j.visres.2021.06.005

22. Chen, Y., Yu, J., Hou, C. & Peng, C. Direct evidence for logarithmic magnitude representation in the central nervous system. *Neuroscience* **561**, 127–138 (2024). https://doi.org/10.1016/j.neuroscience.2024.10.025

23. Merkt, B., Timmer, J. & Kaschek, D. Higher-order Lie symmetries in identifiability and predictability analysis of dynamic models. *Physical Review E* **92**, 012920 (2015). https://doi.org/10.1103/PhysRevE.92.012920

24. Popovych, R. O. & Bihlo, A. Symmetry preserving parameterization schemes. *Journal of Mathematical Physics* **53**, 073102 (2012). https://doi.org/10.1063/1.4734344

25. Kern, S., Nagel, J., Gerchen, M. F., Gürsoy, Ç., Meyer-Lindenberg, A., Kirsch, P., Dolan, R. J., Gais, S. & Feld, G. B. Reactivation strength during cued recall is modulated by graph distance within cognitive maps. *eLife* **12**, RP93357 (2024). https://doi.org/10.7554/eLife.93357

26. Hotelling, H. The generalization of Student’s ratio. *Annals of Mathematical Statistics* **2**, 360–378 (1931). https://doi.org/10.1214/aoms/1177732979

27. de Bruijn, N. G. A combinatorial problem. *Proceedings of the Section of Sciences of the Koninklijke Nederlandse Akademie van Wetenschappen te Amsterdam* **49**, 758–764 (1946). https://research.tue.nl/en/publications/a-combinatorial-problem/
