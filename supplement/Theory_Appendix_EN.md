---
title: "Theory and statistical supplement: from overlapping SR families to property and parameter identification"
lang: en
date: "2026-11-25"
---

# Theory and statistical supplement

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
