### 4. What determines the scope of a mechanistic guarantee?

Direction now has a precise inferential meaning: within a specified task, reference relation and candidate family, do the explanations that remain all have a forward or backward property? This question can be answered because the experimental relation and the model constraints jointly rule out some possibilities. Understanding how those constraints work also explains the scope of the guarantee. We next change the task, the weights and the recording process in turn, and examine how the conclusion changes with them.

#### 4.1 Direction needs a shared task coordinate system

“Forward” is always relative to a particular relation between states. If the participant's internal transition structure is entirely free, a forward mechanism using one transition can make the same predictions as a backward mechanism using its time reversal. Formally, if the allowed internal transitions are closed under time reversal, the forward and backward families have identical observation sets (Supplementary Proposition T05). What is missing in this case is an anchor for direction, beyond a precise estimate of the neural patterns.

Experimental transition manipulations, sequence-memory tests and behavioural choices can provide independent evidence for that anchor. They constrain which state relations the participant uses. Cross-condition comparisons also require a shared neural mapping: the same state components must enter the recordings through a common set of underlying relations. Fixed electrodes and stimulus sets contribute experimental control; the common mapping is supported by these controls together with reference measurements. Both the task and the mapping therefore have explicit empirical roles in the model.

#### 4.2 Freedom in the current response changes the reciprocal guarantee

Complete SR binds the current and predictive terms by one rule. The reciprocal guarantee in the previous section, and several extensions with an additional nonnegative current component, use that structure. Allowing the current weight and predictive gain to vary independently expands the candidates to

$$
G_\sigma=cI+\alpha N_\gamma(\sigma J),
\qquad N_\gamma(A)=A(I-\gamma A)^{-1}.
$$

Here, $\sigma=+1$ and $-1$ denote forward and backward action relative to the task.

With two learned conditions, this expansion can restore an exact direction-reversal alias. Take $\gamma=1/2$ and $\alpha=5$. A forward explanation with current weight $c=1$ produces $-I+4J$ and $-I-4J$ in the two reciprocal conditions. A backward explanation with $c=3$ produces $I-4J$ and $I+4J$. If the second explanation uses the negative of the first explanation's common mapping, all plane patterns in both conditions remain unchanged (Supplementary Counterexample T07). Each explanation has an internally shared mapping, although the explanations assign opposite mechanism directions.

This result has a straightforward relation to the directly referenced half-planes. A direct reference anchors the underlying mapping. Comparing two learned conditions still allows different explanations to use different common mappings. The available information has changed. Complete SR with a nonnegative additional current term retains stronger constraints; its reciprocal sign guarantee extends to two-dimensional blocks, normal blocks of any dimension and square-zero blocks (Supplementary Proposition T09). How the current component enters a model should therefore be specified alongside its discount domain.

#### 4.3 Information about the past can originate in the recording process

A backward pattern has another possible source. Suppose the neural system encodes only the current state at each moment, but the recorded response retains activity from previous stimuli. One response has not disappeared when the next stimulus arrives. Conditional on the current state, previous states are distributed according to the time reversal $P^*$. The mean recorded pattern can consequently take the form

$$
\mathbb E[Y_t\mid S_t]
=\sum_{\ell\geq0}h_\ell(P^*)^\ell B,
$$

where $h_\ell$ weights the response from $\ell$ steps ago. Geometrically decaying weights produce a backward-SR form (Supplementary Theorem T23). The same matrix expression now describes the recording's history rather than an internal computation about the past.

This result separates two notions of time: propagation through a task is indexed by state transitions, whereas measurement persistence depends on elapsed physical time. Varying presentation intervals, estimating an independent response kernel and arranging independent probe trials can provide information that distinguishes them. For a stable, absolutely summable linear history kernel, independent and identically distributed probes make the average contribution of past states independent of current identity. Removing the common state component then removes this source of directional structure. This requires actual independence between past and current identities; a randomly ordered sequence does not establish it by itself. Reciprocal tasks also reverse the conditional relations in measurement history, so history control and reciprocal manipulation serve distinct purposes.

#### 4.4 A property may be identified while other model claims remain testable

On a single plane, every target vector $y$ can be written uniquely as $y=(uI+vJ)b$ for any nonzero reference $b$. When the current term and predictive gain are free, the union of the forward and backward families can describe every such pair of means. The sign of $v$ retains its directional meaning within that family, but the pair of means leaves no residual constraint on the fit of the expanded model as a whole (Supplementary Theorem T06).

Two scientific questions therefore need to be distinguished: identifying a property within a model family, and testing whether that family accounts for the data. Additional independent templates can test whether the same action holds across readouts. Task blocks with different propagation properties can test whether one parameter set explains several relations simultaneously. The latter also offer a route out of the single-scale/multiscale alias: operators that are exactly equal on one plane will remain equal as repeated measurement of that plane becomes more precise.

Structural identification also makes demands on precision. A weak reference, a predictive component close to zero or nearly dependent propagation shapes can place competitors close together. More effective independent repetitions can reduce uncertainty, but a decision still depends on the distance from the present observation to the competing family. Intersection elsewhere in two families does not prevent this observation from lying in an identifiable region. Conversely, distinct ideal directions need not be separated well enough in every recording. Continuous families and joint uncertainty make these two questions explicit (Supplementary Propositions T25–T28).

### 5. From directional properties to discounts and temporal weights

The preceding results leave a constructive question. If an observation supports a direction, what more is needed to determine how far it predicts, or how it combines different future distances? The answer depends on the resolution of the desired explanation. Some observations first constrain a parameter range. Other tasks make different discounts produce independent shapes. A finite-memory design makes individual lags visible. We consider these possibilities in turn.

#### 5.1 Directional regions already contain partial quantitative information

For two pure complete forward SRs on a reciprocal quarter-turn plane,

$$
\phi=\arctan\gamma+\arctan\eta.
$$

Suppose the relative angle is exactly $70^\circ$. Because both discounts are less than one, each contribution is less than $45^\circ$. To sum to $70^\circ$, each must also exceed $25^\circ$. Hence

$$
\gamma>\tan25^\circ\simeq0.4663,
\qquad \eta>\tan25^\circ\simeq0.4663.
$$

This reasoning requires no choice of a best-fitting point. The compatible discounts still form a continuous curve, but every point on it excludes short horizons (Figure 5a). The bounds are sharp: as one discount approaches one, the other can approach its lower bound. Supplementary Theorem T15 gives the full parameter set for a general angle and a common discount ceiling.

For measured data, the joint compatibility set can be projected onto each discount axis to report the values that remain possible. If the candidates are multiscale mixtures, the same angular constraints concern effective phase ratios. Constraining each constituent horizon requires observations that distinguish kernel shape. Partial identification therefore states both what is known and where further information is needed.

#### 5.2 Three propagation shapes determine when parameters are unique

On a single plane, changing the discount can be compensated by the current weight and predictive gain. To prevent that compensation, a task must make different propagation steps leave independent shapes in the recordings. For a specified direction branch, a known task $P$ and a population direct reference $B$, consider

$$
Y=cB+\alpha P(I-\gamma P)^{-1}B,
\qquad c\in\mathbb R,\quad\alpha>0,\quad0\leq\gamma<1.
$$

$B$ is the direct shape, $PB$ is its one-step propagation and $P^2B$ is its two-step propagation. We prove that $c,\alpha,\gamma$ are globally structurally identifiable if and only if these three matrices are linearly independent (Supplementary Theorem T17). Independence is assessed by vectorizing each matrix and comparing the resulting three shapes.

Why three? Comparing the observations from two parameter sets and clearing the two inverse matrices gives a linear combination of $B$, $PB$ and $P^2B$. If these shapes are independent, each coefficient must vanish separately, which forces the parameters to agree. If the shapes are dependent, propagation occupies too small a space, and changes in the current term and predictive gain can absorb a change in discount. The complete proof establishes both the necessary and the sufficient direction.

This condition gives experimental design a concrete target. Additional channels may expose more propagation shapes, but channel count is not the condition itself. A reference with few channels that preserves independent propagation shapes can carry more parameter information than a high-dimensional reference whose shapes repeat. Parameter recovery depends on those shapes and their stability.

#### 5.3 Two different planes provide a shared constraint

Two task planes with different rotation properties give an explicit realization. On the quarter-turn plane, many parameter sets can produce $I+J$. Requiring the same parameters to explain a second plane, with a different rotation angle, adds a new demand: the original compensation must now preserve the second observation as well.

Let the two task blocks be $A_j=a_jI+\beta_jJ$, with $a_j^2+\beta_j^2=1$, $\beta_j>0$ and $a_1\neq a_2$. When $c,\alpha,\gamma$ are shared, their rotational coefficients are

$$
v_j=\frac{\alpha\beta_j}{1-2a_j\gamma+\gamma^2}.
$$

The ratio of $v_1/\beta_1$ to $v_2/\beta_2$ removes the common predictive gain. The remaining ratio uniquely determines the discount in its allowed domain. We then recover $\alpha$ and obtain $c$ from the non-rotational component of each plane (Supplementary Theorem T18; Figure 5b). The $90^\circ$ and $45^\circ$ Fourier planes of an eight-state cycle realize this structure.

The second plane supplies another opportunity: the current weight calculated independently from the two blocks must agree. This is an additional prediction that can fail. More information therefore sharpens parameter estimation and places the model under further test. Shared parameters are central to this design. If each plane can select all its own parameters, the cross-plane constraint disappears.

<!--FIG05-->

#### 5.4 A finite-memory task separates future lags

If the question concerns the shape of the temporal weights themselves, we can arrange a task that separates the lags directly. Let each state be a four-bit word, such as $0101$. At each step, remove the first bit and append a random bit, producing either $1010$ or $1011$ with equal probability. After four steps, every original bit has been replaced.

The sixteen words form a binary de Bruijn transition graph [27]. Take four prespecified contrasts, each comparing states in which one bit is zero with states in which it is one. The task acts on these contrasts as a lower-shift matrix $N$, with $N_{i+1,i}=1$ and all other entries zero. It satisfies $N^4=0$: after four steps, the initial information has disappeared from this contrast space, while the stochastic task continues to evolve (Figure 6a; Supplementary Proposition T19).

We call the sequence of weights across prediction steps a temporal kernel. On a cycle, a later lag can return to the same propagation shape after another circuit, so coefficients in a finite basis may combine contributions from the more distant tail. The finite-memory block separates this relation explicitly: after its memory length, the initial information disappears from the selected contrasts. Any temporal kernel on the observable block becomes

$$
f(N)=w_0I+w_1N+w_2N^2+w_3N^3.
$$

Current, one-step, two-step and three-step contributions now occupy different matrix positions. One temporal weight can no longer be changed and concealed solely by adjusting another while preserving the entire observable action. A qualified direct reference and visibility of the corresponding operator directions allow these coefficients to be read separately.

#### 5.5 Temporal kernels directly test single-scale and multiscale relations

A single complete SR applies one geometric rule across lags. Its first three positive-lag weights must satisfy

$$
w_2^2=w_1w_3.
$$

This equality removes overall scale and leaves the current term aside. For same-direction nonnegative mixtures whose scales share one neural mapping, $w_h=\int\gamma^h\mu(d\gamma)$ instead gives

$$
w_2^2\leq w_1w_3.
$$

When $w_1>0$, the difference has a direct interpretation:

$$
\frac{w_1w_3-w_2^2}{w_1^2}
=\operatorname{Var}_{\nu}(\gamma),
\qquad\nu(d\gamma)=\frac{\gamma\mu(d\gamma)}{w_1}.
$$

It is the variance of discounts reweighted by their contribution at positive lags. Greater dispersion of the contributing predictive horizons increases this quantity (Figure 6c; Supplementary Theorem T20). Relative to the effective discount available on one rotation plane, the finite-memory block provides a direct test of temporal-weight shape.

The task's memory length determines the observable range. A four-bit task exposes three positive lags; a longer binary-memory task extends that range. A zero-discount component contributes no positive-lag weight and remains invisible to these quantities. Inferences about the infinite tail require additional kernel structure or more observable lags.

#### 5.6 A learned reference can still permit recovery of the target kernel

The preceding kernel readout used a direct reference. We now return to the difficulty in Section 2: can the target weights be recovered when the reference has itself been learned? In the finite-memory design, the triangular structure of a reverse reference provides a route.

Let the target be $f(N)$ and the unknown learned reference be $q(N^\top)$, with current coefficient $q_0\neq0$. Their relative action is

$$
T=af(N)q(N^\top)^{-1}.
$$

Because $N^\top$ shifts upwards, it annihilates the first basis vector. The reference inverse therefore only rescales this vector by $1/q_0$. The first column of the relative operator satisfies

$$
T_{h+1,1}=\frac{af_h}{q_0}.
$$

The reference changes overall scale but preserves the target's relative lag weights in this column. Whenever these entries are observable, the relative shape of the target's finite kernel can be recovered; a directional sign interpretation uses $a/q_0>0$. For two pure complete SRs, additional matrix relations separately determine the target and reference discounts (Figure 6d; Supplementary Theorem T21).

A learned reference thus has two consequences. In some tasks, it makes opposite mechanisms produce the same relative action. In others, known structure allows target and reference to leave separately readable traces. The decisive issue is which relations the observation preserves.

<!--FIG06-->

#### 5.7 Zero net direction still leaves mechanistic questions to answer

A zero directional statistic also invites finer questions. One mechanism forms forward and backward maps at the same discount and then mixes them equally. Another first symmetrizes the transition relation and forms an SR from that transition at a positive discount. The relation between symmetric learning and symmetrized transitions has been examined in earlier modelling work [10]. Both constructions have zero net direction, but they perform different computations.

On a fourth-order finite-memory block, their diagonal predictions differ. Mixing the two maps leaves a uniform main diagonal. Propagating through a symmetrized transition includes return paths and creates a different diagonal structure. A common current term and a positive overall gain preserve this distinction (Supplementary Proposition T22). A directional statistic without a positive or negative answer therefore need not exhaust the mechanistic information available in richer observations.

These examples give a route to increasing explanatory resolution: report the properties shared by compatible explanations, determine the parameter ranges that remain, and use independent propagation information to recover point parameters or finite temporal kernels. This is a progression of scientific questions. Statistically, properties and parameter ranges can be derived together from one joint compatibility set. When the scientific target requires the complete operator, a full-rank reference and valid perturbation bounds provide a corresponding recovery route (Methods; Supplementary Theorem T24).
