### 1. Two mechanistic questions about predictive maps

#### 1.1 From state transitions to future occupancy

On a familiar journey from home to the office, we pass a junction and a station, perhaps stopping at a shop. At the junction, we know where we are. We also know which places are likely to occur next in our experience. The successor representation (SR) incorporates this second kind of knowledge into the representation of the current state. It asks: “Starting from here, what might I encounter in the future?”

Dayan introduced the SR in reinforcement learning in 1993, using relationships between successor states to organize representations for long-term reward prediction [1]. The model can be understood by beginning with one transition. Let $P_{ij}$ be the probability of reaching state $j$ on the next step when the current state is $i$. Each row of the transition matrix $P$ describes the possible destinations from one state. Suppose that, from the junction, the probability of going next to the station is $0.75$ and that of going to the shop is $0.25$. From the station and the shop, the probabilities of reaching the office on the following step are $0.8$ and $0.4$, respectively. The probability of reaching the office from the junction in two steps is the sum over these two routes: $0.75\times0.8+0.25\times0.4=0.70$.

Matrix multiplication performs this calculation for all pairs of states:

$$
(P^2)_{ik}=\sum_jP_{ij}P_{jk}.
$$

Thus, $P$ gives the distribution one step ahead, $P^2$ the distribution two steps ahead, and $P^3$ the distribution three steps ahead. A matrix power represents successive transitions: each multiplication sums over all possible intermediate states. For a stable task or behavioural policy, one-step transitions can therefore be unfolded into predictions of the future.

To combine these predictions into a representation of the present, the SR assigns weight $\gamma^t$ to step $t$. We consider continuing tasks with finitely many states and take $0\leq\gamma<1$, giving

$$
M_\gamma(P)=I+\gamma P+\gamma^2P^2+\gamma^3P^3+\cdots.
$$

Here, $I$ represents the current state, $\gamma P$ the weighted distribution one step ahead, and $\gamma^2P^2$ the weighted distribution two steps ahead. A larger $\gamma$ gives more weight to distant events. Each row of the matrix describes discounted occupancy of all states, starting from one current state. A state can be visited again at several future steps, so the row sums occupancy across time and its total can exceed one.

The series also has a compact expression. Multiplying it by $\gamma P$ produces the same sequence of terms beginning with $\gamma P$. Subtracting this expression from the original series cancels every later term, leaving $I$. Hence,

$$
(I-\gamma P)M_\gamma=I,
\qquad
M_\gamma(P)=(I-\gamma P)^{-1}.
$$

The inverse matrix is therefore a concise way of writing the accumulation of occupancy across steps. It describes a relationship among representations that a learning system can approach through successive updates. If immediate rewards at the states form a vector $r$, then $V=M_\gamma r$ converts predicted occupancy into long-term value. Predictions of states and the rewards assigned to those states can thereby be organized separately.

Stachenfeld, Botvinick and Gershman developed this structure into the predictive-map account of the hippocampus [2]. On this account, population activity at the current state contains information about states that may be visited. Two locations can consequently have similar population representations when they predict similar futures, even if the locations are physically far apart. The SR thus provides a concrete hypothesis linking learning, cognitive maps and neural representation.

#### 1.2 The direction and weighting of a cognitive map

Once a representation contains information beyond the current state, its temporal direction becomes a separate question. At the junction, information about home can arise from the route just travelled, whereas information about the station can concern the route ahead. A mechanism can integrate forwards, integrate backwards, or combine both. Research on temporal context, predictive hierarchies and multistep anticipation gives each possibility a substantive cognitive interpretation [3,6,8,9,11].

To define direction, a forward model uses transitions $P$ from the present to the future. A backward model uses the distribution of preceding states conditional on the present. If the task has stationary state distribution $\pi$, its backward transition matrix is

$$
P^*=D_\pi^{-1}P^\top D_\pi,
$$

where the diagonal entries of $D_\pi$ are the stationary state probabilities. The weighting reflects a simple fact: states that occur more often can also be more frequent sources of a given present state. When stationary occupancy is uniform, the expression reduces to $P^*=P^\top$. We first use such tasks to develop the intuition; the observation coordinates and metrics for nonuniform occupancy are given in Methods and the Supplementary Information.

Temporal weighting presents a second question. Even when observations support information about the future, they may leave open how that information is combined. The SR uses the geometric sequence $1,\gamma,\gamma^2,\ldots$, linking the current response, one-step prediction and more distant prediction. This is a substantive model constraint. The sensory response to the current state can vary independently, the overall gain of predictive information can change, and a population can combine several predictive horizons.

A slightly larger family makes this distinction explicit. Define

$$
N_\gamma(P)=P+\gamma P^2+\gamma^2P^3+\cdots
=P(I-\gamma P)^{-1},
$$

and consider $G(P)=cI+\alpha N_\gamma(P)$. The coefficient $c$ controls the current component, $\alpha$ the overall strength of future propagation, and $\gamma$ the relative weighting of successive future steps. A scaled complete SR, $aM_\gamma(P)$, satisfies $c=a$ and $\alpha=a\gamma$. Here $N_0(P)=P$ still contains one-step prediction, whereas the complete SR at zero discount, $M_0(P)=I$, contains only the current component. The complete SR therefore occupies a constrained part of this larger family. We will see that releasing this coupling changes what the observations can determine.

We focus on these two commitments of a predictive representation: the temporal direction of its information and the weighting of different temporal distances. They admit answers at different levels of precision. Data may determine direction while remaining compatible with several horizons, or may constrain direction together with a range of parameter values. We begin with these separately testable questions and build towards more specific mechanistic interpretations.

#### 1.3 How neural measurements carry these relationships

The mathematical model uses states as coordinates, whereas experiments usually record neurons, voxels or electrodes. An electrode does not arrive labelled as a measurement of the station, the office, or a combination of several sources. To connect the two descriptions, write the state-averaged patterns as

$$
Y=G(P)B+E.
$$

Rows of $Y$ correspond to states and columns to recording channels. The matrix $B$ maps state components into those channels, and $E$ represents measurement error. The mechanism $G(P)$ specifies how information about the present, past and future propagates among states. Across experimental conditions, a common $B$ means that the same underlying neural mapping is maintained. Positive overall gains, state-independent offsets and nonuniform state weights are treated together in Methods (Supplementary Result T01).

The empirical question is now more specific. Two mechanisms can differ in their internal operation; do they also leave different relationships after measurement? We begin with a commonly used summary of representation. Even before measurement error is introduced, a mechanistic difference can disappear at this step.

### 2. Two exact ways in which different mechanisms produce the same observation

#### 2.1 From activity patterns to representational distances

Suppose that the activity patterns recorded at locations A and B are

$$
x_A=(1,0.7,0.4),\qquad x_B=(1,0.6,0.5).
$$

The patterns are similar. At a third location, C, the pattern is

$$
x_C=(0.1,0.2,1.5),
$$

which differs more markedly from A. We can then ask whether A and B are similar because they predict similar futures, whereas C predicts a different future.

Representational similarity analysis makes such questions quantitatively accessible [12]. For example, the squared Euclidean distance between two patterns is

$$
d(A,B)^2=\sum_j(x_{Aj}-x_{Bj})^2.
$$

In this illustrative example, $d(A,B)^2=0.02$ and $d(A,C)^2=2.27$. Computing distances between every pair of states gives a representational dissimilarity matrix (RDM). A model predicts which states should be close and which should be separated, and these predictions can be compared with the distance structure of neural activity.

This operation compresses complete activity vectors into distances between states. Rotating all the vectors together, or reflecting them, preserves every pairwise distance. The same RDM can therefore arise from different coordinate relationships. Whether the direction of a predictive map falls within this invariance requires its own analysis.

We obtain an exact answer in canonical coordinates whose features correspond to task states. Let $P$ be doubly stochastic, so that every row and every column sums to one, and let the forward and backward SRs use the same nonzero discount $0<\gamma<1$. Their complete Euclidean row distances satisfy

$$
\operatorname{RDM}\{M_\gamma(P)\}
=\operatorname{RDM}\{M_\gamma(P^\top)\}
\quad\Longleftrightarrow\quad
PP^\top=P^\top P.
$$

The condition on the right is called normality. It permits strongly directed transitions. Moving repeatedly around a directed cycle and moving around it in reverse are different trajectories, yet these transitions can satisfy the condition. Their forward and backward predictive maps then have identical distances between states (Fig. 1; Supplementary Result T02).

Loss of direction in geometry therefore does not require a nearly directionless behavioural process. On a directed cycle, the probability of taking a reverse edge can approach zero, making the process increasingly irreversible while preserving forward–backward geometric equivalence at every admissible discount. Several independent rotation blocks also allow multiple labelled tasks with different directions and identical geometry (Supplementary Result T03). The full proof gives upper and lower bounds on geometric differences in terms of $PP^\top-P^\top P$.

This result characterizes the identification boundary of canonical Euclidean observations. Other feature metrics, neural mappings and partial observations can be analysed through their corresponding observation maps. The starting point is that a task with a clear direction need not preserve that direction in its pairwise distances. When similarity agrees with a predictive map, the relevant competitors include every mechanism that predicts the same similarity.

#### 2.2 Why a reference matters when complete patterns are retained

A natural next step is to retain complete patterns and examine relationships that distances discard. Complete patterns, however, still pass through an unknown neural mapping. In $Y=MB$, a single observation contains both $M$ and $B$; a change in mechanism can be accompanied by a change in the mapping paired with it. A reference measurement provides a common basis for this decomposition.

The most direct case is a reference pattern $B$ and a target pattern $MB$. They pass through the same recording channels, allowing us to ask how state information has propagated from the reference to the target. When the reference contains the information needed to recover the relevant action, we can obtain a relative operator. An operator is the rule that transforms one pattern into another. It retains how coordinates combine and is consequently richer than pairwise distances.

The scientific meaning of the reference is therefore central. A direct reference can anchor the underlying state patterns. A reference that has already learned the relevant transitions carries its own predictive structure into the comparison. Both can yield clear and reproducible recordings, but they define different inferential objects.

To express this distinction, we first remove the mean component shared by all states. The remaining coordinates describe differences between states and are called state contrasts. Let $A$ denote the task action in these coordinates, with forward SRs $M_\gamma(A)$ and $M_\eta(A)$ in the target and reference conditions. If the conditions share a neural mapping and the relative action is observable, then

$$
T=aM_\gamma(A)M_\eta(A)^{-1},\qquad a>0,
$$

where $a$ is the ratio of overall gains. The recovered object describes the relationship between two learned representations.

#### 2.3 A learned reference can create a complete directional alias

Two forward SRs on the same task can have different predictive horizons. Rearranging their relative operator gives

$$
T=a\{I+(\gamma-\eta)A(I-\gamma A)^{-1}\}.
$$

Directionality is thereby linked to the relationship between the two horizons. A directional marker formed from the directed part of the task is positive for a forward SR measured against a direct reference. Against a forward reference learned on the same task, its sign becomes the sign of $\gamma-\eta$ (Supplementary Result T04). A target that predicts a shorter distance into the future can therefore produce a backward relative marker when compared with a reference that predicts farther ahead. Both underlying conditions are forward; the comparison expresses their difference.

The effect can extend beyond a sign reversal. We construct a seven-state task in which each state transitions with equal probability to states 1, 2 or 4 steps ahead in its cyclic numbering. Its entire six-dimensional state-contrast space admits the exact relation

$$
M_{0.3}(A)M_{0.8}(A)^{-1}
=aM_\beta(A^\top),
$$

with $a=634/519$ and $\beta=75/173$, both valid positive gain and discount parameters. The left side arises from two forward learning conditions with different horizons; the right side is a backward SR. Equality holds throughout the six-dimensional contrast operator, so every state-contrast prediction determined by that operator is the same (Fig. 1; Supplementary Result T04).

We call exact equality between different mechanisms in a specified observation an alias. The two examples reach it by different routes. Geometric aliasing arises when activity patterns are compressed into distances. Learned-reference aliasing persists after a relative operator has been retained, because the reference itself forms part of the mechanism being compared. The first directs attention to what the summary preserves; the second to what the baseline defines.

Exact aliasing also distinguishes two experimental difficulties. Uncertainty caused by noise can shrink as effective repeated measurement increases. Explanations that are identical at the observed level converge on the same result as precision improves. New mechanistic information then requires additional constraints among tasks, references or readouts. The question develops from how accurately a quantity can be estimated to which properties it already determines and how further properties can be made identifiable.

<!--FIG01-->

### 3. Identifying properties shared by all compatible explanations

#### 3.1 Continuous prediction families turn counterexamples into an identification problem

A pair of exact counterexamples shows that two explanations can be confused. We also want to know how much of the prediction space is ambiguous and whether there are regions in which a conclusion is possible. We therefore let discounts, gains and permitted reference parameters vary over their entire admissible ranges and consider every observation they can generate. A mechanistic label then corresponds to a continuous prediction set, extending beyond the point obtained from a particular fit.

Write a complete explanation as $\theta$, including its direction, numerical parameters and permitted measurement factors. For a given experiment and readout, $\Phi(\theta)$ is the observation that explanation produces. Given an observed value $y$, define

$$
\mathcal C(y)=\{\theta:\Phi(\theta)=y\}.
$$

This compatible set retains every candidate capable of explaining the same observation. If the property of interest is $q(\theta)$, such as forward versus backward direction, net direction or membership in a parameter range, the question becomes whether the explanations in $\mathcal C(y)$ agree about that property.

This yields a general identification principle: a property is determined by observation exactly when explanations producing the same observation share that property. Its global form is

$$
\Phi(\theta_1)=\Phi(\theta_2)
\Longrightarrow
q(\theta_1)=q(\theta_2).
$$

The reason is direct. If the same observation remains compatible with opposing properties, a rule that reads only that observation receives the same input under the two true mechanisms and cannot be guaranteed to classify both correctly. If every compatible explanation shares the property, reading it from any one explanation gives the same answer; that answer is determined by the observation. This is the standard principle of identifiability. Here we use it to derive the overlap and distinguishable regions of specific continuous SR families (Supplementary Result T05).

The geometric and reference aliases now have a common interpretation: the observation retains candidates with different directions. Conversely, a parameter may remain nonunique because the compatible set contains a curve along which its value changes, even though every explanation on that curve shares the same direction. We next show this situation explicitly. It gives qualitative properties a distinct status that can be established mathematically.

#### 3.2 Direction can be determined while discount varies continuously

Consider a task action that is easy to draw. Summarize the state patterns of a cyclic task along two predetermined contrast axes, such that one transition acts on the axes as a counterclockwise rotation through $90^\circ$. Write

$$
J=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad J^2=-I.
$$

The axes describe differences between states. They can be measured jointly through many neural channels; they do not require two neurons dedicated to rotation. We give a concrete eight-state realization below.

With a suitable direct reference, consider the family in which current and predictive contributions vary independently. The current coefficient in state-contrast coordinates is $c\in\mathbb R$, the predictive gain is $\alpha\geq0$, and the discount is $0\leq\gamma<1$. Forward action is

$$
G_+=cI+\alpha J(I-\gamma J)^{-1}=uI+vJ,
$$

where

$$
u=c-\frac{\alpha\gamma}{1+\gamma^2},\qquad
v=\frac{\alpha}{1+\gamma^2}.
$$

The coefficient $u$ gives the contribution along the original direction, and $v$ the contribution along the known task rotation. As all admissible parameters vary, the forward family covers exactly the half-plane $v\geq0$. The backward model for the same task uses the time-reversed action $-J$, covering exactly $v\leq0$ (Fig. 2; Supplementary Result T06). A nonzero observable $v$ therefore determines direction within this expanded family. The line $v=0$ is the current-only boundary shared by both families.

Discount has a different status. For every $0\leq\gamma<1$, choosing

$$
c=1+\gamma,\qquad \alpha=1+\gamma^2
$$

produces

$$
G_+=I+J.
$$

For example, $\gamma=0$ gives $c=1$ and $\alpha=1$, whereas $\gamma=0.5$ gives $c=1.5$ and $\alpha=1.25$. The models differ in their current contribution and future decay but produce exactly the same action. Varying discount continuously traces an entire curve of such equivalents.

The same observation can therefore establish a forward property while leaving discount undetermined along this curve. The distinction remains even when observation is exact. For observations within the admissible prediction set of the complete SR, coupling current and future weights can reduce the compatible set, sometimes to a unique parameter value. That precision depends jointly on the observations and the coupling rule.

This example makes the relationship between general and special models concrete. A numerical parameter can appear well determined within one model and lose uniqueness when other admissible weights are released. A property such as direction can nevertheless remain unchanged throughout the expanded set. Property identification thus supplies a scientific statement that can survive parameter variation.

#### 3.3 A multiscale representation can have the appearance of a single discount

The same plane yields another result. Suppose a common neural mapping carries several forward SRs with nonnegative weights, together with an optional nonnegative current component. Their combined action is

$$
H(J)=r_0I+\int M_\gamma(J)\,\mu(d\gamma),
\qquad r_0\geq0,\quad\mu\geq0.
$$

Every component can be expressed as $uI+vJ$ on this plane. Adding them still leaves only these two coefficients. Whenever the total action is nonzero, it can therefore be written exactly as

$$
H(J)=a_{\mathrm{eff}}M_{\gamma_{\mathrm{eff}}}(J),
\qquad \gamma_{\mathrm{eff}}=v/u.
$$

A representation with several genuine horizons can thus appear on this observation plane as a complete SR with one effective discount (Fig. 2; Supplementary Result T16). The equality concerns the entire action, so greater precision on the same plane continues to recover the same effective model.

The effective discount has a clear descriptive role: it summarizes the relative rotational component on this plane. For mechanistic interpretation, the compatible set retains different compositions of horizons. All can be forward while differing about how many scales exist and how large they are. Continuous-family analysis therefore provides a useful approach to multiscale SRs: establish the properties preserved across possible compositions, then seek additional observations that reveal differences between them.

<!--FIG02-->

#### 3.4 Reciprocal tasks add directional constraints between learned conditions

A direct reference supported the preceding half-plane result. We now return to a reference that has itself been learned. We keep the underlying neural mapping common across two conditions, change their transition structures, and ask what information this experimental constraint adds. We first consider complete SRs, then extend the analysis to nonnegative mixtures and specified additional current components.

Again use a quarter-turn plane. Because $J^2=-I$, a complete SR is

$$
M_\gamma(J)=\frac{I+\gamma J}{1+\gamma^2}.
$$

This is a positive scaling combined with a rotation through $\arctan\gamma$. Two forward SRs on the same task have relative angle $\arctan\gamma-\arctan\eta$, whose sign depends on their horizons.

Now let the target task act as $J$ on the common plane and the reference task as $-J$. Both conditions still predict forwards relative to their own tasks, with target discount $\gamma$ and reference discount $\eta$. The action from reference to target is

$$
T=aM_\gamma(J)M_\eta(-J)^{-1}.
$$

The target has a positive angle and the reference a negative angle. Taking their relative relationship therefore adds the angles:

$$
\phi=\arctan\gamma+\arctan\eta.
$$

As long as at least one discount is nonzero, the angle is positive. Two backward mechanisms give the opposite sign. The order of target and reference is fixed when interpreting the sign; reversing the comparison reverses the relative angle.

The directional guarantee also extends to general real contraction blocks of task contrasts. Denote the block by $A$, and let the two conditions implement $A$ and $A^\top$ in common coordinates while maintaining positive gains and a shared neural mapping. The antisymmetric part $K=(A-A^\top)/2$ captures the directed action of the task. When $K\neq0$ and at least one discount is nonzero, two complete forward SRs satisfy $\langle T,K\rangle_F>0$, whereas two backward SRs give a negative value. This result permits nonnormal tasks and different horizons in the two conditions (Supplementary Result T08). The scope of extensions with nonnegative additional current components is described in Results 4.

Implementation can be simplified further. Two complete transition graphs need not reverse everywhere; the state relationships used for analysis need to carry the reciprocal action. Consider eight states arranged in a cycle. The first task advances one state at a time, and the second advances three states. Their complete transition matrices are $P$ and $P^3$. Along a predetermined pair of contrasts, the phases of the eight states repeat four directions. A one-state advance then rotates the plane through $90^\circ$, whereas a three-state advance rotates it through $270^\circ$, or $90^\circ$ in reverse. If $C$ extracts those two contrasts, then

$$
CP=JC,\qquad CP^3=-JC.
$$

Both tasks preserve this analysis space, allowing its action to be described separately. This is a common reducing subspace (Fig. 3; Supplementary Result T11). It provides a way to obtain reciprocal information using the same stimuli organized into different sequences. Multisequence tasks already offer an experimental setting for studying context-dependent representations of the past and future [8]. Common-subspace analysis identifies explicit task relationships capable of supporting the present identification guarantee.

#### 3.5 Which observations remain exclusive when mechanisms can switch direction?

Changing the task can also change the mechanism participants use. We therefore expand the candidates to four combinations: both target and reference are forward (FF), both are backward (BB), the target is forward and the reference backward (FB), or the target is backward and the reference forward (BF). Forward and backward are defined relative to each condition's own task.

In the $J,-J$ tasks, let each condition use a complete SR, allow its discount to range over $[0,g]$ with $g<1$, and let positive overall gains vary freely. Set $h=\arctan g$. Calculating the relative angles generated by all parameters gives the exact ranges

| Mechanisms in the two conditions | Complete relative-angle range |
|---|---|
| FF | $[0,2h]$ |
| BB | $[-2h,0]$ |
| FB or BF | $[-h,h]$ |

The table translates each mechanistic label into its complete empirical predictions (Fig. 3; Supplementary Result T13). Over the full domain of discounts below one, $30^\circ$ can be produced by two forward angles, $15^\circ+15^\circ$, or by a directional switch, $40^\circ-10^\circ$. With admissible positive gains, the explanations can produce the same relative operator. A comparison restricted to FF and BB would omit the second compatible mechanism.

The expanded families also leave clear exclusive regions. Each complete SR rotates by less than $45^\circ$, so a combination of one forward and one backward mechanism has relative angle within $(-45^\circ,45^\circ)$. A positive angle strictly greater than $45^\circ$ excludes both switching explanations and BB; the negative outer region supports a common backward direction. With a shared upper bound $g<1$, the boundaries contract to $\pm h$. Requiring a strict crossing of the boundary also defines a clear geometric target for inference with finite error.

The role of continuous families is particularly visible here. Mechanisms need not be separated everywhere: the location of the current observation determines the available conclusion. Overlap retains several explanations, whereas an outer region supports a shared property. Deriving this distinction over complete parameter ranges makes it directly useful for experimental design and uncertainty analysis.

The result extends to a condition containing several forward and backward complete SRs, provided their weights and any additional current contribution are nonnegative and all components share the neural mapping. For a nonzero total action, the result remains $uI+vJ$, with $u>0$ and $|v|\leq gu$. The same angle ranges can therefore be organized by each condition's net direction. The positive outer region requires both conditions to be net forward; the negative outer region requires both to be net backward (Supplementary Result T14). This identifiable property permits components pointing in the opposite direction. It concerns their combined action. Multiscale composition remains in the compatible set, while the common directional constraint has become stronger.

<!--FIG03-->

#### 3.6 A task-compatible readout preserves the information needed for direction

These results identify useful relationships between actions. Measurement must next preserve them across many neural channels. Recovering an arbitrary complete operator requires the reference to span every state direction to be recovered. When the question concerns only a sign or angle, less information may suffice (Supplementary Result T10).

We first extract the common plane from the task structure, then select a template using training data and freeze it before formal comparison. The readout turns patterns on that plane into two numbers. Its defining requirement is compatibility with the task: applying a task action and then reading it out gives the same relationship as reading out first and then applying the corresponding two-dimensional action.

Writing the linear readout as $z$, compatibility is

$$
z\{(uI+vJ)X\}=(uI+vJ)z(X).
$$

One construction takes inner products of a pattern with a training template and with its $J$-rotated version. The templates supply two comparison axes separated by $90^\circ$, preserving the scaling and rotation present in the candidate action. The full construction and its commutation identity are given in Methods and Supplementary Result T12.

Let the readout of the population reference mean be $b$ and that of the target be $y$. If the original patterns share the relation $T=uI+vJ$, then the readouts still satisfy

$$
y=(uI+vJ)b.
$$

Their signed area is

$$
\det[b,y]=v\|b\|^2.
$$

Whenever the true reference projection $b$ is nonzero, the area preserves the sign of $v$. Together with the dot product $b^\top y=u\|b\|^2$, it also yields the relative angle. A pair of task-compatible coordinates can therefore preserve direction, while full reference rank belongs to the stronger goal of recovering an arbitrary operator (Fig. 4).

The readout selects information around the scientific question. The training template need not equal the true reference, but the population reference must have a nonzero projection onto the chosen comparison axes. Additional templates can improve visibility or test a shared action, with all templates retained within the same repeated unit. Mechanisms that are already equal throughout the task plane remain equal under additional readouts of that plane; separating them requires another propagation relationship.

Candidate comparison can also be restricted to the common prediction space actually used by the candidates. Directions unrelated to their differences need not consume the same error budget. A common task plane specifies the experimental action, a compatible readout preserves that action, and a common prediction space focuses comparisons on the relevant components. Each aligns the information requirement with the question being asked.

#### 3.7 Joint uncertainty brings compatible sets to finite measurements

So far, the relationships have been exact. With finite measurements, reference and target patterns both vary around their population means, and errors within the same repeated unit can be correlated. Directional inference needs to carry these variations into a common compatible set.

Let the formal observations be repeated pairs $(b_i,y_i)$. A candidate action $T_\theta$ predicts that the target mean equals the transformed reference mean. We can therefore examine the residual

$$
r_i(\theta)=y_i-T_\theta b_i.
$$

Under the candidate relationship, its population mean is zero. Residual covariance includes target error, reference error and their within-pair association. If the reference varies substantially, the candidate action multiplying it also changes the resulting uncertainty. Joint treatment retains this dependence.

Given independently trained and frozen tasks, templates and candidate spaces, we invert a Hotelling mean test to retain every compatible continuous candidate. The stated sampling conditions are independent and identically distributed Gaussian pairs, a positive-definite residual covariance, and a sample size exceeding residual dimension. Write the retained set as $\widehat{\mathcal C}_\alpha$. When the true mechanism belongs to the prespecified candidate space, its probability of retention is $1-\alpha$ (Methods; Supplementary Result T26).

A property is declared using the same principle as before: the retained set must be nonempty, and all its explanations must agree about the property. If the true candidate is retained, their shared property necessarily agrees with it. Consequently,

$$
\Pr\{\text{a property is declared and is wrong}\}\leq\alpha.
$$

The guarantee connects continuous model families directly to a scientific conclusion. It controls the probability that the procedure makes an erroneous declaration; error conditional on having made a declaration is a separate statistical target. Different directional labels are retained even when they generate the same operator, so an exact alias naturally prevents a unique directional declaration.

The analysis admits several informative outcomes. If every compatible candidate is forward, a forward property can be declared. If opposing directions remain, directional ambiguity is retained. If no candidate is compatible, the specified mechanism and observation relationship require examination. If the reference has no visible projection or the repetition structure is insufficient, the corresponding information gap is reported. Distinguishing these outcomes allows even an inconclusive experiment to identify what information is needed next.

Error control also does not require complete separation between the families everywhere. Their global minimum distance may be zero because some boundaries intersect, while the current observation remains far from every opposing property. Our separation results express such local judgments through error neighbourhoods and candidate distances (Supplementary Result T25). In application, windows, templates and task choices are fixed using independent training, and repeated units are defined by acquisition. Results 6 connects these requirements to actual measurements.

Continuous families thus specify the explanations to retain; experiments supply constraints across conditions; task-compatible readouts preserve the relevant action; and joint uncertainty determines which candidates finite data can exclude. Direction and horizon composition already have different identification status. We next examine how the guarantees change with model and measurement assumptions, before seeking information that can recover the temporal weights themselves.

<!--FIG04-->
