# Supplementary empirical methods and results: reference estimation and conditional directional calibration with real repeated measurements

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
