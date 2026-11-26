#!/usr/bin/env python3
"""Portable deterministic algebra checks for the SR manuscript.

Run from any directory:
    python /path/to/SR_manuscript/code/math_checks.py

The only third-party dependency is SymPy. Results are written to
``source_data/math_checks.json`` relative to this source package.

These checks reproduce specified identities and representative exact examples
from Supplementary Results T04, T08, T11 and T19--T24. They do not replace the
general proofs, generate neural data, simulate participants, or estimate power.
T07, T12, T17, T18 and partial identification are additionally checked by the
package's existing ``code/math_additions.py``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
Q = sp.Rational
ASSERTION_COUNT = 0


def require(condition: Any, label: str) -> None:
    """Fail even when Python is run with -O; record completed assertions."""
    global ASSERTION_COUNT
    if not bool(condition):
        raise AssertionError(label)
    ASSERTION_COUNT += 1


def reduced(expr: Any) -> sp.Expr:
    return sp.cancel(sp.simplify(expr))


def require_zero(expr: Any, label: str) -> None:
    if isinstance(expr, sp.MatrixBase):
        require(all(reduced(value) == 0 for value in expr), label)
    else:
        require(reduced(expr) == 0, label)


def inner(left: sp.MatrixBase, right: sp.MatrixBase) -> sp.Expr:
    require(left.shape == right.shape, "Frobenius inner-product dimensions")
    return reduced(sum(left[i, j] * right[i, j]
                       for i in range(left.rows) for j in range(left.cols)))


def frobenius(matrix: sp.MatrixBase) -> sp.Expr:
    return sp.sqrt(inner(matrix, matrix))


def sr(discount: Any, matrix: sp.MatrixBase) -> sp.MatrixBase:
    return (sp.eye(matrix.rows) - discount * matrix).inv()


def kernel(coefficients: list[Any], matrix: sp.MatrixBase) -> sp.MatrixBase:
    result = sp.zeros(matrix.rows)
    for lag, coefficient in enumerate(coefficients):
        result += coefficient * matrix**lag
    return result


def strings(matrix: sp.MatrixBase) -> list[list[str]]:
    return [[str(reduced(matrix[i, j])) for j in range(matrix.cols)]
            for i in range(matrix.rows)]


def fraction_and_decimal(expr: Any) -> dict[str, Any]:
    expr = reduced(expr)
    return {"exact": str(expr), "decimal": float(sp.N(expr, 18))}


def require_doubly_stochastic(matrix: sp.MatrixBase, label: str) -> None:
    one = sp.ones(matrix.rows, 1)
    require(matrix.rows == matrix.cols, f"{label}: square")
    require(all(value >= 0 for value in matrix), f"{label}: nonnegative")
    require_zero(matrix * one - one, f"{label}: row sums")
    require_zero(matrix.T * one - one, f"{label}: column sums")


def check_paley() -> dict[str, Any]:
    n = 7
    transition = sp.zeros(n)
    for state in range(n):
        for lag in (1, 2, 4):
            transition[state, (state + lag) % n] = Q(1, 3)
    require_doubly_stochastic(transition, "Paley")
    h = sp.eye(n) - sp.ones(n) / n
    k = (transition - transition.T) / 2
    require(h.rank() == 6, "Paley: contrast dimension")
    require_zero(h * transition * h + h / 6 - k, "Paley: contrast block")
    require_zero(k * k + Q(7, 36) * h, "Paley: skew-block square")

    gamma, eta, beta, gain = Q(3, 10), Q(4, 5), Q(75, 173), Q(634, 519)
    relative_full = sr(gamma, transition) * (sp.eye(n) - eta * transition)
    backward_full = gain * sr(beta, transition.T)
    relative = h * relative_full * h
    backward = h * backward_full * h
    require_zero(relative - backward, "Paley: complete contrast alias")
    require_zero(relative - (Q(53, 48) * h - Q(25, 56) * k),
                 "Paley: exact contrast coefficients")
    require(relative_full != backward_full, "Paley: constant direction excluded")

    return {
        "theorem": "T04",
        "passed": True,
        "arithmetic": "exact rational",
        "states": n,
        "contrast_dimension": 6,
        "transition": strings(transition),
        "target_discount": str(gamma),
        "reference_discount": str(eta),
        "equivalent_backward_discount": str(beta),
        "equivalent_positive_gain": str(gain),
        "relative_contrast_coefficients_H_K": ["53/48", "-25/56"],
        "scope": "Equality holds after projection by H on both sides. The constant-state direction is deliberately not included.",
    }


def check_eight_cycle() -> dict[str, Any]:
    n = 8
    transition = sp.zeros(n)
    for state in range(n):
        transition[state, (state + 1) % n] = 1
    require_doubly_stochastic(transition, "Eight-cycle")
    j = sp.Matrix([[0, -1], [1, 0]])
    c = sp.Matrix([[sp.cos(sp.pi * state / 2) / 2 for state in range(n)],
                   [-sp.sin(sp.pi * state / 2) / 2 for state in range(n)]])
    other = transition**3
    require_zero(c * sp.ones(n, 1), "Eight-cycle: contrasts")
    require_zero(c * c.T - sp.eye(2), "Eight-cycle: orthonormal rows")
    require_zero(c * transition - j * c, "Eight-cycle: CP=JC")
    require_zero(c * other + j * c, "Eight-cycle: CP^3=-JC")
    require_zero(c * transition.T + j * c, "Eight-cycle: transpose block")
    require_zero(c * other.T - j * c, "Eight-cycle: other transpose block")
    require(other != transition.T, "Eight-cycle: whole graphs are not reversals")
    w = (transition - transition.T - other + other.T) / 2
    require_zero(w - c.T * (2 * j) * c, "Eight-cycle: rank-two marker")
    require(w.rank() == 2, "Eight-cycle: marker rank")
    require_zero(inner(w, w) - 8, "Eight-cycle: marker norm")

    gamma, eta, gain = Q(3, 10), Q(4, 5), Q(7, 6)
    full_relative = gain * sr(gamma, transition) * (sp.eye(n) - eta * other)
    plane_relative = gain * sr(gamma, j) * (sp.eye(2) + eta * j)
    require_zero(c * full_relative - plane_relative * c,
                 "Eight-cycle: relative action commutes with restriction")
    score = inner(full_relative, w) / 8
    expected = gain * (gamma + eta) / (2 * (1 + gamma**2))
    require_zero(score - expected, "Eight-cycle: normalized marker formula")
    require(score > 0, "Eight-cycle: forward score")
    return {
        "theorem": "T11",
        "passed": True,
        "arithmetic": "exact rational",
        "contrast_matrix": strings(c),
        "marker_rank": 2,
        "marker_squared_frobenius_norm": "8",
        "discounts": [str(gamma), str(eta)],
        "gain": str(gain),
        "relative_plane_operator": strings(plane_relative),
        "normalized_forward_score": str(reduced(score)),
        "scope": "P and P^3 are reciprocal on the prescribed reducing plane, although P^3 is not P^T.",
    }


def check_reciprocal_signs() -> dict[str, Any]:
    j = sp.Matrix([[0, -1], [1, 0]])
    cycle = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    swap = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    transition = sp.eye(3) / 2 + cycle / 3 + swap / 6
    require_doubly_stochastic(transition, "Nonnormal reciprocal task")
    u = sp.Matrix([[1 / sp.sqrt(2), 1 / sp.sqrt(6)],
                   [-1 / sp.sqrt(2), 1 / sp.sqrt(6)],
                   [0, -2 / sp.sqrt(6)]])
    require_zero(u.T * u - sp.eye(2), "Nonnormal task: contrast coordinates")
    require_zero(u.T * sp.ones(3, 1), "Nonnormal task: constant removed")
    task_block = sp.simplify(u.T * transition * u)
    # Its parent is doubly stochastic, so its Euclidean contrast action is a
    # contraction. Verify this instance independently by exact positive minors.
    contraction_gap = sp.eye(2) - task_block.T * task_block
    require(contraction_gap[0, 0] > 0 and contraction_gap.det() > 0,
            "Nonnormal task: exact contraction certificate")
    require(task_block * task_block.T != task_block.T * task_block,
            "Nonnormal task: not normal")

    three = sp.Matrix([[Q(1, 5), Q(1, 3), 0],
                       [0, -Q(1, 4), Q(1, 6)],
                       [0, 0, Q(1, 7)]])
    require(inner(three, three) < 1,
            "Three-dimensional example: Frobenius contraction certificate")
    require(three * three.T != three.T * three,
            "Three-dimensional example: not normal")

    examples = [
        ("quarter_turn", j, "A^T A=I"),
        ("nonnormal_stochastic_task_contrast", task_block,
         "I-A^T A has positive exact leading principal minors"),
        ("three_dimensional_nonnormal_contraction", three,
         "squared Frobenius norm is strictly below one"),
    ]
    rows = []
    gamma, eta = Q(1, 3), Q(4, 5)
    for name, a, certificate in examples:
        identity = sp.eye(a.rows)
        s, k = (a + a.T) / 2, (a - a.T) / 2
        require(k != sp.zeros(a.rows), f"{name}: nonzero direction")
        f = sr(gamma, a)
        forward = f * (identity - eta * a.T)
        backward = sr(gamma, a.T) * (identity - eta * a)
        symmetric = sr(gamma, s) * (identity - eta * s)
        skew_formula = f * ((gamma + eta) * k
                            - gamma * eta * (s * k + k * s)) * f.T
        require_zero((forward - forward.T) / 2 - skew_formula,
                     f"{name}: exact skew identity")
        ff, bb, symmetric_score = inner(forward, k), inner(backward, k), inner(symmetric, k)
        require(ff > 0, f"{name}: FF trace sign")
        require(bb < 0, f"{name}: BB trace sign")
        require_zero(symmetric_score, f"{name}: specified symmetric control")
        # The theorem keeps strict sign with one zero discount, but not with two.
        one_zero = identity - eta * a.T
        require(inner(one_zero, k) > 0, f"{name}: one-zero-discount boundary")
        require_zero(inner(identity, k), f"{name}: both-zero-discount boundary")
        rows.append({
            "name": name,
            "matrix": strings(a),
            "contraction_certificate": certificate,
            "discounts": [str(gamma), str(eta)],
            "FF_score": str(ff),
            "BB_score": str(bb),
            "symmetric_score": str(symmetric_score),
            "one_zero_discount_score": str(inner(one_zero, k)),
            "both_zero_discount_score": "0",
        })
    return {
        "theorem": "T08",
        "passed": True,
        "arithmetic": "exact rational or algebraic",
        "examples": rows,
        "scope": "Representative exact instances, including nonnormal blocks; the general sign proof remains in T08. The arbitrary three-dimensional block is an algebraic example, not an experimental or participant dataset.",
    }


def make_debruijn(k: int = 4) -> tuple[sp.MatrixBase, sp.MatrixBase, sp.MatrixBase]:
    n = 2**k
    transition = sp.zeros(n)
    for state in range(n):
        for bit in (0, 1):
            transition[state, ((state << 1) & (n - 1)) | bit] += Q(1, 2)
    c = sp.Matrix([[(-1 if (state >> (k - 1 - bit)) & 1 else 1) / sp.sqrt(n)
                    for state in range(n)] for bit in range(k)])
    shift = sp.zeros(k)
    for index in range(k - 1):
        shift[index + 1, index] = 1
    return transition, c, shift


def check_debruijn(transition: sp.MatrixBase, c: sp.MatrixBase,
                   shift: sp.MatrixBase) -> dict[str, Any]:
    k, n = c.shape
    require_doubly_stochastic(transition, "de Bruijn")
    require_zero(c * sp.ones(n, 1), "de Bruijn: zero-mean contrasts")
    require_zero(c * c.T - sp.eye(k), "de Bruijn: orthonormal contrasts")
    require_zero(c * transition - shift * c, "de Bruijn: lower-shift action")
    require_zero(c * transition.T - shift.T * c, "de Bruijn: reverse action")
    require_zero(shift**k, "de Bruijn: nilpotent contrast block")
    require(shift**(k - 1) != sp.zeros(k), "de Bruijn: exact nilpotency index")
    require_zero(transition**k - sp.ones(n) / n,
                 "de Bruijn: full task remains stochastic after memory replacement")
    gamma = Q(3, 5)
    require_zero(c * sr(gamma, transition) - sr(gamma, shift) * c,
                 "de Bruijn: resolvent restricted to contrast block")
    return {
        "theorem": "T19",
        "passed": True,
        "arithmetic": "exact rational",
        "memory_bits": k,
        "states": n,
        "contrast_matrix": strings(c),
        "contrast_shift": strings(shift),
        "nilpotency_index": k,
        "full_task_k_step_action": "ones(n,n)/n, not zero",
        "resolved_discount_example": str(gamma),
        "scope": "A lawful doubly stochastic finite-memory task; only the prescribed contrast block is nilpotent.",
    }


def psi(matrix: sp.MatrixBase, lag: int) -> sp.Expr:
    k = matrix.rows
    return reduced(sum(matrix[index + lag, index] - matrix[index, index + lag]
                       for index in range(k - lag)) / (k - lag))


def check_direct_kernel(shift: sp.MatrixBase) -> dict[str, Any]:
    k = shift.rows
    coefficients = list(sp.symbols(f"w0:{k}", real=True))
    gain = sp.symbols("a", positive=True)
    current = sp.symbols("c", real=True)
    forward = gain * kernel(coefficients, shift) + current * sp.eye(k)
    backward = gain * kernel(coefficients, shift.T) + current * sp.eye(k)
    for lag in range(1, k):
        require_zero(psi(forward, lag) - gain * coefficients[lag],
                     f"Direct kernel: positive lag {lag}")
        require_zero(psi(backward, lag) + gain * coefficients[lag],
                     f"Direct kernel: reverse lag {lag}")
        require_zero(psi((forward + backward) / 2, lag),
                     f"Direct kernel: symmetric lag {lag}")

    gamma = sp.symbols("gamma", real=True)
    geometric = gain * sr(gamma, shift) + current * sp.eye(k)
    require_zero(psi(geometric, 2)**2 - psi(geometric, 1) * psi(geometric, 3),
                 "Direct kernel: geometric moment equality")

    discounts, weights = [Q(1, 4), Q(3, 4)], [Q(2, 5), Q(3, 5)]
    moments = [sum(weight * discount**lag for weight, discount in zip(weights, discounts))
               for lag in range(k)]
    gap = reduced(moments[1] * moments[3] - moments[2]**2)
    require(gap > 0, "Direct kernel: positive nondegenerate mixture gap")
    nu = [weight * discount / moments[1] for weight, discount in zip(weights, discounts)]
    require_zero(sum(nu) - 1, "Direct kernel: lag-reweighted probability")
    nu_mean = sum(weight * discount for weight, discount in zip(nu, discounts))
    variance = sum(weight * (discount - nu_mean)**2
                   for weight, discount in zip(nu, discounts))
    require_zero(gap / moments[1]**2 - variance,
                 "Direct kernel: reweighted variance identity")
    scaled = Q(7, 5) * kernel(moments, shift) + Q(11, 10) * sp.eye(k)
    recovered_variance = (psi(scaled, 1) * psi(scaled, 3) - psi(scaled, 2)**2) / psi(scaled, 1)**2
    require_zero(recovered_variance - variance, "Direct kernel: gain and current invariance")

    zero_mass_moments = [Q(2, 5) * sp.Integer(0)**lag + Q(3, 5) * Q(1, 2)**lag
                         for lag in range(k)]
    require_zero(zero_mass_moments[1] * zero_mass_moments[3] - zero_mass_moments[2]**2,
                 "Direct kernel: zero-discount mass remains invisible")
    return {
        "theorem": "T20",
        "passed": True,
        "arithmetic": "symbolic identities and exact rational mixture",
        "positive_lag_functionals": [str(psi(forward, lag)) for lag in range(1, k)],
        "geometric_gap": "0",
        "mixture": {
            "discounts": [str(value) for value in discounts],
            "weights": [str(value) for value in weights],
            "moments_w0_to_w3": [str(value) for value in moments],
            "w1_w3_minus_w2_squared": str(gap),
            "lag_reweighted_weights": [str(value) for value in nu],
            "lag_reweighted_discount_variance": str(reduced(variance)),
        },
        "zero_discount_mass_example_gap": "0",
        "scope": "A direct reference and observable operator functionals are assumed. Nonnegative scales share one mapping. Equality does not exclude zero-discount mass or identify an infinite tail.",
    }


def check_learned_reverse_kernel(shift: sp.MatrixBase) -> dict[str, Any]:
    k = shift.rows
    target_coefficients = list(sp.symbols(f"f0:{k}", real=True))
    q0 = sp.symbols("q0", real=True, nonzero=True)
    reference_coefficients = [q0, *sp.symbols(f"q1:{k}", real=True)]
    gain = sp.symbols("a", positive=True)
    target = kernel(target_coefficients, shift)
    reference = kernel(reference_coefficients, shift.T)
    relative = gain * target * reference.inv()
    first = sp.eye(k)[:, 0]
    require_zero(shift.T * first, "Learned reference: first-vector annihilation")
    require_zero(reference.inv() * first - first / q0,
                 "Learned reference: inverse reference on first vector")
    for lag in range(k):
        require_zero(relative[lag, 0] - gain * target_coefficients[lag] / q0,
                     f"Learned reference: first-column coefficient {lag}")

    gamma, eta = sp.symbols("gamma eta", real=True)
    pure = gain * sr(gamma, shift) * (sp.eye(k) - eta * shift.T)
    require_zero(pure[0, 0] - gain, "Learned pure SR: gain recovery")
    require_zero(pure[1, 0] / pure[0, 0] - gamma,
                 "Learned pure SR: target discount recovery")
    require_zero(-pure[0, 1] / pure[0, 0] - eta,
                 "Learned pure SR: reference discount recovery")
    require_zero(pure[2, 0]**2 - pure[1, 0] * pure[3, 0],
                 "Learned pure SR: first-column geometric relation")

    numeric_exact = pure.subs({gain: 1, gamma: Q(3, 5), eta: Q(4, 5)})
    wrong_probe_gap = reduced(psi(numeric_exact, 2)**2
                             - psi(numeric_exact, 1) * psi(numeric_exact, 3))
    require_zero(wrong_probe_gap - Q("-0.18607104"),
                 "Learned reference: direct-reference probe is not transferable")
    return {
        "theorem": "T21",
        "passed": True,
        "arithmetic": "symbolic identities and exact rational example",
        "first_column": [str(reduced(relative[lag, 0])) for lag in range(k)],
        "pure_SR_inverse": {"a": "T[0,0]", "gamma": "T[1,0]/T[0,0]", "eta": "-T[0,1]/T[0,0]"},
        "array_indices": "zero-based in this JSON; one-based in the manuscript",
        "wrong_direct_probe_example": {
            "target_discount": "3/5",
            "reference_discount": "4/5",
            "gain": "1",
            "psi2_squared_minus_psi1_psi3": fraction_and_decimal(wrong_probe_gap),
        },
        "scope": "q0 is nonzero; directional sign also requires a/q0>0. The relevant first-column entries must be observable. This is a different readout from the directly referenced subdiagonal differences.",
    }


def check_zero_direction_mechanisms(shift: sp.MatrixBase) -> dict[str, Any]:
    gamma, eta = sp.symbols("gamma eta", real=True)
    mixed_maps = (sr(gamma, shift) + sr(gamma, shift.T)) / 2
    symmetrized_transition = (shift + shift.T) / 2
    symmetric_sr = sr(eta, symmetrized_transition)
    require_zero(mixed_maps - mixed_maps.T, "Zero direction: map mixture symmetric")
    require_zero(symmetric_sr - symmetric_sr.T, "Zero direction: transition SR symmetric")
    for index in range(shift.rows):
        require_zero(mixed_maps[index, index] - 1,
                     f"Zero direction: map-mixture diagonal {index}")
    x = eta / 2
    expected = x**2 / (1 - 3 * x**2 + x**4)
    difference = reduced(symmetric_sr[1, 1] - symmetric_sr[0, 0])
    require_zero(difference - expected, "Zero direction: exact diagonal contrast")
    value = reduced(difference.subs(eta, Q(3, 5)))
    require(value > 0, "Zero direction: positive-discount exact example")
    scaled = Q(7, 5) * symmetric_sr + Q(11, 10) * sp.eye(shift.rows)
    require_zero(scaled[1, 1] - scaled[0, 0] - Q(7, 5) * difference,
                 "Zero direction: common current cancels and positive gain preserves contrast")
    return {
        "theorem": "T22",
        "passed": True,
        "arithmetic": "symbolic identity and exact rational positive-discount example",
        "map_mixture_diagonal": ["1"] * shift.rows,
        "symmetric_transition_diagonal_difference": str(difference),
        "eta_example": "3/5",
        "positive_diagonal_difference_example": fraction_and_decimal(value),
        "scope": "The equally mixed forward/backward maps use the same discount. A positive discount in the symmetrized-transition SR yields a different diagonal prediction; at zero discount the strict distinction disappears.",
    }


def check_operator_error_bound() -> dict[str, Any]:
    # Exact two-row, three-column reference with orthonormal rows. No random
    # samples or empirical observations are generated.
    observed_reference = sp.Matrix([[Q(3, 5), Q(4, 5), 0],
                                    [-Q(12, 25), Q(9, 25), Q(4, 5)]])
    require_zero(observed_reference * observed_reference.T - sp.eye(2),
                 "Operator bound: orthonormal reference rows")
    require(observed_reference.rank() == observed_reference.rows,
            "Operator bound: observed reference full row rank")
    true_operator = sp.diag(Q(3, 2), Q(1, 2))
    reference_error = sp.diag(Q(1, 20), 0) * observed_reference
    true_reference = observed_reference - reference_error
    true_target = true_operator * true_reference
    target_error = sp.diag(Q(1, 100), 0) * observed_reference
    observed_target = true_target + target_error
    pseudoinverse = observed_reference.T * (observed_reference * observed_reference.T).inv()
    observed_operator = observed_target * pseudoinverse
    e_b, e_y = frobenius(reference_error), frobenius(target_error)
    s = sp.Integer(1)
    require_zero(e_b - Q(1, 20), "Operator bound: actual reference perturbation")
    require_zero(e_y - Q(1, 100), "Operator bound: actual target perturbation")
    require(s > e_b, "Operator bound: singular-value margin")
    require_zero(observed_operator - sp.diag(Q(287, 200), Q(1, 2)),
                 "Operator bound: exact recovered operator")
    operator_norm = Q(287, 200)  # Largest singular value of the verified diagonal matrix.
    actual_error = frobenius(observed_operator - true_operator)
    bound = reduced((e_y + operator_norm * e_b) / (s - e_b))
    require(actual_error <= bound, "Operator bound: deterministic T24 inequality")
    w = sp.diag(1, 0)
    functional_error = abs(inner(observed_operator - true_operator, w))
    require(functional_error <= bound * frobenius(w),
            "Operator bound: prescribed functional error")
    return {
        "theorem": "T24",
        "passed": True,
        "arithmetic": "exact rational, including norms in this constructed example",
        "observed_reference_shape": [2, 3],
        "observed_reference_rank": 2,
        "observed_reference": strings(observed_reference),
        "true_operator": strings(true_operator),
        "recovered_operator": strings(observed_operator),
        "e_B": str(e_b),
        "e_Y": str(e_y),
        "reference_minimum_singular_value": str(s),
        "recovered_operator_spectral_norm": str(operator_norm),
        "actual_frobenius_error": fraction_and_decimal(actual_error),
        "T24_upper_bound": fraction_and_decimal(bound),
        "functional_error": str(functional_error),
        "scope": "The population relation Y=T B, full row rank, and actual perturbation bounds hold by construction. This is a deterministic inequality check, not a statistical confidence interval or a neural simulation.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "source_data" / "math_checks.json",
                        help="JSON output path (default: package source_data/math_checks.json)")
    args = parser.parse_args()
    transition, contrasts, shift = make_debruijn(4)
    checks = {
        "paley_complete_contrast_alias": check_paley(),
        "eight_cycle_common_reducing_plane": check_eight_cycle(),
        "reciprocal_trace_sign_examples": check_reciprocal_signs(),
        "debruijn_nilpotent_contrast_block": check_debruijn(transition, contrasts, shift),
        "direct_finite_kernel_and_positive_mixture": check_direct_kernel(shift),
        "learned_reverse_reference_kernel": check_learned_reverse_kernel(shift),
        "zero_direction_competing_mechanisms": check_zero_direction_mechanisms(shift),
        "complete_operator_deterministic_error": check_operator_error_bound(),
    }
    result = {
        "schema_version": 1,
        "manuscript_revision": "2026-11-26",
        "generator": "code/math_checks.py",
        "dependency_versions": {"sympy": sp.__version__},
        "scope": "Deterministic symbolic/rational checks of stated manuscript algebra and specified examples. No neural data, simulated participants, random draws, or power calculations.",
        "proof_status": "Representative computational checks accompany, and do not replace, the general proofs in the theoretical supplement.",
        "complementary_checks": "code/math_additions.py",
        "checks": checks,
        "check_groups": len(checks),
        "completed_assertions": ASSERTION_COUNT,
        "all_checks_passed": all(check["passed"] for check in checks.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"all_checks_passed": result["all_checks_passed"],
                      "check_groups": len(checks),
                      "completed_assertions": ASSERTION_COUNT,
                      "output": str(args.output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
