#!/usr/bin/env python3
"""Exact algebra added for the 2026-11-26 SR manuscript audit.

No neural data, virtual participants, or stochastic power simulations are used.
Run from any working directory; outputs go to ../source_data/math_additions.json.
"""
from pathlib import Path
import csv
import json
import math
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source_data"
OUT.mkdir(parents=True, exist_ok=True)


def zero(expr):
    return sp.cancel(sp.expand(expr)) == 0


def matrix_zero(matrix):
    return all(zero(value) for value in matrix)


results = {
    "scope": "Deterministic symbolic audit; no participant or neural observations",
    "date": "2026-11-26",
    "checks": {},
}

# Archived Paley alias, kept here as standalone exact figure source data.
P7 = sp.zeros(7)
for i in range(7):
    for lag in (1, 2, 4):
        P7[i, (i+lag) % 7] = sp.Rational(1, 3)
H7 = sp.eye(7)-sp.ones(7)/7
K7 = (P7-P7.T)/2
assert matrix_zero(H7*P7*H7-(-H7/6+K7))
assert matrix_zero(K7*K7+7*H7/36)
F = lambda gamma, matrix: (sp.eye(matrix.rows)-gamma*matrix).inv()
paley_relative = H7*F(sp.Rational(3,10),P7)*(sp.eye(7)-sp.Rational(4,5)*P7)*H7
paley_backward = sp.Rational(634,519)*H7*F(sp.Rational(75,173),P7.T)*H7
assert paley_relative == paley_backward
assert paley_relative == sp.Rational(53,48)*H7-sp.Rational(25,56)*K7
with (OUT/"paley_transition_matrix.csv").open("w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["state"]+[f"state_{j}" for j in range(7)])
    for i in range(7):
        writer.writerow([i]+[str(P7[i,j]) for j in range(7)])
results["checks"]["paley_six_dimensional_operator_alias"] = {
    "exact": True,
    "states": 7,
    "contrast_dimension": 6,
    "forward_target_discount": "3/10",
    "forward_reference_discount": "4/5",
    "apparent_backward_discount": "75/173",
    "apparent_positive_gain": "634/519",
    "restriction": "Both sides projected by H7; constant state direction is not aliased.",
    "matrix_source": "paley_transition_matrix.csv",
}

# Multiplication of two resolvents turns equality into a degree-two polynomial.
x, g, d, c, cp, alpha, alphap = sp.symbols("x g d c cp alpha alphap", real=True)
lhs = (c + alpha*x/(1-g*x) - cp - alphap*x/(1-d*x))*(1-g*x)*(1-d*x)
rhs = (c-cp) + (alpha-alphap-(c-cp)*(g+d))*x + ((c-cp)*g*d-alpha*d+alphap*g)*x**2
assert zero(lhs-rhs)
results["checks"]["krylov_collision_polynomial"] = {
    "exact": True,
    "coefficients": [str(sp.expand(rhs).coeff(x, j)) for j in range(3)],
}

# If P^2 B = s B + t P B, the cyclic span is invariant.  In its basis,
# P has the following exact 2 x 2 representation, for all such P and B.
s, t = sp.symbols("s t", real=True)
C = sp.Matrix([[0, s], [1, t]])
e1 = sp.Matrix([1, 0])
Dg = 1-t*g-s*g*g
Dd = 1-t*d-s*d*d
ng = C*(sp.eye(2)-g*C).inv()*e1
assert matrix_zero(ng-sp.Matrix([g*s, 1])/Dg)
alpha_new = alpha*Dd/Dg
c_new = c+alpha*s*(g-d)/Dg
mean_old = c*e1+alpha*ng
mean_new = c_new*e1+alpha_new*C*(sp.eye(2)-d*C).inv()*e1
assert matrix_zero(mean_old-mean_new)
results["checks"]["krylov_rank_two_compensation"] = {
    "exact": True,
    "resolvent_denominator": str(Dg),
    "alpha_at_new_discount": str(alpha_new),
    "c_at_new_discount": str(c_new),
}

# A rank-one channel map can nevertheless expose three propagation shapes.
J3 = sp.zeros(3)
J3[1, 0] = J3[2, 1] = 1
B3 = sp.Matrix([1, 0, 0])
Krylov3 = sp.Matrix.hstack(B3, J3*B3, J3**2*B3)
assert Krylov3.det() == 1
Y3 = c*B3+alpha*J3*(sp.eye(3)-g*J3).inv()*B3
assert matrix_zero(Y3-sp.Matrix([c, alpha, alpha*g]))
results["checks"]["krylov_rank_three_observation"] = {
    "exact": True,
    "neural_map_rank": int(B3.rank()),
    "propagation_shape_rank": int(Krylov3.rank()),
    "mean": [str(v) for v in Y3],
    "interpretation": "One channel suffices here for gamma after three distinct propagation shapes; this is not full-operator recovery.",
}

# The necessity construction for a genuine quarter-turn is exact.
J = sp.Matrix([[0, -1], [1, 0]])
B2 = sp.Matrix([2, 1])
g0, d0, alpha0, c0 = sp.Rational(1, 3), sp.Rational(4, 5), sp.Rational(5, 2), sp.Rational(7, 3)
alpha1 = alpha0*(1+d0*d0)/(1+g0*g0)
c1 = c0-alpha0*(g0-d0)/(1+g0*g0)
Y0 = c0*B2+alpha0*J*(sp.eye(2)-g0*J).inv()*B2
Y1 = c1*B2+alpha1*J*(sp.eye(2)-d0*J).inv()*B2
assert Y0 == Y1 and alpha1 > 0
results["checks"]["quarter_turn_discount_alias"] = {
    "exact": True,
    "old": {"gamma": str(g0), "alpha": str(alpha0), "c": str(c0)},
    "new": {"gamma": str(d0), "alpha": str(alpha1), "c": str(c1)},
    "identical_mean": [str(v) for v in Y0],
}

# Two different task frequencies: all coefficients must share c, alpha, gamma.
a1, a2 = sp.symbols("a1 a2", real=True)
D1, D2 = 1-2*a1*g+g*g, 1-2*a2*g+g*g
q = D2/D1
q_derivative = 2*(a1-a2)*(1-g*g)/D1**2
assert zero(sp.diff(q, g)-q_derivative)
t_recovered = (q-1)/(2*(q*a1-a2))
assert zero(t_recovered-g/(1+g*g))
results["checks"]["two_plane_inverse"] = {
    "exact_ratio_derivative": True,
    "exact_intermediate_inverse": True,
    "derivative": str(q_derivative),
}

# A realizable two-plane lack-of-fit example from the eight-state cycle.
# Using sqrt(2)/2, rather than an arbitrary unit-circle cosine, ensures that
# both blocks are Fourier modes of a finite stochastic task.
aa2 = sp.sqrt(2)/2
bb2 = sp.sqrt(2)/2
discount = sp.Rational(1,2)
dd1 = 1+discount**2
dd2 = 1-2*aa2*discount+discount**2
vv1 = 1/dd1
vv2 = bb2/dd2
cc1 = 1+discount/dd1
cc2 = 1-(aa2-discount)/dd2
assert sp.simplify(vv2-(8+10*sp.sqrt(2))/17) == 0
assert sp.simplify(cc2-(19-6*sp.sqrt(2))/17) == 0
assert sp.simplify(cc1-cc2-(24+30*sp.sqrt(2))/85) == 0
results["checks"]["eight_cycle_two_plane_lack_of_fit"] = {
    "exact": True,
    "gamma_forced_by_rotation_coefficients": "1/2",
    "alpha_forced_by_rotation_coefficients": 1,
    "u1": 1, "u2": 1,
    "v1": str(vv1), "v2": str(sp.radsimp(vv2)),
    "c1": str(cc1), "c2": str(sp.radsimp(cc2)),
    "realizability": "pi/2 and pi/4 Fourier blocks of an eight-state cycle",
}

# Exact full-pattern alias for free current + pure-future components.
Nplus = J*(sp.eye(2)-sp.Rational(1, 2)*J).inv()
Nminus = -J*(sp.eye(2)+sp.Rational(1, 2)*J).inv()
Gfp, Gfm = sp.eye(2)+5*Nplus, sp.eye(2)+5*Nminus
Gbp, Gbm = 3*sp.eye(2)+5*Nminus, 3*sp.eye(2)+5*Nplus
assert Gbp == -Gfp and Gbm == -Gfm
T_alias = Gfp*Gfm.inv()
assert T_alias == Gbp*Gbm.inv()
assert T_alias == (-15*sp.eye(2)-8*J)/17
results["checks"]["free_current_full_pattern_alias"] = {
    "exact": True,
    "discount": "1/2",
    "future_gain": 5,
    "forward_current": 1,
    "backward_current": 3,
    "reference_change": "B_prime = -B",
    "relative_operator": [[str(T_alias[i,j]) for j in range(2)] for i in range(2)],
}

# Task-compatible readout: verify intertwining for arbitrary 2 x 2 templates.
r11,r12,r21,r22,x11,x12,x21,x22 = sp.symbols("r11 r12 r21 r22 x11 x12 x21 x22", real=True)
R = sp.Matrix([[r11,r12],[r21,r22]])
X = sp.Matrix([[x11,x12],[x21,x22]])
inner = lambda A,B: sum(A[i,j]*B[i,j] for i in range(A.rows) for j in range(A.cols))
z = lambda A: sp.Matrix([inner(R,A), inner(J*R,A)])
assert matrix_zero(z(J*X)-J*z(X))
results["checks"]["compatible_readout_intertwining"] = {"exact": True, "normalization_omitted_from_both_sides": True}

# Partial identification: deterministic boundary table for manuscript figures.
rows = []
for phi_degrees in (0, 20, 40, 45, 50, 60, 70, 80):
    phi = math.radians(phi_degrees)
    lower = math.tan(max(0.0, phi-math.pi/4))
    upper = min(1.0, math.tan(phi)) if phi < math.pi/4 else 1.0
    if phi_degrees < 45:
        interval = "closed lower and upper"
    else:
        interval = "open lower and upper"
    rows.append({
        "phi_degrees": phi_degrees,
        "discount_infimum": lower,
        "discount_supremum": upper,
        "endpoint_status": interval,
        "unique_FF_label_under_open_unit_discount": phi_degrees >= 45,
        "positive_distance_from_switching_closure": phi_degrees > 45,
    })
with (OUT/"partial_discount_bounds.csv").open("w", encoding="utf-8", newline="") as f:
    writer=csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
results["checks"]["partial_discount_bound_table"] = {
    "source": "Analytic sharp identified interval, evaluated deterministically",
    "table": "partial_discount_bounds.csv",
    "rows": len(rows),
    "seventy_degree_lower_infimum": math.tan(math.radians(25)),
}

# The covariance-dimension distinction is arithmetic, not empirical calibration.
from scipy.stats import f as f_dist
n, error = 17, .05
K2=2*(n-1)/(n*(n-2))*f_dist.ppf(1-error,2,n-2)
K4=4*(n-1)/(n*(n-4))*f_dist.ppf(1-error,4,n-4)
assert K4 > K2
results["checks"]["hotelling_dimension_example"] = {
    "n": n, "alpha": error, "K2": float(K2), "K4": float(K4),
    "scope": "Threshold illustration only; n=17 does not certify actual repeat independence.",
}

results["all_checks_passed"] = True
(OUT/"math_additions.json").write_text(json.dumps(results, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
print(json.dumps(results, indent=2, ensure_ascii=False))
