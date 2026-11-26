"""Targeted analytic regression for the repaired reference-ball solver.

This is numerical software verification, not simulated participants or an
additional neural experiment. Existing calibration CSVs are never changed.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import warnings

import numpy as np

from sr_geometry import reference_ball_residual

HERE = Path(__file__).resolve().parent
checks = []


def check(name, T, Bhat, Yhat, e0, expected):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = reference_ball_residual(T, Bhat, Yhat, e0)
    direct = float(np.linalg.norm(Yhat-T@(Bhat+result["delta"])))
    error = abs(result["residual"]-expected)
    # Relative to the expected residual, rather than to the large input
    # energy: the 1e-9 perpendicular component must not disappear.
    tolerance = max(2e-14*abs(expected), 5e-26)
    radius_tolerance = 4*np.finfo(float).eps*max(e0, np.finfo(float).tiny)
    record = dict(case=name, radius=e0, expected_residual=expected,
                  residual=result["residual"], absolute_error=error,
                  residual_tolerance=tolerance, direct_residual=direct,
                  radius_used=result["reference_error_used"],
                  warnings=[str(w.message) for w in caught])
    record["passed"] = bool(
        error <= tolerance and not caught
        and np.isfinite(result["residual"])
        and result["reference_error_used"] <= e0+radius_tolerance
        and abs(direct-result["residual"]) <= tolerance)
    checks.append(record)
    assert record["passed"], record


for shape in ("square", "rectangular"):
    T = np.diag([1., 0.]) if shape == "square" else np.array([[1.], [0.]])
    Bhat = np.zeros((T.shape[1], 1))
    Yhat = np.array([[1.], [1e-9]])
    for e0 in (0., .5, 1.-1e-12, 1., 2.):
        expected = math.hypot(max(1.-e0, 0.), 1e-9)
        check(f"{shape}_rank_deficient_radius_{e0:.13g}", T, Bhat, Yhat, e0, expected)

for e0 in (0., 2.):
    check(f"zero_operator_radius_{e0:g}", np.zeros((2, 3)), np.zeros((3, 1)),
          np.array([[3.], [4.]]), e0, 5.)

# Small but genuinely positive singular mode. Construct the radius from a
# known positive multiplier; the analytic optimal residual is 5e-11.
T = np.diag([1., 1e-10])
Yhat = np.array([[1.], [1e-10]])
lam = 1e-20
delta = np.array([[1./(1.+lam)], [.5]])
check("small_positive_singular_mode", T, np.zeros((2, 1)), Yhat,
      float(np.linalg.norm(delta)), float(np.linalg.norm(Yhat-T@delta)))

spec = importlib.util.spec_from_file_location(
    "reference_ball_v1_snapshot", HERE/"v1_before_reference_ball_fix/sr_geometry.py")
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
legacy = []
for name, T, e0 in (
    ("square_rank_deficient_exception", np.diag([1., 0.]), .5),
    ("rectangular_false_compatibility", np.array([[1.], [0.]]), 1.-1e-12),
):
    Bhat = np.zeros((T.shape[1], 1))
    Yhat = np.array([[1.], [1e-9]])
    record = dict(case=name)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            result = old.reference_ball_residual(T, Bhat, Yhat, e0)
            record.update(old_reported_residual=float(result["residual"]),
                          old_direct_residual=float(np.linalg.norm(Yhat-T@result["delta"])),
                          target_radius=1e-10,
                          old_reported_compatible=bool(result["residual"] <= 1e-10))
        except Exception as exc:
            record["old_exception"] = f"{type(exc).__name__}: {exc}"
        record["old_warnings"] = [str(w.message) for w in caught]
    fixed = reference_ball_residual(T, Bhat, Yhat, e0)
    record.update(repaired_residual=fixed["residual"],
                  repaired_compatible=bool(fixed["residual"] <= 1e-10))
    legacy.append(record)

output = dict(status="PASS", checks=checks, passed_checks=len(checks),
              legacy_counterexamples=legacy,
              solver_sha256=hashlib.sha256((HERE/"sr_geometry.py").read_bytes()).hexdigest(),
              scope="Targeted analytic numerical regression; no participant simulation; no statistical coverage claim.",
              independent_audit="../reviews/ekman_recheck_math_solver.json",
              real_result_comparison="../reviews/compatibility_after_solver.json")
(HERE/"reference_ball_v2_checks.json").write_text(json.dumps(output, indent=2)+"\n")
print(json.dumps({"status":output["status"], "passed_checks":len(checks),
                  "max_absolute_error":max(x["absolute_error"] for x in checks),
                  "solver_sha256":output["solver_sha256"]}, indent=2))
