"""Exact-algebra scope check: Krylov coordinates need not be SR lag moments."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import helmert

P=np.roll(np.eye(5),1,axis=1)
C=helmert(5,full=False)
A=C@P@C.T
B=np.eye(4)
g=.5
M=np.linalg.solve(np.eye(4)-g*A,B)
b=np.array([(g**k-g**4)/(1-g**5) for k in range(4)])
Y=sum(b[k]*np.linalg.matrix_power(A,k)@B for k in range(4))
D=np.column_stack([(np.linalg.matrix_power(A,k)@B).ravel() for k in range(4)])
positive=D[:,1:]-np.outer(D[:,0],D[:,0]@D[:,1:]/(D[:,0]@D[:,0]))
out={
  'model':'pure single-discount SR on a centered 5-cycle',
  'gamma':g,
  'C_P_minus_A_C_norm':float(np.linalg.norm(C@P-A@C)),
  'A4_plus_A3_plus_A2_plus_A_plus_I_norm':float(np.linalg.norm(sum(np.linalg.matrix_power(A,k) for k in range(5)))),
  'finite_basis_fit_residual':float(np.linalg.norm(M-Y)),
  'finite_basis_coefficients':b.tolist(),
  'positive_design_rank_after_current':int(np.linalg.matrix_rank(positive)),
  'fitted_coordinate_moment_determinant':float(b[1]*b[3]-b[2]**2),
  'analytic_fitted_determinant':float(-g**5*(1-g)**2/(1-g**5)**2),
  'true_infinite_lag_coefficients':[g,g*g,g**3],
  'true_lag_moment_determinant':float(g*g**3-(g*g)**2),
  'interpretation':'The finite regression coordinates are identifiable but are not the actual first three infinite-kernel lag weights. Applying their moment determinant would falsely reject a genuine single-discount SR.'}
assert out['positive_design_rank_after_current'] == 3
assert out['finite_basis_fit_residual'] < 1e-12
assert abs(out['fitted_coordinate_moment_determinant'] + 8/961) < 1e-14
assert out['true_lag_moment_determinant'] == 0
p=Path(__file__).resolve().parents[1]/'empirical/finite_kernel_moment_scope.json'
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
