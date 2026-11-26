"""Reference-error compatibility on fixed candidate transforms, real repeats.
Error radii are sensitivity multipliers of split-repeat differences, NOT CIs.
"""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
from scipy.optimize import brentq,minimize
from sklearn.covariance import LedoitWolf
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"math"))
from sr_geometry import reference_ball_residual,residualize_current
z=np.load(ROOT.parent/"restored/complete/analysis/D05/analysis_input.npz")
x=z["patterns"]; x=x-x.mean(1,keepdims=True); n=100
P=np.roll(np.eye(n),-1,axis=1);I=np.eye(n);S=(P+P.T)/2;g=.5
ops=[P@np.linalg.inv(I-g*P),P.T@np.linalg.inv(I-g*P.T),S@np.linalg.inv(I-g*S)]
groups=[np.arange(3*i,3*i+3) for i in range(4)]
out=[];summary=[];checks=[]
for fold in range(4):
 ids=[(fold+k)%4 for k in range(4)]
 src=x[groups[ids[0]]]
 cov=LedoitWolf().fit((src-src.mean(0)).reshape(-1,72)).covariance_
 ev,V=np.linalg.eigh(cov);W=(V/np.sqrt(ev))@V.T
 A,B,C,D=[x[groups[i]].mean(0)@W for i in ids]
 a=np.linalg.norm(A)/np.linalg.norm(residualize_current(ops[0]@A,A))
 candidates=[I,I+a*ops[0],I+a*ops[1],I+a*ops[2]]
 for truth_id in [0,1]:
  Y=candidates[truth_id]@B
  e0=np.linalg.norm(A-C)/np.sqrt(2);e1=np.linalg.norm(candidates[truth_id]@(B-D))/np.sqrt(2)
  for name,T in zip(["current","forward","backward","symmetric"],candidates):
   def fun(t):return reference_ball_residual(T,A,Y,t*e0)["residual"]-t*e1
   high=1.
   while fun(high)>0:high*=2
   threshold=brentq(fun,0,high,xtol=1e-10)
   summary.append(dict(fold=fold,truth=["current","forward"][truth_id],candidate=name,minimum_error_multiplier=threshold,empirical_e0=e0,empirical_e1=e1,gamma=g,predictive_gain=a,scope="fixed candidates only, empirical noise scales not confidence bounds"))
   for t in [0,.25,.5,.75,1,1.5,2]:
    r=reference_ball_residual(T,A,Y,t*e0)
    necessary=max(0,np.linalg.norm(Y-T@A)-np.linalg.norm(T,2)*t*e0)
    out.append(dict(fold=fold,truth=["current","forward"][truth_id],candidate=name,multiplier=t,exact_min_residual=r["residual"],target_error_radius=t*e1,compatible=r["residual"]<=t*e1,triangle_lower_bound=necessary,primal_check=np.linalg.norm(Y-T@(A+r["delta"])),reference_error_used=r["reference_error_used"],reference_error_radius=t*e0))
# independent low-dimensional constrained optimizer, actual measured entries
T=np.array([[1.,.4,0],[0,1.,.4],[0,0,1.]])
A=x[:3,:3,:2].mean(0);Y=x[3:6,:3,:2].mean(0)
for radius in [.01,.1,1.]:
 ans=reference_ball_residual(T,A,Y,radius)
 f=lambda z:np.sum((Y-T@(A+z.reshape(A.shape)))**2)
 opt=minimize(f,np.zeros(A.size),method="SLSQP",constraints=[{"type":"ineq","fun":lambda z:radius**2-z@z}],options={"maxiter":3000,"ftol":1e-13})
 checks.append(dict(radius=radius,optimizer_success=bool(opt.success),reference_ball_sse=ans["residual"]**2,independent_slsqp_sse=float(opt.fun),sse_difference=abs(ans["residual"]**2-opt.fun),constraint_violation=max(0,float(np.linalg.norm(opt.x)-radius))))
pd.DataFrame(out).to_csv(ROOT/"calibration/reference_compatibility_curves.csv",index=False)
pd.DataFrame(summary).to_csv(ROOT/"calibration/reference_compatibility_thresholds.csv",index=False)
(ROOT/"reviews/reference_ball_independent_checks.json").write_text(json.dumps(checks,indent=2))
print(pd.DataFrame(summary).groupby(["truth","candidate"]).minimum_error_multiplier.agg(["min","max"]).to_string())
print(json.dumps(checks,indent=2))

