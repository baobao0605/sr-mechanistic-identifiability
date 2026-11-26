"""Cross-repeat checks for noise-inflated SR template geometry.
Independent held-out reference groups C,D; feature metric trained only on A.
Cross energies are unbiased under independent, zero-mean repeat noise and a
stable latent reference. They are not population inference from one person.
"""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.covariance import LedoitWolf
R=Path(__file__).resolve().parents[1]
z=np.load(R.parent/"restored/complete/analysis/D05/analysis_input.npz")
groups=[np.arange(3*i,3*i+3) for i in range(4)]
rows=[]
for n in [4,8,16,32,64,100]:
 I=np.eye(n);P=np.roll(I,-1,axis=1)
 F=P@np.linalg.inv(I-.5*P);B=P.T@np.linalg.inv(I-.5*P.T);G=P@np.linalg.inv(I-.7*P)
 for nc in [6,12,24,48,72]:
  ix=np.r_[tuple(np.where(z["voxel_roi"]==roi)[0][:nc//3] for roi in ["V1","lFFA","rFFA"])]
  x=z["patterns"][:,:n,ix];x=x-x.mean(1,keepdims=True)
  for fold in range(4):
   ids=[(fold+k)%4 for k in range(4)];tr=x[groups[ids[0]]]
   cov=LedoitWolf().fit((tr-tr.mean(0)).reshape(-1,nc)).covariance_
   ev,V=np.linalg.eigh(cov);W=(V/np.sqrt(ev))@V.T
   A,C,D=[x[groups[ids[k]]].mean(0)@W for k in [0,2,3]]
   aa=A.ravel();gg=(G@A).ravel();ff=(F@A).ravel()
   rg=gg-aa*(aa@gg)/(aa@aa);rf=ff-aa*(aa@ff)/(aa@aa)
   a=max(0,(rg@rf)/(rg@rg));c=(aa@(ff-a*gg))/(aa@aa)
   for name,T in [("current_reference",I),("fixed_forward_minus_backward",F-B),("gamma_point5_vs_point7_after_train_gain",F-c*I-a*G)]:
    U=T@C;Vv=T@D
    cross=float(np.sum(U*Vv));naive=float((np.sum(U*U)+np.sum(Vv*Vv))/2)
    rows.append(dict(states=n,channels=nc,fold=fold,contrast=name,cross_repeat_energy=cross,within_repeat_energy=naive,cross_to_within_fraction=cross/naive,source_train_sessions=",".join(map(str,groups[ids[0]]+1)),heldout_C_sessions=",".join(map(str,groups[ids[2]]+1)),heldout_D_sessions=",".join(map(str,groups[ids[3]]+1))))
df=pd.DataFrame(rows);df.to_csv(R/"calibration/cross_repeat_contrast_energy.csv",index=False)
print(df.query("states==100 and channels==72").groupby("contrast").agg(cross_min=("cross_repeat_energy","min"),cross_max=("cross_repeat_energy","max"),fraction_min=("cross_to_within_fraction","min"),fraction_max=("cross_to_within_fraction","max")).to_string())

