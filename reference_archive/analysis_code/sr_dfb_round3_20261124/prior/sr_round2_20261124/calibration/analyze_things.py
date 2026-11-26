"""Frozen SR-discrimination calibration on real THINGS repeated responses.

All target transformations are imposed computations, not observed cognition.
The complete observed repeat is transformed (its noise is transformed too).
Four disjoint groups provide source-train,target-train,source-test,target-test.
Metric estimated from the source-train group's 3 repeats for every repeat budget.
"""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
from pathlib import Path
import json,hashlib,sys
import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf
from scipy.linalg import helmert
ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT.parent/"restored/complete"
sys.path.insert(0,str(ROOT/"math"))
from sr_geometry import residualize_current,ray_geometry,kernel_design
OUT=Path(__file__).resolve().parent
GAMMAS=np.round(np.arange(0,.951,.05),8)
KINDS=("forward","backward","symmetric_transition","balanced_fb")
STATES=(4,8,16,32,64,100)
CHANNELS=(6,12,24,48,72)
REPEATS=(1,2,3)
TRUTH_GAMMAS=(.25,.5,.8)
STRENGTHS=(.1,.3,1.)
z=np.load(PRIOR/"analysis/D05/analysis_input.npz")
Y=z["patterns"];roi=z["voxel_roi"]
assert Y.shape==(12,100,72)
groups=[np.arange(3*i,3*i+3) for i in range(4)]
folds=[(0,1,2,3),(1,2,3,0),(2,3,0,1),(3,0,1,2)]
def ops(n):
 P=np.roll(np.eye(n),-1,axis=1)
 I=np.eye(n);S=(P+P.T)/2
 ev,V=np.linalg.eigh(S)
 F=[];B=[];SY=[];BA=[]
 powers=[np.linalg.matrix_power(P,k+1) for k in range(n)]
 for g in GAMMAS:
  f=np.einsum("k,kij->ij",g**np.arange(n),powers)/(1-g**n)
  b=f.T
  F.append(f);B.append(b);SY.append((V*(ev/(1-g*ev)))@V.T)
  gb=g/(1+np.sqrt(1-g*g))
  fb=np.einsum("k,kij->ij",gb**np.arange(n),powers)/(1-gb**n)
  BA.append((fb+fb.T)/2)
 return P,np.array([F,B,SY,BA])
def fit_all(b,x,y):
 # b=(n*p,), x=(K,G,n*p), y=(Ntruth,n*p)
 xx=x.reshape(-1,len(b));bb=b@b;bx=xx@b/bb;by=y@b/bb
 r=xx-bx[:,None]*b;yr=y-by[:,None]*b;rr=(r*r).sum(1)
 alpha=np.maximum(0,np.divide(r@yr.T,rr[:,None],out=np.zeros((len(xx),len(y))),where=rr[:,None]>1e-20))
 c=by[None,:]-alpha*bx[:,None]
 # exact residual formula is stable enough, clamp roundoff
 sse=np.maximum(0,(yr*yr).sum(1)[None,:]-2*alpha*(r@yr.T)+alpha*alpha*rr[:,None])
 return sse.reshape(len(KINDS),len(GAMMAS),len(y)),c.reshape(len(KINDS),len(GAMMAS),len(y)),alpha.reshape(len(KINDS),len(GAMMAS),len(y))
rows=[];geo=[];rankrows=[];coeffrows=[];negative=[];foldmeta=[]
for n in STATES:
 P,O=ops(n)
 for nc in CHANNELS:
  ix=np.concatenate([np.where(roi==r)[0][:nc//3] for r in ("V1","lFFA","rFFA")])
  assert len(ix)==nc
  data=Y[:,:n,ix]
  data=data-data.mean(1,keepdims=True)
  H=helmert(n,full=False)
  for fold,(ga,gb,gc,gd) in enumerate(folds):
   resid=data[groups[ga]]-data[groups[ga]].mean(0)
   cov=LedoitWolf().fit(resid.reshape(-1,nc))
   ev,V=np.linalg.eigh(cov.covariance_)
   W=(V/np.sqrt(np.maximum(ev,1e-20)))@V.T
   d=data@W
   foldmeta.append(dict(n=n,channels=nc,fold=fold,metric_source_sessions=",".join(map(str,groups[ga]+1)),source_train_group=ga,target_train_group=gb,source_test_group=gc,target_test_group=gd,shrinkage=cov.shrinkage_))
   for nr in REPEATS:
    A,B,C,D=[d[groups[g][:nr]].mean(0) for g in (ga,gb,gc,gd)]
    XA=np.einsum("kgij,jp->kgip",O,A);XC=np.einsum("kgij,jp->kgip",O,C)
    RA=residualize_current(XA,A)
    obsrank=int(np.linalg.matrix_rank(H@A));sv=np.linalg.svd(H@A,compute_uv=False)
    kr=kernel_design(P,A,n_lags=3,drop_current=False)
    sk=np.linalg.svd(kr,compute_uv=False)
    rankrows.append(dict(states=n,channels=nc,repeats=nr,fold=fold,contrast_dimension=n-1,rank=obsrank,missing_contrasts=n-1-obsrank,smin_full=sv[-1] if obsrank==n-1 else 0,kernel4_rank=int(np.linalg.matrix_rank(kr)),kernel4_condition=sk[0]/sk[-1] if sk[-1]>sk[0]*1e-12 else np.inf,reference_split_frobenius=float(np.linalg.norm(A-C))))
    for aa,bb in [(0,1),(0,2),(0,3),(1,2),(2,3)]:
     r=ray_geometry(RA[aa],RA[bb])
     k=np.unravel_index(np.nanargmin(r["residual_per_unit"]),r["residual_per_unit"].shape)
     geo.append(dict(states=n,channels=nc,repeats=nr,fold=fold,pair=KINDS[aa]+"__"+KINDS[bb],grid_min_positive_ray_residual=float(r["residual_per_unit"][k]),gamma_a=GAMMAS[k[0]],gamma_b=GAMMAS[k[1]]/(1+np.sqrt(1-GAMMAS[k[1]]**2)) if bb==3 else GAMMAS[k[1]],description="grid minimum only"))
    truths=[("current",0.,0.,np.eye(n),np.nan)]
    for tk in (0,1,2):
     for tg in TRUTH_GAMMAS:
      gi=int(np.argmin(abs(GAMMAS-tg)));rr=np.linalg.norm(RA[tk,gi])
      for strength in STRENGTHS:
       a=strength*np.linalg.norm(A)/rr
       T=np.eye(n)+a*O[tk,gi]
       truths.append((KINDS[tk],tg,strength,T,a))
    T=np.stack([t[3] for t in truths])
    ytr=np.einsum("tij,jp->tip",T,B).reshape(len(T),-1)
    yte=np.einsum("tij,jp->tip",T,D).reshape(len(T),-1)
    sse,cs,alphas=fit_all(A.ravel(),XA,ytr)
    idx=sse.argmin(1)
    b=A.ravel();ct=(ytr@b)/(b@b)
    pred0=ct[:,None]*C.ravel()
    err0=((pred0-yte)**2).sum(1)
    scores=[err0];fitinfo=[]
    den=(yte*yte).sum(1)
    for ki,name in enumerate(KINDS):
     g=idx[ki];ar=np.arange(len(T));c=cs[ki,g,ar];a=alphas[ki,g,ar]
     x=XC[ki,g].reshape(len(T),-1)
     pred=c[:,None]*C.ravel()+a[:,None]*x
     err=((pred-yte)**2).sum(1);scores.append(err)
     fitinfo.append((g,c,a))
    scores=np.array(scores);winners=scores.argmin(0);labels=("current",)+KINDS
    for t,(tk,tg,strength,_,truea) in enumerate(truths):
     order=np.argsort(scores[:,t]);topgap=(scores[order[1],t]-scores[order[0],t])/den[t]
     tied=np.where(scores[:,t]<=scores[order[0],t]+1e-10*max(den[t],1e-30))[0]
     rawset=[labels[q] for q in tied]
     classes=[]
     for q in tied:
      if q==0:cls="current"
      else:
       gg,cv,av=fitinfo[q-1]
       activation=av[t]*np.linalg.norm(RA[q-1,gg[t]])/np.linalg.norm(A)
       cls="current" if activation<1e-10 else "symmetric" if q in (3,4) else labels[q]
      classes.append(cls)
     classes=sorted(set(classes));truth_class="symmetric" if tk=="symmetric_transition" else tk
     unique_correct=len(classes)==1 and classes[0]==truth_class
     rows.append(dict(states=n,channels=nc,repeats=nr,fold=fold,truth=tk,true_gamma=tg,predictive_strength_ratio=strength,true_predictive_gain=truea,winner=labels[winners[t]],correct_family=unique_correct,winner_set="|".join(rawset),winner_classes="|".join(classes),unique_class=len(classes)==1,truth_in_winner_classes=truth_class in classes,heldout_topgap=topgap,current_sse=err0[t]/den[t],full_row_rank=obsrank==n-1,empirical_reference_split_ratio=np.linalg.norm(A-C)/np.linalg.norm(A)))
     for ki,name in enumerate(KINDS):
      g,c,a=fitinfo[ki]
      gp=GAMMAS[g[t]];actual_gamma=gp/(1+np.sqrt(1-gp*gp)) if name=="balanced_fb" else gp
      active=a[t]*np.linalg.norm(RA[ki,g[t]])/np.linalg.norm(A)>=1e-10
      flat_profile=np.ptp(sse[ki,:,t])<=1e-10*max(np.sum(ytr[t]**2),1e-30)
      na_reason="no_visible_prediction" if not active else "flat_gamma_profile" if flat_profile else ""
      coeffrows.append(dict(states=n,channels=nc,repeats=nr,fold=fold,truth=tk,true_gamma=tg,predictive_strength_ratio=strength,candidate=name,estimated_gamma=actual_gamma if active and not flat_profile else np.nan,selected_grid_gamma=gp,active_predictive_fit=active,flat_gamma_profile_numerical=flat_profile,gamma_na_reason=na_reason,current_gain=c[t],predictive_gain=a[t],train_relative_sse=sse[ki,g[t],t]/np.sum(ytr[t]**2),heldout_relative_sse=scores[ki+1,t]/den[t],heldout_gain_over_current=(err0[t]-scores[ki+1,t])/den[t]))
    # Exact reversible-control geometry: P_rev=P_rev.T. F=B for any B.
    S=(P+P.T)/2
    NR=S@np.linalg.solve(np.eye(n)-.5*S,A)
    NB=S.T@np.linalg.solve(np.eye(n)-.5*S.T,A)
    negative.append(dict(states=n,channels=nc,repeats=nr,fold=fold,reversible_F_minus_B=float(np.linalg.norm(NR-NB))))
  print("finished",n,nc,flush=True)
pd.DataFrame(rows).to_csv(OUT/"recovery_trials.csv",index=False)
pd.DataFrame(coeffrows).to_csv(OUT/"model_scores.csv",index=False)
pd.DataFrame(geo).to_csv(OUT/"model_geometry.csv",index=False)
pd.DataFrame(rankrows).to_csv(OUT/"rank_and_kernel.csv",index=False)
pd.DataFrame(negative).to_csv(OUT/"reversible_control.csv",index=False)
pd.DataFrame(foldmeta).to_csv(OUT/"split_and_metric.csv",index=False)
df=pd.DataFrame(rows);cc=pd.DataFrame(coeffrows)
summary={"data":{"subjects":1,"sessions":12,"images":100,"voxels":72},"n_calibration_cases":len(df),"n_model_scores":len(cc),"scope":"Real-repeat imposed transforms; no inferred neural SR; correlated cases do not count as independent people. Tie-aware unique class recovery; symmetric_transition and balanced_fb share exact reparameterized family; fitted alpha=0 is current-only.","source_sha256":hashlib.sha256((PRIOR/"analysis/D05/analysis_input.npz").read_bytes()).hexdigest(),"primary100_72_rep3":df.query("states==100 and channels==72 and repeats==3").groupby(["truth","predictive_strength_ratio"]).agg(correct=("correct_family","mean"),n=("correct_family","size")).reset_index().to_dict("records"),"reversible_exact_max":max(r["reversible_F_minus_B"] for r in negative)}
(OUT/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2))
print(json.dumps(summary,ensure_ascii=False,indent=2))

