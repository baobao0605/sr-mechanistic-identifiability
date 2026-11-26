"""Retrospective, held-out event-lag encoding on real Allen natural scene counts.
No synthetic neural observations or assumed exponential response kernels.
"""
from pathlib import Path
import numpy as np, pandas as pd, json, time, hashlib
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
threadpool_limits(2)
BASE=Path(__file__).resolve().parent/'inputs'
OUT=Path(__file__).resolve().parent
z=np.load(BASE/'natural_scene_counts.npz')
y=z['state']; t=z['start']; n=len(y); groups=np.minimum(np.arange(n)*4//n,3)
cont=np.r_[False,(np.diff(t)<.3)&(np.diff(z['block'])==0)]
# Uniform three-event guard covers every compared predictor and group boundary.
elig=z['valid'].copy()&(y>=0)
for k in range(1,4):
 elig &= np.r_[np.zeros(k,bool), np.logical_and.reduce([cont[j:n-k+j] for j in range(1,k+1)])]
 elig &= np.r_[groups[:k]!=groups[:k],groups[k:]==groups[:-k]]
 elig &= np.r_[np.logical_and.reduce([cont[j:n-k+j] for j in range(1,k+1)]),np.zeros(k,bool)]
 elig &= np.r_[groups[:-k]==groups[k:],groups[-k:]!=groups[-k:]]
idx=np.where(elig)[0]; labels=y+1
YY=z['response'].astype(float); rois={'VISp':z['areas']=='VISp','all_visual':np.ones(len(z['areas']),bool)}
alpha_grid=np.array([1.,10.,100.,1000.,10000.,100000.])
models={'C':[0],'H1':[0,-1],'H3':[0,-1,-2,-3],'F1':[0,1]}

def design(lags):
 return np.concatenate([np.eye(119)[labels[idx+k]] for k in lags],axis=1)

def fit_grid(Xtr,Ytr,Xva,Yva,history_scale=1.):
 xm=Xtr.mean(0); xs=Xtr.std(0); xs[xs<1e-12]=1
 ym=Ytr.mean(0); ys=np.sqrt(Ytr.var(0)+1)
 A=(Xtr-xm)/xs; V=(Xva-xm)/xs; A[:,119:]*=history_scale; V[:,119:]*=history_scale; B=(Ytr-ym)/ys; Z=(Yva-ym)/ys
 vals,U=eigh(A.T@A,check_finite=False); vals=np.maximum(vals,0)
 UTXY=U.T@(A.T@B); VU=V@U
 sse=[]
 for a in alpha_grid:
  pred=VU@(UTXY/(vals[:,None]+a)); sse.append(np.sum((Z-pred)**2,axis=0))
 return np.array(sse),np.sum(Z*Z,axis=0),np.sqrt(vals[::-1]),(xm,xs,ym,ys,U,vals,UTXY)

rows=[]; inner=[]; spectra=[]; coef_records={}; st=time.time()
for name,lags in models.items():
 X=design(lags); Y=YY[idx]; g=groups[idx]
 for test in range(4):
  tr=g!=test; te=g==test; training_groups=np.unique(g[tr]); inners=[]
  scales=[0.] if name=='C' else [0.,.01,.1,1.]
  inner_by_scale=[]
  for hs in scales:
   inners=[]
   for val in training_groups:
    aa=tr&(g!=val); bb=g==val
    ss,base,sv,_=fit_grid(X[aa],Y[aa],X[bb],Y[bb],hs); inners.append(ss)
   inner_by_scale.append(np.sum(inners,axis=0))
  total=np.concatenate(inner_by_scale,axis=0)
  chosen_by_roi={roi:int(np.argmin(total[:,sel].sum(1))) for roi,sel in rois.items()}
  outer_by_scale={}
  for roi,sel in rois.items():
   chosen=chosen_by_roi[roi]; hs=scales[chosen//len(alpha_grid)]; ai=chosen%len(alpha_grid)
   if hs not in outer_by_scale: outer_by_scale[hs]=fit_grid(X[tr],Y[tr],X[te],Y[te],hs)
   ss,base,sv,fit=outer_by_scale[hs]
   for ii in range(len(total)):
    inner.append({'model':name,'test_group':test,'roi':roi,'history_scale':scales[ii//len(alpha_grid)],'alpha':alpha_grid[ii%len(alpha_grid)],'inner_SSE':total[ii,sel].sum(),'selected':ii==chosen})
   rows.append({'model':name,'test_group':test,'roi':roi,'n_test':int(te.sum()),'n_train':int(tr.sum()),'n_units':int(sel.sum()),'history_scale':hs,'alpha':alpha_grid[ai],'SSE_standardized':float(ss[ai,sel].sum()),'baseline_SSE_standardized':float(base[sel].sum()),'SSE_over_baseline':float(ss[ai,sel].sum()/base[sel].sum()),'R2_to_training_mean':float(1-ss[ai,sel].sum()/base[sel].sum())})
  # Unpenalized predictor geometry; centered categorical redundancies are expected.
  xs=X[tr].std(0); xs[xs<1e-12]=1
  sv=np.linalg.svd((X[tr]-X[tr].mean(0))/xs,compute_uv=False)
  rank=int((sv>sv[0]*1e-10).sum())
  for j,v in enumerate(sv): spectra.append({'model':name,'test_group':test,'index':j,'singular_value':v,'rank':rank})
  print(name,test,round(time.time()-st,1),flush=True)
 pd.DataFrame(rows).to_csv(OUT/'allen_history_outer_folds.csv',index=False)
 pd.DataFrame(inner).to_csv(OUT/'allen_history_inner_selection.csv',index=False)
 pd.DataFrame(spectra).to_csv(OUT/'allen_history_design_spectra.csv',index=False)
res=pd.DataFrame(rows)
summ=res.groupby(['roi','model'],as_index=False)[['SSE_standardized','baseline_SSE_standardized','n_test']].sum()
summ['R2_to_training_mean']=1-summ.SSE_standardized/summ.baseline_SSE_standardized
for roi in rois:
 c=float(summ[(summ.roi==roi)&(summ.model=='C')].SSE_standardized.iloc[0]);summ.loc[summ.roi==roi,'SSE_reduction_vs_current_pct']=100*(1-summ.loc[summ.roi==roi,'SSE_standardized']/c)
summ.to_csv(OUT/'allen_history_summary.csv',index=False)
audit={'source_sha256':hashlib.sha256((BASE/'natural_scene_counts.npz').read_bytes()).hexdigest(),'n_events_all':n,'n_eligible':len(idx),'group_counts':dict(zip(*[a.tolist() for a in np.unique(groups[idx],return_counts=True)])),'excluded_by_uniform_guard_from_previous_nonblank':int((z['valid']&(y>=0)).sum()-len(idx)),'lag_units':'one actual event, median spacing seconds '+str(np.median(np.diff(t))), 'alpha_grid':alpha_grid.tolist(),'history_scales':[0.,.01,.1,1.],'amendment':'regularization_amendment.json','methods':'All models use exact same outer observations and three-event guard. Each inner fit standardizes X and Y from that inner training split only. No p-values: one mouse; four folds are not four independent animals. Regularized categorical encoding in real post-stimulus counts, not an SR fitted brain mechanism.'}
(OUT/'allen_history_audit.json').write_text(json.dumps(audit,indent=2))
print(summ.to_string(index=False),flush=True)
