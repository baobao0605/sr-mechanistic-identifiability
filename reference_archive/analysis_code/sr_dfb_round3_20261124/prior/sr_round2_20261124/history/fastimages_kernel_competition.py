"""Independently calibrated real sensory kernels, fair SR template competition.
All labels/times are actual stimuli. No simulated subjects and no mechanism truth inferred.
"""
from pathlib import Path
import json,hashlib,time
import numpy as np,pandas as pd
from scipy.optimize import nnls
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
threadpool_limits(2)
OUT=Path(__file__).resolve().parent
meta=json.loads((OUT/'fastimages_extraction_audit.json').read_text())
z=np.load(OUT/'fastimages_kernel_input.npz')
E=z['localizer_epochs'];lab=z['localizer_labels'];block=z['localizer_blocks'];sf=meta['sfreq']
trials=[]
for sid in range(len(z['trial_blocks'])):
 a,b=z['trial_indptr'][sid:sid+2]
 trials.append({'sequence_id':sid,'block':int(z['trial_blocks'][sid]),'onset_samples':z['trial_onset_samples'][sid], 'labels':z['trial_labels'][sid], 'sample_ids':z['trial_sample_ids_flat'][a:b], 'Y':z['trial_y_flat'][a:b], 'SOA_ms':float(z['trial_SOA_ms'][sid])})

def fit_nonnegative(G,b):
 v,U=eigh((G+G.T)/2,check_finite=False); use=v>max(v.max()*1e-12,1e-18)
 A=np.sqrt(v[use])[:,None]*U[:,use].T; yy=(U[:,use].T@b)/np.sqrt(v[use])
 coef,_=nnls(A,yy,maxiter=1000)
 return coef,float(coef@G@coef-2*coef@b)

def centered(s):
 G=s['G']-s['sx']@s['sx'].T/s['n']; b=s['b']-s['sx']@s['sy']/s['n']; return G,b

def aggregate(S,ii):
 return {k:sum(S[j][k] for j in ii) for k in ['G','b','yy','sx','sy','n']}

def score(s,c,off):
 return float(s['yy']-2*c@s['b']+c@s['G']@c-2*off@s['sy']+2*off@(c@s['sx'])+s['n']*(off@off))

def resid_distance(G,delta,nuis):
 # Schur complement gives squared distance after best unconstrained shared nuisance refit.
 A=G[np.ix_(nuis,nuis)]; v=G[nuis,:]@delta
 return max(0,float(delta@G@delta-v@np.linalg.pinv(A,rcond=1e-10)@v))

rows=[]; profiles=[]; geom=[]; kernel_rows=[]; kernels=[]; scales=[]; stats_records=[]; t0=time.time()
gammas=np.arange(1,10)/10
for fold in range(8):
 trainloc=block!=fold
 H=np.array([E[trainloc&(lab==j)].mean(0) for j in range(5)])
 residual=E[trainloc]-H[lab[trainloc]]
 sc=np.sqrt(np.mean(residual**2,axis=(0,1))); sc=np.maximum(sc,np.median(sc)*1e-6)
 H=H/sc[None,None,:]; generic=H.mean(0); specific=H-generic
 kernels.append(H);scales.append(sc)
 # Examine temporal-contrast rank and positive exponential separability later using saved H.
 S=[]
 for z in trials:
  Y=z['Y']/sc; ids=z['sample_ids']; on=z['onset_samples']; labels=z['labels']; N=len(ids)
  X=np.zeros((28,N,102)) # full generic5, full identity5, F1..F4/B1..B4, memoryless generic5/identity5
  for i in range(5):
   age=ids-on[i];use=(age>=0)&(age<=100); memuse=use.copy()
   if i<4:memuse &= ids<on[i+1]
   X[i,use]=generic[age[use]];X[5+i,use]=specific[labels[i],age[use]]
   X[18+i,memuse]=generic[age[memuse]];X[23+i,memuse]=specific[labels[i],age[memuse]]
   for k in range(1,5):
    if i+k<5:X[9+k,use]+=specific[labels[i+k],age[use]]
    if i-k>=0:X[13+k,use]+=specific[labels[i-k],age[use]]
  flat=X.reshape(28,-1);s={'G':flat@flat.T,'b':flat@Y.ravel(),'yy':float(np.sum(Y*Y)),'sx':X.sum(1),'sy':Y.sum(0),'n':N}
  S.append(s)
 S_tr=aggregate(S,[i for i,z in enumerate(trials) if z['block']!=fold]); S_te=aggregate(S,[i for i,z in enumerate(trials) if z['block']==fold]); G,b=centered(S_tr)
 coeffs={};chosen_g={};trainloss={}
 for name,inds in [('memoryless',np.arange(18,28)),('sensory_history',np.arange(10))]:
  c,l=fit_nonnegative(G[np.ix_(inds,inds)],b[inds]);cfull=np.zeros(28);cfull[inds]=c;coeffs[name]=cfull;chosen_g[name]=None;trainloss[name]=l
 for family,code in [('F',0),('B',1),('S',2)]:
  best=None
  for gamma in gammas:
   A=np.zeros((28,11));A[:10,:10]=np.eye(10); w=gamma**np.arange(1,5)
   if code in [0,2]:A[10:14,10]=w/(2 if code==2 else 1)
   if code in [1,2]:A[14:18,10]=w/(2 if code==2 else 1)
   c,l=fit_nonnegative(A.T@G@A,A.T@b); cf=A@c
   profiles.append({'fold':fold,'family':family,'gamma':gamma,'train_excess_quadratic_loss':l,'SR_gain':c[-1]})
   if best is None or l<best[0]:best=(l,cf,gamma,float(c[-1]))
  trainloss[family],coeffs[family],chosen_g[family],bg=best
 for name,c in coeffs.items():
  off=(S_tr['sy']-c@S_tr['sx'])/S_tr['n']; base=S_tr['sy']/S_tr['n']; overall=score(S_te,c,off); baseline=score(S_te,np.zeros(28),base)
  rows.append({'fold':fold,'SOA_ms':'all','model':name,'gamma':chosen_g[name],'SSE':overall,'baseline_SSE':baseline,'R2_vs_training_mean':1-overall/baseline,'train_excess_quadratic_loss':trainloss[name]})
  for nominal in [132,164,228,612]:
   ii=[i for i,z in enumerate(trials) if z['block']==fold and min([132,164,228,612],key=lambda a:abs(a-z['SOA_ms']))==nominal];st=aggregate(S,ii);ss=score(st,c,off); bb=score(st,np.zeros(28),base)
   rows.append({'fold':fold,'SOA_ms':nominal,'model':name,'gamma':chosen_g[name],'SSE':ss,'baseline_SSE':bb,'R2_vs_training_mean':1-ss/bb,'train_excess_quadratic_loss':trainloss[name]})
 for gamma in gammas:
  delta=np.zeros(28);delta[10:14]=gamma**np.arange(1,5);delta[14:18]=-gamma**np.arange(1,5)
  d0=resid_distance(G,delta,np.arange(18,28));d1=resid_distance(G,delta,np.r_[np.arange(10),np.arange(18,28)])
  geom.append({'fold':fold,'gamma':gamma,'distance_squared_after_memoryless_nuisance':d0,'distance_squared_after_physical_history_nuisance':d1,'remaining_distance_fraction':np.sqrt(d1/d0) if d0>0 else np.nan,'removed_information_fraction':1-d1/d0 if d0>0 else np.nan})
 stats_records.append({k:np.array([s[k] for s in S]) for k in S[0]})
 print('FAST fold',fold,'seconds',round(time.time()-t0,1),flush=True)
 pd.DataFrame(rows).to_csv(OUT/'fastimages_kernel_outer_folds.csv',index=False)
 pd.DataFrame(profiles).to_csv(OUT/'fastimages_kernel_gamma_profiles.csv',index=False)
 pd.DataFrame(geom).to_csv(OUT/'fastimages_kernel_discrimination_geometry.csv',index=False)
np.savez_compressed(OUT/'fastimages_kernel_sufficient_stats.npz',**{k:np.array([s[k] for s in stats_records]) for k in stats_records[0]},kernels=np.array(kernels),noise_scales=np.array(scales),trial_blocks=np.array([z['block'] for z in trials]),trial_SOA_ms=np.array([z['SOA_ms'] for z in trials]),trial_labels=np.array([z['labels'] for z in trials]))
r=pd.DataFrame(rows);s=r.groupby(['SOA_ms','model'],as_index=False)[['SSE','baseline_SSE']].sum();s['R2_vs_training_mean']=1-s.SSE/s.baseline_SSE
for soa in s.SOA_ms.unique():
 base=float(s[(s.SOA_ms==soa)&(s.model=='sensory_history')].SSE.iloc[0]);s.loc[s.SOA_ms==soa,'SSE_reduction_vs_sensory_history_pct']=100*(1-s.loc[s.SOA_ms==soa,'SSE']/base)
s.to_csv(OUT/'fastimages_kernel_summary.csv',index=False)
audit={'n_people':1,'n_blocks':8,'n_sequences':64,'n_ordinary_localizer':len(lab),'n_magnetometers':102,'sfreq':sf,'kernel_seconds':[0,1.0],'sequence_window':'first physical stimulus onset to final physical stimulus onset +1.0 s, own trial prestimulus baseline -0.2..0 s','source_first_time_seconds':meta['source_first_time_seconds'],'source_highpass_Hz':meta['source_highpass_Hz'],'source_lowpass_Hz':meta['source_lowpass_Hz'],'nuisance':'five position-specific generic evoked gains + five position-specific category contrast gains, constrained nonnegative; shared sensor offsets estimated only on training sequences. All candidate families use same nuisances and independent localizer kernels.','gamma_grid':gammas.tolist(),'family_scope':'Added activation of earlier/later physical sequence labels passed through same measured sensory kernels; this is operational model competition, not a claim of learned SR or definitive causal absence. Gamma trained only on seven sequence blocks; held block untouched.','geometry':'Unconstrained nuisance projection is conservative observation-space distance for fixed gamma and shared estimated kernel. Neural noise and kernel uncertainty are not certified by this number.','input_cache':'fastimages_kernel_input.npz contains real measured epochs in numeric arrays; no raw FIF or pickle needed to rerun.','source_sha256':meta['source_sha256']}
(OUT/'fastimages_kernel_audit.json').write_text(json.dumps(audit,indent=2))
print(s.to_string(index=False),flush=True)
