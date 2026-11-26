"""Separate numerical formulations for the new Ekman computations."""
from pathlib import Path
import ast,json,hashlib,io
import numpy as np,pandas as pd
from scipy.linalg import null_space
from scipy.optimize import lsq_linear
from scipy.stats import gamma
from scipy.signal import fftconvolve
from threadpoolctl import threadpool_limits
threadpool_limits(2)
O=Path(__file__).resolve().parent;A=O/'input/analysis_scripts';F=A/'02_Analyses_and_article_Figures'
# Execute only the isolated, inspected author's pure mapping function.
src=(A/'fx__PONS_helper_functions.py').read_text();start=src.index('def return_sequence_idxp_idxc');end=src.find('\ndef ',start+5)
ns={};exec(compile(src[start:end],'<author_sequence_mapping>','exec'),ns)
mapping=ns['return_sequence_idxp_idxc']
subject=np.loadtxt(A/'subjects_final.txt',dtype=str);seq=dict(np.loadtxt(A/'sequence_ids.txt',dtype=str))
files={'HP_Figure6':F/'Figure_06/data/N35_HPa_localizer_LR_loo_probs.npy','V1_Figure6':F/'Figure_06/data/N35_V1_localizer_LR_loo_probs.npy','HP_Figure5_sensitivity':F/'Figure_05/data/N35_HP_localizer_LR_loo_probs.npy'}
co=pd.read_csv(O/'localizer_unique_subject.csv');er=[]
for name,path in files.items():
    arr=np.load(path)
    for s,sub in enumerate(subject):
        idx,_=mapping(int(seq[sub]));i,j=np.where(~np.eye(8,dtype=bool));pos=np.full(8,-1);pos[idx]=range(4)
        isseq=(pos[i]>=0)&(pos[j]>=0);dist=np.minimum(abs(i-j),8-abs(i-j));lag=pos[j]-pos[i]
        Z=np.column_stack([np.ones(56)]+[(i==k)*1. for k in range(1,8)]+[(j==k)*1. for k in range(1,8)]+[(dist==k)*1. for k in range(1,4)]+[isseq*1.])
        T=-abs(lag)*isseq;U=np.sign(lag)*isseq;y=arr[s,i,j]
        for metric,x,other in [('temporal_proximity_unique',T,U),('direction_unique',U,T)]:
            N=np.column_stack([Z,other]);Q=null_space(N.T)
            xx=Q.T@x;yy=Q.T@y;beta=xx@yy/(xx@xx)
            old=co.query('source==@name and subject==@sub and metric==@metric').value.iloc[0]
            er.append(dict(source=name,subject=sub,metric=metric,beta_nullspace=beta,beta_OLS=old,absolute_difference=abs(beta-old)))
pd.DataFrame(er).to_csv(O/'independent_coefficient_checks.csv',index=False)

# Bounded-variable least squares at the selected grid value, contrasting the
# primary residualization implementation with a joint optimization formulation.
pr=pd.read_csv(O/'localizer_heldcue_predictions.csv');bv=[]
for name,path in files.items():
    arr=np.load(path)
    for si in [0,34]:
        sub=subject[si];idx,_=mapping(int(seq[sub]));i,j=np.where(~np.eye(8,dtype=bool));pos=np.full(8,-1);pos[idx]=range(4)
        isseq=(pos[i]>=0)&(pos[j]>=0);dist=np.minimum(abs(i-j),8-abs(i-j));lag=pos[j]-pos[i]
        Z=np.column_stack([np.ones(56)]+[(j==k)*1. for k in range(1,8)]+[(dist==k)*1. for k in range(1,4)]+[isseq*1.]);y=arr[si,i,j]
        for held in [0,7]:
            tr=i!=held;te=~tr
            for model in ['forward','backward','symmetric']:
                rr=pr.query('source==@name and subject==@sub and held_cue==@held and model==@model').iloc[0];g=rr.gamma
                active=isseq & ((lag>0) if model=='forward' else (lag<0) if model=='backward' else np.ones(56,bool))
                x=np.where(active,g**np.maximum(abs(lag)-1,0),0);X=np.column_stack([Z,x]);p=X.shape[1]
                res=lsq_linear(X[tr],y[tr],bounds=(np.r_[np.full(p-1,-np.inf),0],np.full(p,np.inf)),method='bvls',tol=1e-12)
                sse=np.sum((y[te]-X[te]@res.x)**2)
                bv.append(dict(source=name,subject=sub,held=held,model=model,SSE_BVLS=sse,SSE_primary=rr.test_SSE,absolute_difference=abs(sse-rr.test_SSE)))
pd.DataFrame(bv).to_csv(O/'independent_bvls_checks.csv',index=False)

# Independently approximate convolution by a 10ms sampled impulse response,
# comparing all measured event onsets with the analytic CDF implementation.
saved=np.load(O/'actual_history_mixing_matrices.npz');conv=[]
for sub in ['S01','S02']:
    ev=O/f'input/analysis_scripts/subject_data/{sub}/evs';dt=.01;tt=np.arange(0,1020,dt);ht=np.arange(0,90,dt);h=gamma.pdf(ht,a=4,scale=1.5)
    G=[];times=np.arange(666)*1.5+.75
    # Matrix row length is taken from actual design (not assumed same across people).
    dmat=np.loadtxt(io.StringIO((O/f'input/data/{sub}/func/localizer_hp100_fwhm5_v3.feat/design.mat').read_text().split('/Matrix')[1]));times=(np.arange(len(dmat))+.5)*1.5
    all_on=[]
    for k in range(8):
        es=np.loadtxt(ev/f'localizer_dot{k+1}.txt');sig=np.zeros(len(tt));all_on.append(es[:,0])
        for on,dur,amp in es:sig[(tt>=on)&(tt<on+dur)]+=amp
        v=fftconvolve(sig,h)[:len(tt)]*dt;G.append(np.interp(times,tt,v))
    G=np.column_stack(G);C=np.zeros((8,8))
    for k in range(8):C[k]=np.mean([G[int(np.floor(t/1.5))+2:int(np.floor(t/1.5))+9].mean(0) for t in all_on[k]],axis=0)
    target=saved[f'{sub}_6.0_unfiltered_LTI_3to13p5s'];err=np.max(abs(C-target))
    conv.append(dict(subject=sub,dt_s=dt,max_absolute_C_error=err))
pd.DataFrame(conv).to_csv(O/'independent_HRF_convolution_checks.csv',index=False)

manifest=O/'repository_MANIFEST.txt';manifest_checks=[]
if manifest.exists():
    lookup={line.split(' ',1)[1]:line.split(' ',1)[0] for line in manifest.read_text().splitlines() if ' ' in line}
    for path in [*files.values(),F/'Figure_06/Figure_06_tuning.ipynb',A/'subjects_final.txt',A/'sequence_ids.txt']:
        rel=str(path.relative_to(O/'input'));actual=hashlib.sha256(path.read_bytes()).hexdigest()
        manifest_checks.append(dict(path=rel,sha256=actual,official_sha256=lookup.get(rel),matches_official=actual==lookup.get(rel)))
out={'coefficient_nullspace_max_error':max(x['absolute_difference'] for x in er),'n_coefficient_checks':len(er),
     'joint_BVLS_max_SSE_error':max(x['absolute_difference'] for x in bv),'n_BVLS_checks':len(bv),
     'convolution_10ms_max_error':max(x['max_absolute_C_error'] for x in conv),'official_manifest_checks':manifest_checks}
assert out['coefficient_nullspace_max_error']<1e-12
assert out['joint_BVLS_max_SSE_error']<1e-10
assert out['convolution_10ms_max_error']<.002
assert all(x['matches_official'] for x in manifest_checks)
(O/'independent_checks.json').write_text(json.dumps(out,ensure_ascii=False,indent=2));print(json.dumps(out,ensure_ascii=False,indent=2))
