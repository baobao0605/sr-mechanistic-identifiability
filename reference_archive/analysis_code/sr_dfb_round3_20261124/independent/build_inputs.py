"""Reconstruct immutable Kern/Turner folds without executing the old analysis."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
import pandas as pd
from scipy.linalg import helmert
from sklearn.covariance import LedoitWolf

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
PRIOR=ROOT/'prior/sr_round2_20261124'
OLD=ROOT/'prior/restored/complete'
sys.path.insert(0,str(PRIOR/'math'))
from sr_geometry import templates

POST=['P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4','PO8','O1','Oz','O2']
KEYS=['dataset','window','fold','space']

def means(x,y,n):
    return np.array([x[y==i].mean(0) for i in range(n)])

def halves(x,y,n):
    ix=[[],[]]
    for j in range(n):
        ind=np.flatnonzero(y==j);mid=len(ind)//2
        ix[0].extend(ind[:mid]);ix[1].extend(ind[mid:])
    return [(x[np.sort(k)],y[np.sort(k)]) for k in ix]

def fourier(n,ks):
    t=2*np.pi*np.arange(n)/n
    return np.array([f(k*t) for k in ks for f in (np.cos,np.sin)])*np.sqrt(2/n)

def whiten_fit(x,y,n):
    mu=means(x,y,n);res=x-mu[y]
    cov=LedoitWolf().fit(res).covariance_
    vals,U=np.linalg.eigh(cov)
    return (U/np.sqrt(np.maximum(vals,max(float(vals.max()),1e-300)*1e-12)))@U.T

def get_data():
    x={};seq=np.array([ord(s)-65 for s in 'ABCDEFGEHIBJ'])
    kp=OLD/'prior/kern/derived/shared_ssp'
    for sess in (1,2):
        z=np.load(kp/f'localizer{sess}_filtered_features.npz')
        d=pd.read_csv(kp/f'localizer{sess}_trials.csv')
        keep=d.standard.to_numpy(bool)&~d.ssp_artifact_flag.to_numpy(bool)
        for win in ('visual_early','visual_late'):
            x[sess,win]=(z[win][keep].astype(float),d.label.to_numpy()[keep])
    yield 'Kern',x,10,np.eye(10)[seq],{'full11':helmert(12,full=False),'design_k4_2D':fourier(12,[4]),'design_k4_k1_4D':fourier(12,[4,1])}
    x={}
    for sess in (1,2):
        z=np.load(OLD/f'analysis/D08/session{sess}_localizer.npz')
        ix=[list(z['channels']).index(c) for c in POST]
        for wi,win in ((1,'early'),(2,'late')):
            x[sess,win]=((z['features'][:,wi]-z['baseline'])[:,ix],z['labels'])
    yield 'Turner',x,40,np.eye(40),{'full39':helmert(40,full=False),'Fourier1_2D':fourier(40,[1]),'Fourier12_4D':fourier(40,[1,2])}

def kernel(A,g,kind):
    if kind=='current':return np.zeros_like(A)
    if kind=='forward':T=A
    elif kind=='backward':T=A.T
    elif kind=='symmetric_transition':T=(A+A.T)/2
    elif kind=='balanced_fb':
        return .5*(kernel(A,g,'forward')+kernel(A,g,'backward'))
    else:raise KeyError(kind)
    return np.linalg.solve(np.eye(len(A))-g*T,T)

def main():
    HERE.mkdir(exist_ok=True)
    old=pd.read_csv(PRIOR/'math/calibration_folds.csv')
    arrays={};meta=[];audit=[]
    for dataset,data,n,L,spaces in get_data():
        P=np.roll(np.eye(len(L)),1,axis=1)
        for win in sorted(set(w for s,w in data)):
            split={s:halves(*data[s,win],n) for s in (1,2)}
            for source,target in ((1,2),(2,1)):
                for trainhalf in (0,1):
                    fold=f'{source}to{target}_train{trainhalf}'
                    xa,ya=split[source][trainhalf];xc,yc=split[source][1-trainhalf]
                    xb,yb=split[target][trainhalf];xd,yd=split[target][1-trainhalf]
                    W=whiten_fit(xa,ya,n)
                    RR=L@means(xa@W,ya,n);SS=L@means(xc@W,yc,n)
                    VV=L@means(xb@W,yb,n);ZZ=L@means(xd@W,yd,n)
                    for space,C in spaces.items():
                        key={"dataset":dataset,"window":win,"fold":fold,"space":space}
                        fid=f'fold{len(meta):02d}';A=C@P@C.T
                        R,V,S,Z=[C@v for v in (RR,VV,SS,ZZ)]
                        for name,v in zip(('A','R','V','S','Z'),(A,R,V,S,Z)):arrays[f'{fid}_{name}']=v
                        key.update(fold_id=fid,n_contrasts=len(A),n_channels=R.shape[1],reduction_error=float(np.linalg.norm(C@P-A@C)))
                        meta.append(key)
                        mask=np.ones(len(old),dtype=bool)
                        for k in KEYS:mask &=old[k].eq(key[k]).to_numpy()
                        part=old[mask]
                        assert len(part)==95
                        for idx,row in part.iterrows():
                            T0=np.eye(len(A))+float(row.true_alpha)*kernel(A,float(row.true_gamma),row.true_kind)
                            Ytr=T0@V;Yte=T0@Z
                            if row.candidate=='current':Tfit=float(row.current_gain)*np.eye(len(A))
                            else:Tfit=float(row.current_gain)*np.eye(len(A))+float(row.predictive_gain)*kernel(A,float(row.gamma_grid_argmin_audit_only),row.candidate)
                            sse=float(np.sum((Yte-Tfit@S)**2));trsse=float(np.sum((Ytr-Tfit@R)**2))
                            audit.append({'old_row':int(idx),'fold_id':fid,'candidate':row.candidate,'test_sse_old':row.test_sse,'test_sse_rebuilt':sse,'train_sse_old':row.train_sse,'train_sse_rebuilt':trsse})
        print('Built',dataset,flush=True)
    np.savez_compressed(HERE/'repeated_inputs.npz',**arrays)
    pd.DataFrame(meta).to_csv(HERE/'input_folds.csv',index=False)
    d=pd.DataFrame(audit);d.to_csv(HERE/'input_reconstruction_audit.csv',index=False)
    tr=float(np.max(np.abs(d.train_sse_old-d.train_sse_rebuilt)));te=float(np.max(np.abs(d.test_sse_old-d.test_sse_rebuilt)))
    assert np.allclose(d.test_sse_old,d.test_sse_rebuilt,atol=1e-9,rtol=1e-10)
    assert np.allclose(d.train_sse_old,d.train_sse_rebuilt,atol=1e-9,rtol=1e-10)
    summary={'folds':len(meta),'model_scores_rebuilt':len(audit),'test_sse_max_abs_error':te,'train_sse_max_abs_error':tr,'status':'PASS','old_scores_sha256':hashlib.sha256((PRIOR/'math/calibration_folds.csv').read_bytes()).hexdigest(),'cache_sha256':hashlib.sha256((HERE/'repeated_inputs.npz').read_bytes()).hexdigest()}
    (HERE/'input_reconstruction_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
