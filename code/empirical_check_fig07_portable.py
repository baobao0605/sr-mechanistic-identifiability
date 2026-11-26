"""Portable Figure 7 source checks and refits from the included measured excerpt.

Run from any directory with Python, NumPy and pandas. No AB archive is needed.
The excerpt is measured, baseline-corrected EEG from the prior supplied caches;
this is not a raw-EEG preprocessing pipeline. The reference CV refit below follows
the explicit leave-event-out independent audit already present in the AB release;
re-executing it is not claimed as a third independent implementation.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]

def center(x):return np.asarray(x)-np.mean(x,axis=0,keepdims=True)

def contrast(n,k=1):
    t=2*np.pi*np.arange(n)/n
    return np.array([g(j*t) for j in range(1,k+1) for g in (np.cos,np.sin)])*np.sqrt(2/n)

def refit(pools):
    n=len(pools);p=pools[0].shape[-1];m=min(map(len,pools))
    ranks=sorted(set(min(n-1,p,r) for r in (1,2,4,8,16,32,64,p)))
    smoothing={'fourier_'+str(k):contrast(n,k).T@contrast(n,k) for k in (1,2,4,8,19) if 2*k<n}
    smoothing['fourier_full']=np.eye(n)-np.ones((n,n))/n
    candidates=[('scalar',0)]+[('svd',r) for r in ranks]+[(name,0) for name in smoothing]
    indices=[np.rint(np.linspace(0,len(a)-1,m)).astype(int) for a in pools]
    folds=[]
    for j in range(m):
        tr=center(np.array([np.delete(a,indices[s][j],axis=0).mean(0) for s,a in enumerate(pools)]))
        va=center(np.array([a[indices[s][j]] for s,a in enumerate(pools)]))
        u,sv,v=np.linalg.svd(tr,full_matrices=False)
        bases={('scalar',0):tr}
        bases.update({('svd',r):(u[:,:r]*sv[:r])@v[:r] for r in ranks})
        bases.update({(name,0):op@tr for name,op in smoothing.items()})
        folds.append((bases,va))
    scores={(family,r,float(lam)):np.mean([np.mean((lam*b[family,r]-va)**2) for b,va in folds])
            for family,r in candidates for lam in np.linspace(0,1,11)}
    mean=center(np.array([a.mean(0) for a in pools]));u,sv,v=np.linalg.svd(mean,full_matrices=False)
    out={'zero':np.zeros_like(mean),'mean':mean}
    for name,keep in [('cv_scalar',lambda k:k[0]=='scalar'),('cv_svd',lambda k:k[0]=='svd'),('cv_fourier',lambda k:k[0].startswith('fourier'))]:
        best=min([k for k in scores if keep(k)],key=lambda k:(scores[k],k[2],k[1],k[0]))
        fam,r,lam=best
        basis=mean if fam=='scalar' else (u[:,:r]*sv[:r])@v[:r] if fam=='svd' else smoothing[fam]@mean
        out[name]=lam*basis
    return out

def main():
    sd=ROOT/'source_data';a=np.load(sd/'fig07_empirical_arrays.npz',allow_pickle=False)
    table=pd.read_csv(sd/'fig07_all_reference_estimators.csv')
    errors=[];nchecks=0
    for session in [1,2]:
        x=a[f'R{session}__measurements'];labels=a[f'R{session}__labels'];t=a[f'R{session}__time'];blocks=a[f'R{session}__block']
        for wi,window in enumerate(['early','late']):
            for b in np.unique(blocks):
                ix=blocks==b;mid=(t[ix].min()+t[ix].max())/2
                tr=ix&(t<=mid-.5);ev=ix&(t>=mid+.5)
                assert not np.any(tr&ev) and t[ev].min()-t[tr].max()>=1
                pools=[x[tr&(labels==j),wi] for j in range(40)]
                target=center(np.array([x[ev&(labels==j),wi].mean(0) for j in range(40)]))
                key=f'Turner_s{session}_b{b}_{window}_all_local';C=contrast(40)
                for method,fit in refit(pools).items():
                    row=table[(table.key==key)&(table.estimator==method)].iloc[0]
                    matrixerr=float(np.max(abs(fit-a['reference__'+key+'__'+method])))
                    vals={'outer_mse':np.mean((fit-target)**2),'zero_mse':np.mean(target**2),
                          'plane_mse':np.mean((C@(fit-target))**2),'plane_zero_mse':np.mean((C@target)**2)}
                    metricerr=max(abs(float(row[k])-v) for k,v in vals.items())
                    assert matrixerr<1e-10 and metricerr<1e-10
                    errors.append((matrixerr,metricerr));nchecks+=1
    # Rebuild frozen linear readouts from the included original event excerpt.
    raw=a['R2__measurements'][:,0];labels=a['R2__labels'];nr=int(np.bincount(labels).min())
    ix=np.array([np.flatnonzero(labels==s)[:nr] for s in range(40)]).T
    p=pd.read_csv(sd/'fig07_conditional_sr.csv');readout_errors=[]
    for context in p.context.unique():
        pref=context+'__B_trained_A__';C,W,H,R=[a[pref+k] for k in ['Ceff','W','H','R']]
        got=np.einsum('aij,rij->ra',H,np.einsum('ij,rjk,kl->ril',C,raw[ix],W))/np.linalg.norm(R)
        err=float(np.max(abs(got-a[pref+'all_readouts'])));assert err<1e-10
        assert np.array_equal(ix[:4],a[pref+'training_event_ids'])
        assert set(ix[:4].ravel()).isdisjoint(ix[4:].ravel())
        readout_errors.append(err)
    # Source-table consistency checks do not infer independent sample counts.
    b=pd.read_csv(sd/'fig07_reference_blocks.csv');g=pd.read_csv(sd/'fig07_reference_aggregate.csv')
    for w,rows in b.groupby('window'):
        value=100*(1-rows.plane_mse.mean()/rows.plane_zero_mse.mean())
        assert abs(value-g[g.window==w].gain_percent.iloc[0])<1e-10
    assert np.max(abs(p.boundary_rational/p.critical-p.ratio))<1e-12
    tr=pd.read_csv(sd/'fig07_transport_conditions.csv');summary=pd.read_csv(sd/'fig07_transport_summary.csv')
    assert len(tr)==summary.computational_conditions.sum()
    assert (tr.decision=='abstain').all()
    for dataset,rows in tr.groupby('dataset'):
        row=summary[summary.dataset==dataset].iloc[0]
        assert len(rows)==row.computational_conditions and len(rows)==row.abstentions
    result={'passed':True,'reference_refits':nchecks,'maximum_reference_matrix_error':max(x[0] for x in errors),
            'maximum_reference_metric_error':max(x[1] for x in errors),'frozen_readout_reconstructions':len(readout_errors),
            'maximum_readout_error':max(readout_errors),'cross_record_computed_conditions':len(tr),
            'scope':'Included derived measured-array refits, frozen-readout reconstruction and table consistency. Not raw EEG preprocessing, iid certification, natural SR identification or independent replication.'}
    out=ROOT/'empirical';out.mkdir(exist_ok=True)
    (out/'fig07_portable_checks.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
