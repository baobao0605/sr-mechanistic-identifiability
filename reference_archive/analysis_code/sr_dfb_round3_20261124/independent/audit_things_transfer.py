"""Independent audit of the THINGS record-disjoint transfer implementation.

No rule is selected here. Recompute frozen rules from calibration rows,
validate recorded splits/hashes, and use direct real matrices + an SVD
projector to check selected raw-record -> index calculations.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'; os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.covariance import LedoitWolf

ROOT=Path(__file__).resolve().parents[1]; HERE=Path(__file__).resolve().parent
EMP=ROOT/'empirical'; DATA=EMP/'session_disjoint'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def kernel(n,kind,gamma):
    P=np.zeros((n,n));P[np.arange(n),(np.arange(n)-1)%n]=1
    if kind=='forward':K=P
    elif kind=='backward':K=P.T
    elif kind=='symmetric_transition':K=(P+P.T)/2
    elif kind=='balanced_fb':
        g=gamma/(1+np.sqrt(1-gamma*gamma))
        return (kernel(n,'forward',g)+kernel(n,'backward',g))/2
    else:raise ValueError(kind)
    return np.linalg.solve(np.eye(n)-gamma*K,K)

def main():
    rules=json.loads((DATA/'calibration_rules.json').read_text())
    summary=json.loads((DATA/'summary.json').read_text())
    cal=pd.read_csv(DATA/'calibration_trials.csv')
    val=pd.read_csv(DATA/'validation_trials.csv')
    assert rules['protocol_sha256']==sha(EMP/'session_disjoint_contract_v2.md')
    assert summary['rules_sha256']==sha(DATA/'calibration_rules.json')
    result={'status':'PASS_WITH_SCOPE','rules_sha256':sha(DATA/'calibration_rules.json'),
            'protocol_sha256':rules['protocol_sha256'],'code_sha256':sha(EMP/'session_disjoint.py'),
            'rules_checked':0,'threshold_rules_checked':0,'isotonic_max_error':0.,
            'record_groups_disjoint':True,'primary_index_max_relative_error':0.,
            'raw_input_conditions_checked':0,'transport_full_max_error':0.,
            'transport_projected_max_error':0.,'q_train_max_error':0.,'matrix_fft_orientation_max_error':0.}
    for name,d,offset in [('calibration',cal,0),('validation',val,6)]:
        sp=pd.read_csv(DATA/(name+'_splits.csv'))
        assert len(d)==840 and len(sp)==30
        assert all(sp.A.eq(offset+1)) and all(sp.B.eq(offset+2))
        assert all(sp.C.eq(offset+3)) and all(sp.D.eq(offset+4))
        assert all(sp.W_sessions.eq(','.join(map(str,[offset+1,offset+5,offset+6]))))
    a=set(range(1,7));b=set(range(7,13));assert a.isdisjoint(b)
    for r in rules['rules'].values():
        d=cal[cal.truth.isin(['forward','backward'])] if r['cohort']=='FB' else cal[cal.truth.ne('current')]
        x=d[r['column']].to_numpy(float);y=d[r['outcome']].to_numpy(bool)
        fit=IsotonicRegression(increasing=True,out_of_bounds='clip').fit(x,y)
        assert len(fit.X_thresholds_)==len(r['isotonic_x'])
        err=max(np.max(np.abs(fit.X_thresholds_-r['isotonic_x'])),np.max(np.abs(fit.y_thresholds_-r['isotonic_y'])))
        assert err<1e-12;result['isotonic_max_error']=max(result['isotonic_max_error'],float(err))
        for frozen in r['risk_thresholds']:
            candidates=[]
            for t in np.unique(x):
                keep=x>=t;count=int(keep.sum());wrong=int((~y[keep]).sum())
                if count>=20 and wrong/count<=frozen['target_empirical_error']+1e-12:
                    candidates.append((count,float(t),wrong))
            if candidates:
                count,t,wrong=max(candidates,key=lambda z:(z[0],z[1]))
                assert frozen['status']=='selected' and abs(t-frozen['threshold'])<1e-12
                assert count==frozen['calibration_selected'] and wrong==frozen['calibration_wrong']
            else:assert frozen['status']=='abstain_all' and frozen['threshold'] is None
            result['threshold_rules_checked']+=1
        result['rules_checked']+=1
    # Real-space evaluation is independent of the empirical FFT Orbit and
    # the spectral geometry implementation. The cyclic matrix orientation
    # must agree with that already fixed in the source task, as audited below.
    z=np.load(ROOT/'prior/restored/complete/analysis/D05/analysis_input.npz')
    raw0=z['patterns'];roi=z['voxel_roi'];checks=[]
    cache={}
    for phase,d,offset in [('calibration',cal,0),('validation',val,6)]:
      scores=pd.read_csv(DATA/(phase+'_model_scores.csv'))
      for n in (4,16,100):
       for nc in (6,72):
        ix=np.concatenate([np.flatnonzero(roi==r)[:nc//3] for r in ('V1','lFFA','rFFA')])
        raw=raw0[:,:n,ix];raw=raw-raw.mean(1,keepdims=True)
        wraw=raw[[offset,offset+4,offset+5]];res=wraw-wraw.mean(0)
        cov=LedoitWolf().fit(res.reshape(-1,nc)).covariance_
        ev,V=np.linalg.eigh(cov);W=(V/np.sqrt(np.maximum(ev,1e-20)))@V.T
        A=raw[offset]@W;B=raw[offset+1]@W
        P=np.zeros((n,n));P[np.arange(n),(np.arange(n)-1)%n]=1
        lam=np.exp(-2j*np.pi*np.arange(n)/n)
        fft_forward=np.fft.ifft((lam/(1-.5*lam))[:,None]*np.fft.fft(A,axis=0),axis=0).real
        orientation_error=float(np.linalg.norm(kernel(n,'forward',.5)@A-fft_forward))
        result['matrix_fft_orientation_max_error']=max(result['matrix_fft_orientation_max_error'],orientation_error)
        blocks=[];v=A.copy()
        for j in range(n):blocks.append(v.ravel());v=P@v
        left,s,_=np.linalg.svd(np.column_stack(blocks),full_matrices=False)
        U=left[:,s>np.max(s)*1e-11]
        subset=d[d.states.eq(n)&d.channels.eq(nc)&d.truth.ne('current')]
        for row in subset.itertuples(index=False):
            def K(kind,g):
                key=(n,kind,float(g))
                if key not in cache:cache[key]=kernel(*key)
                return cache[key]
            T=np.eye(n)+row.true_gain*K(row.truth,row.true_gamma)
            E=T@(B-A);e=float(np.linalg.norm(E));eu=float(np.linalg.norm(U.T@E.ravel()))
            r=scores[scores.states.eq(n)&scores.channels.eq(nc)&scores.truth.eq(row.truth)&scores.true_gamma.eq(row.true_gamma)&scores.predictive_strength_ratio.eq(row.predictive_strength_ratio)&scores.candidate.eq(row.truth)]
            assert len(r)==1;f=r.iloc[0]
            pred=f.current_gain*A+f.predictive_gain*K(row.truth,f.selected_grid_gamma)@A
            q=float(np.linalg.norm(pred-T@A))
            result['transport_full_max_error']=max(result['transport_full_max_error'],abs(e-row.repeat_error_train_full))
            result['transport_projected_max_error']=max(result['transport_projected_max_error'],abs(eu-row.repeat_error_train_krylov))
            result['q_train_max_error']=max(result['q_train_max_error'],abs(q-row.q_train))
            for mode in (('binary','multi') if row.truth in ('forward','backward') else ('multi',)):
                m=getattr(row,'margin_'+mode+'_lower');index=m/(2*eu+q)
                saved=getattr(row,'I_pair_q_'+mode+'_krylov')
                re=abs(index-saved)/max(1.,abs(saved))
                result['primary_index_max_relative_error']=max(result['primary_index_max_relative_error'],re)
                checks.append({'phase':phase,'states':n,'channels':nc,'truth':row.truth,'gamma':row.true_gamma,
                               'rho':row.predictive_strength_ratio,'mode':mode,'index':index,'saved_index':saved,'difference':index-saved})
            result['raw_input_conditions_checked']+=1
    assert max(result[k] for k in ('transport_full_max_error','transport_projected_max_error','q_train_max_error'))<1e-7
    assert result['primary_index_max_relative_error']<1e-10
    assert result['matrix_fft_orientation_max_error']<1e-8
    result['scope']='Retrospective same-person/same-images transfer; record independence is a partition fact, not an inferential independence guarantee. Freeze order verified in code and run log; SHA verifies unchanged saved rule, not an external pre-registration timestamp.'
    pd.DataFrame(checks).to_csv(HERE/'review_things_raw_checks.csv',index=False)
    (HERE/'review_things_numeric_audit.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
