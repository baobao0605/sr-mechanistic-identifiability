"""Decision-rule ablation: nearest candidate family within the held-out target.

This is explicitly NOT held-out prediction. It aligns the comparator with
the deterministic separation theorem and quantifies an oracle bound.
"""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
OLD=ROOT/'prior/sr_round2_20261124/math'
sys.path.insert(0,str(OLD))
from sr_geometry import fit_affine,residualize_current
from build_inputs import kernel
from analyze_margins import krylov_basis,fnorm,KEYS,TRUTH

KINDS=['forward','backward','symmetric_transition']
GRID=np.linspace(0,.95,96)

def choose(records,kinds,norm2):
    subset=[r for r in records if r['candidate'] in kinds]
    best=min(r['test_profile_sse'] for r in subset);tol=1e-10*norm2
    classes=sorted(set(r['observed_class'] for r in subset if r['test_profile_sse']<=best+tol))
    return '|'.join(classes)

def main():
    cache=np.load(HERE/'repeated_inputs.npz');meta=pd.read_csv(HERE/'input_folds.csv')
    margin=pd.read_csv(HERE/'margin_folds.csv');rows=[];modelrows=[]
    for _,m in meta.iterrows():
        key={k:m[k] for k in KEYS};fid=m.fold_id
        A,S,Z=[cache[f'{fid}_{k}'] for k in ('A','S','Z')]
        U,_=krylov_basis(A,S);ss=float(np.sum(S*S))
        templ={k:np.stack([kernel(A,g,k)@S for g in GRID]) for k in KINDS}
        cur_eq={k:bool(np.linalg.norm(residualize_current(templ[k],S))<=1e-12*max(np.linalg.norm(templ[k]),1e-300)) for k in KINDS}
        part=margin[margin.fold_id.eq(fid)]
        for _,r in part.iterrows():
            tkey={k:r[k] for k in TRUTH};T0=np.eye(len(A))+float(r.true_alpha)*kernel(A,float(r.true_gamma),r.true_kind)
            Y=T0@Z;mu=T0@S;E=Y-mu
            project=lambda v:(U@(U.T@v.ravel())).reshape(v.shape)
            remove=lambda v:v-S*(float(np.sum(S*v))/ss)
            eQC=fnorm(remove(project(E)));yQC=remove(project(Y));yU=project(Y)
            norm2=float(np.sum(Y*Y));c=float(np.sum(S*Y)/ss)
            crecord={'candidate':'current','observed_class':'current','test_profile_sse':fnorm(Y-c*S)**2,'gamma_hat':np.nan,'gamma_argmin':np.nan,'current_gain':c,'predictive_gain':0.,'gamma_profile_flat':True}
            records=[crecord]
            for kind in KINDS:
                fit=fit_affine(S,templ[kind],Y);j=int(np.argmin(fit['sse']))
                tol=1e-10*max(norm2,1.)
                ties=np.asarray(fit['sse'])<=fit['sse'][j]+tol
                a=float(fit['predictive_gain'][j]);flat=bool(np.all(ties))
                obs='current' if a<=1e-12 or cur_eq[kind] else ('symmetric' if kind=='symmetric_transition' else kind)
                records.append({'candidate':kind,'observed_class':obs,'test_profile_sse':float(fit['sse'][j]),'gamma_hat':float(GRID[j]) if not flat and obs!='current' else np.nan,'gamma_argmin':float(GRID[j]),'current_gain':float(fit['current_gain'][j]),'predictive_gain':a,'gamma_profile_flat':flat})
            for rr in records:modelrows.append({**key,**tkey,'fold_id':fid,**rr})
            row={**key,**tkey,'fold_id':fid,'oracle_error_Q_U':eQC,'oracle_error_U':float(r.oracle_test_error_U),'oracle_error_full':float(r.oracle_test_error_full),
                 'binary_winners':choose(records,['forward','backward'],norm2),'multi_winners':choose(records,['forward','backward','symmetric_transition','current'],norm2)}
            row['binary_correct']=r.true_kind!='current' and row['binary_winners']==r.true_kind
            row['multi_correct']=row['multi_winners']==r.true_kind
            truth=next(v for v in records if v['candidate']==r.true_kind)
            # All candidate profiles fit the same current component. The
            # residual in Q U is its full SSE minus common orthogonal noise.
            if r.true_kind=='current':pred=truth['current_gain']*S
            else:pred=truth['current_gain']*S+truth['predictive_gain']*templ[r.true_kind][int(round(truth['gamma_argmin']*100))]
            rr=fnorm(yQC-remove(pred))
            row['correct_family_profile_residual_Q_U']=rr
            row['gamma_hat_true_family']=truth['gamma_hat'];row['gamma_reportable']=np.isfinite(truth['gamma_hat'])
            row['gamma_error_lt_point2']=bool(np.isfinite(truth['gamma_hat']) and abs(truth['gamma_hat']-r.true_gamma)<.2-1e-10)
            row['profile_sse_identity_error']=truth['test_profile_sse']-(rr**2+fnorm(Y-yU)**2)
            if r.true_kind!='current':
                for mode in ('binary','multi','gamma_far'):
                    lo=float(r[f'test_{mode}_margin_lower']);hi=float(r[f'test_{mode}_margin_upper'])
                    row[f'{mode}_margin_lower']=lo;row[f'{mode}_margin_upper']=hi
                    row[f'{mode}_oracle_Q_U_bound_pass']=lo>2*eQC+1e-10
                    row[f'{mode}_oracle_Q_U_tight_pass']=lo>eQC+rr+1e-10
                    row[f'{mode}_oracle_full_bound_pass']=lo>2*float(r.oracle_test_error_full)+1e-10
                    row[f'{mode}_oracle_U_bound_pass']=lo>2*float(r.oracle_test_error_U)+1e-10
            rows.append(row)
    d=pd.DataFrame(rows);d.to_csv(HERE/'profiled_decision_folds.csv',index=False)
    pd.DataFrame(modelrows).to_csv(HERE/'profiled_decision_models.csv',index=False)
    x=d[d.true_kind!='current'];summary={'n_fold_conditions':len(d),'n_direction_conditions':len(x),'sse_identity_max_abs':float(d.profile_sse_identity_error.abs().max()),'binary_correct':int(x.binary_correct.sum()),'multi_correct':int(x.multi_correct.sum()),'bound_results':{}}
    for mode in ('binary','multi','gamma_far'):
        truth_col='gamma_error_lt_point2' if mode=='gamma_far' else f'{mode}_correct'
        for suf in ('full_bound','U_bound','Q_U_bound','Q_U_tight'):
            col=f'{mode}_oracle_{suf}_pass';a=x[col].astype(bool)
            wrong=int((a&~x[truth_col]).sum())
            summary['bound_results'][f'{mode}_{suf}']={'selected':int(a.sum()),'wrong':wrong}
            assert wrong==0,(mode,suf)
    assert summary['sse_identity_max_abs']<1e-8
    summary['scope']='Target-level nearest-family decision; not held-out parameter transfer, and oracle errors use known calibration transform plus observed paired means.'
    (HERE/'profiled_decision_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
