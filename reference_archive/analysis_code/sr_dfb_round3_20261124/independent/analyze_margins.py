"""Point-to-SR-family margins and real-repeat errors: Kern/Turner.

No new model fits to held-out targets are used for prediction. Post-hoc
oracle error decompositions are explicitly separated from pretest indices.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import pandas as pd
from scipy.linalg import norm

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
OLD=ROOT/'prior/sr_round2_20261124/math'
sys.path.insert(0,str(ROOT/'theory'))
from continuous_geometry import SpectralSRGeometry
from build_inputs import kernel

KEYS=['dataset','window','fold','space']
TRUTH=['true_kind','true_gamma','prediction_residual_strength']
GAMMAS=(.2,.5,.8)

def fnorm(x):return float(norm(np.asarray(x).ravel()))

def krylov_basis(A,B):
    cols=[];x=B.copy()
    for _ in range(len(A)):
        cols.append(x.ravel());x=A@x
    D=np.array(cols).T
    U,s,_=np.linalg.svd(D,full_matrices=False)
    keep=s>max(D.shape)*np.finfo(float).eps*s[0]
    return U[:,keep],s

def residual(B,X):
    return X-B*(np.sum(B*X)/np.sum(B*B))

def winner(rows,candidates):
    r=rows[rows.candidate.isin(candidates)]
    best=float(r.test_sse.min());eps=1e-10*float(r.test_target_norm2.iloc[0])
    tied=r.loc[r.test_sse<=best+eps,'candidate'].tolist()
    tied=sorted(set('symmetric' if x in ('balanced_fb','symmetric_transition') else x for x in tied))
    return '|'.join(tied)

def margin_call(geo,kind,g,rival,intervals,atol,max_evals):
    if rival=='current':return {'lower':1.,'upper':1.,'status':'exact_current','evals':0}
    out=geo.point_margin(kind,g,rival,intervals=intervals,atol=atol,max_evals=max_evals)
    if not isinstance(out,dict):out=vars(out)
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--max-folds',type=int);ap.add_argument('--atol',type=float,default=2e-4);ap.add_argument('--max-evals',type=int,default=2049)
    args=ap.parse_args();cache=np.load(HERE/'repeated_inputs.npz')
    meta=pd.read_csv(HERE/'input_folds.csv')
    if args.max_folds:meta=meta.iloc[:args.max_folds]
    scores=pd.read_csv(OLD/'calibration_folds.csv')
    rows=[];geometry=[];scorechecks=[];started=time.time()
    for mi,m in meta.iterrows():
        key={k:m[k] for k in KEYS};fid=m.fold_id
        A,R,V,S,Z=[cache[f'{fid}_{k}'] for k in ('A','R','V','S','Z')]
        UR,sr=krylov_basis(A,R);US,ss=krylov_basis(A,S)
        geos={'train':SpectralSRGeometry.from_reference(A,R),'test':SpectralSRGeometry.from_reference(A,S)}
        margins={}
        for phase,geo in geos.items():
            for truekind in ('forward','backward'):
                other='backward' if truekind=='forward' else 'forward'
                for g in GAMMAS:
                    for rival in (other,'symmetric_transition','current','gamma_far'):
                        intervals=[(0.,.95)]
                        kind=rival
                        if rival=='gamma_far':
                            kind=truekind;intervals=[]
                            if g-.2>=-1e-10:intervals.append((0.,max(0.,g-.2)))
                            if g+.2<=.95+1e-10:intervals.append((min(.95,g+.2),.95))
                        mr=margin_call(geo,truekind,g,kind,intervals,args.atol,args.max_evals)
                        margins[phase,truekind,g,rival]=mr
                        geometry.append({**key,'fold_id':fid,'phase':phase,'true_kind':truekind,'true_gamma':g,'rival':rival,'lower':mr['lower'],'upper':mr['upper'],'interval_width':mr['upper']-mr['lower'],'status':mr.get('status','unspecified'),'evals':mr.get('evaluations',mr.get('evals',np.nan)),'gamma_rival':mr.get('gamma_witness',np.nan)})
        part=scores.copy()
        for k in KEYS:part=part[part[k].eq(key[k])]
        for tkey,sub in part.groupby(TRUTH,sort=False):
            truekind,g0,rho=tkey;g0=float(g0);rho=float(rho)
            a0=float(sub.true_alpha.iloc[0]);T0=np.eye(len(A))+a0*kernel(A,g0,truekind)
            muR=T0@R;muS=T0@S;Ytr=T0@V;Yte=T0@Z
            e_train=(Ytr-muR).ravel();e_ref=(T0@(S-R)).ravel()/np.sqrt(2);e_test=(Yte-muS).ravel()
            pred={};qR=qS=rtrue=np.nan
            projectY=(US@(US.T@Yte.ravel())).reshape(Yte.shape)
            spectral_projection_error=fnorm(projectY-geos['test'].project_krylov(Yte))
            discarded=fnorm(Yte-projectY)**2
            for idx,fit in sub.iterrows():
                if fit.candidate=='current':Tf=float(fit.current_gain)*np.eye(len(A))
                else:Tf=float(fit.current_gain)*np.eye(len(A))+float(fit.predictive_gain)*kernel(A,float(fit.gamma_grid_argmin_audit_only),fit.candidate)
                pred[fit.candidate]=Tf@S
                dfit=fnorm(Tf@S-(US@(US.T@(Tf@S).ravel())).reshape(S.shape))
                sseproj=fnorm(projectY-Tf@S)**2
                scorechecks.append({**key,**dict(zip(TRUTH,tkey)),'candidate':fit.candidate,'old_sse':float(fit.test_sse),'projected_sse':sseproj,'discarded_energy':discarded,'sse_identity_error':float(fit.test_sse)-sseproj-discarded,'prediction_outside_U_norm':dfit})
                if fit.candidate==truekind:
                    qR=fnorm(Tf@R-muR);qS=fnorm(Tf@S-muS);rtrue=fnorm(projectY-Tf@S)
                    ghat=float(fit.gamma_hat) if np.isfinite(fit.gamma_hat) else np.nan
            row={**key,'fold_id':fid,**dict(zip(TRUTH,tkey)),'true_alpha':a0,
                 'krylov_rank_train':UR.shape[1],'krylov_rank_test':US.shape[1],
                 'n_observations':R.size,'train_reference_norm':fnorm(R),'test_reference_norm':fnorm(S),
                 'spectral_vs_svd_projection_error':spectral_projection_error,
                 'train_current_residual_strength':fnorm(residual(R,muR)),
                 'test_current_residual_strength':fnorm(residual(S,muS)),
                 'train_transport_error_full':fnorm(e_train),'train_transport_error_U':fnorm(UR.T@e_train),
                 'reference_repeat_error_full':fnorm(e_ref),'reference_repeat_error_U':fnorm(UR.T@e_ref),
                 'oracle_test_error_full':fnorm(e_test),'oracle_test_error_U':fnorm(US.T@e_test),
                 'correct_fit_error_train':qR,'correct_fit_error_test':qS,'correct_fit_test_residual_U':rtrue,
                 'binary_winners':winner(sub,['forward','backward']),
                 'multi_winners':winner(sub,['forward','backward','symmetric_transition','current']),
                 'gamma_hat_true_family':ghat,
                 'gamma_reportable':bool(np.isfinite(ghat)),
                 'gamma_error_lt_point2':bool(np.isfinite(ghat) and abs(ghat-g0)<.2-1e-10)}
            row['binary_correct']=truekind!='current' and row['binary_winners']==truekind
            row['multi_correct']=row['multi_winners']==truekind
            if truekind!='current':
                other='backward' if truekind=='forward' else 'forward'
                for phase in ('train','test'):
                    strength=row[f'{phase}_current_residual_strength']
                    for label,rival in [('binary',other),('symmetric','symmetric_transition'),('current','current'),('gamma_far','gamma_far')]:
                        out=margins[phase,truekind,g0,rival]
                        for bound in ('lower','upper'):row[f'{phase}_{label}_margin_{bound}']=strength*float(out[bound])
                    for bound in ('lower','upper'):
                        row[f'{phase}_multi_margin_{bound}']=min(row[f'{phase}_{j}_margin_{bound}'] for j in ('binary','symmetric','current'))
                for mode in ('binary','multi'):
                    dm=row[f'test_{mode}_margin_lower']
                    row[f'oracle_{mode}_full_bound_pass']=bool(dm>2*row['oracle_test_error_full']+qS+1e-10)
                    row[f'oracle_{mode}_U_bound_pass']=bool(dm>2*row['oracle_test_error_U']+qS+1e-10)
                    row[f'oracle_{mode}_U_tight_pass']=bool(dm>row['oracle_test_error_U']+rtrue+1e-10)
                    for scale in ('transport','repeat'):
                        ev=row['train_transport_error_U' if scale=='transport' else 'reference_repeat_error_U']
                        row[f'pretest_{mode}_{scale}_ratio']=row[f'train_{mode}_margin_lower']/max(2*ev+qR,np.finfo(float).tiny)
                    row[f'pretest_{mode}_unit_geometry']=row[f'train_{mode}_margin_lower']/max(row['train_current_residual_strength'],np.finfo(float).tiny)
            rows.append(row)
        pd.DataFrame(rows).to_csv(HERE/'margin_folds.csv',index=False)
        pd.DataFrame(geometry).to_csv(HERE/'continuous_geometry.csv',index=False)
        pd.DataFrame(scorechecks).to_csv(HERE/'projection_checks.csv',index=False)
        print(f'Completed {mi+1}/{len(meta)} {key}; {time.time()-started:.1f}s',flush=True)
    checks=pd.DataFrame(scorechecks);df=pd.DataFrame(rows)
    maxerr=float(checks.sse_identity_error.abs().max());outside=float(checks.prediction_outside_U_norm.max())
    assert maxerr<1e-8
    for mode in ('binary','multi'):
        for bound in ('full_bound','U_bound','U_tight'):
            mask=df.get(f'oracle_{mode}_{bound}_pass',pd.Series(False,index=df.index)).fillna(False)
            assert bool(df.loc[mask,f'{mode}_correct'].all()),f'failed oracle implication {mode} {bound}'
    summary={'folds':len(meta),'conditions':len(df),'geometry_evaluations':len(geometry),'score_checks':len(checks),'projection_sse_identity_max_abs':maxerr,'prediction_outside_U_max_norm':outside,'seconds':time.time()-started,'atol':args.atol,'max_evals':args.max_evals,'oracle_implications':'PASS','status':'COMPLETE'}
    (HERE/'margin_run_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
