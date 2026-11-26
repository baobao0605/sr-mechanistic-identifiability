"""Record-disjoint 6+6 transfer of empirical SR margin scores.

The first six sessions train the score calibration and abstention rules.
Only after those rules are written/frozen are sessions 7--12 analyzed.
This is a retrospective within-participant transfer analysis, not a new
participant sample or a confidence-coverage guarantee.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score,brier_score_loss
ROOT=Path(__file__).resolve().parents[1];HERE=Path(__file__).resolve().parent
OUT=HERE/'session_disjoint';OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'theory'))
from continuous_geometry import SpectralSRGeometry
from analyze_margin_inputs import Orbit,multiplier,apply_coeff,KINDS,GOLD,STATES,CHANNELS,TG

def fit_all(b,x,y):
    xx=x.reshape(-1,len(b));bb=b@b;bx=xx@b/bb;by=y@b/bb
    r=xx-bx[:,None]*b;yr=y-by[:,None]*b;rr=(r*r).sum(1)
    dot=r@yr.T
    alpha=np.maximum(0,np.divide(dot,rr[:,None],out=np.zeros_like(dot),where=rr[:,None]>1e-20))
    c=by[None,:]-alpha*bx[:,None]
    sse=np.maximum(0,(yr*yr).sum(1)[None,:]-2*alpha*dot+alpha*alpha*rr[:,None])
    shape=(len(KINDS),len(GOLD),len(y))
    return sse.reshape(shape),c.reshape(shape),alpha.reshape(shape)

def run_half(offset,label):
    z=np.load(ROOT/'prior/restored/complete/analysis/D05/analysis_input.npz')
    data0=z['patterns'];roi=z['voxel_roi'];rows=[];scores=[];geometry=[];splits=[]
    for n in STATES:
      coeffs=np.array([multiplier(n,k,GOLD) for k in KINDS])
      for nc in CHANNELS:
        ix=np.concatenate([np.where(roi==r)[0][:nc//3] for r in ('V1','lFFA','rFFA')])
        raw=data0[:,:n,ix];raw=raw-raw.mean(1,keepdims=True)
        wi=np.array([offset,offset+4,offset+5]);rr=raw[wi]-raw[wi].mean(0)
        cov=LedoitWolf().fit(rr.reshape(-1,nc));ev,V=np.linalg.eigh(cov.covariance_)
        W=(V/np.sqrt(np.maximum(ev,1e-20)))@V.T
        A,B,C,D=raw[offset:offset+4]@W
        splits.append(dict(phase=label,states=n,channels=nc,A=offset+1,B=offset+2,C=offset+3,D=offset+4,
                           W_sessions=','.join(map(str,wi+1)),whitening_shrinkage=float(cov.shrinkage_)))
        oa=Orbit(A);geom=SpectralSRGeometry(np.exp(-2j*np.pi*np.arange(n)/n),oa.w)
        lookup={}
        for kind in KINDS[:3]:
          for gamma in TG:
            margins={'current':1.}
            for rival in KINDS[:3]:
              if rival==kind:continue
              res=geom.point_margin(kind,gamma,rival,atol=1e-4,max_evals=4097)
              geometry.append(dict(phase=label,states=n,channels=nc,truth=kind,true_gamma=gamma,rival=rival,**res))
              margins[rival]=res['lower']
            lookup[(kind,gamma)]=margins
        XA=apply_coeff(coeffs,A);XC=apply_coeff(coeffs,C)
        truths=[('current',0.,0.,0.,np.ones(n))]
        for kind in KINDS[:3]:
          for gamma in TG:
            norm=geom.residual_norm(kind,gamma)
            for rho in (.1,.3,1.):
              alpha=rho*oa.norm/norm
              truths.append((kind,gamma,rho,alpha,1+alpha*multiplier(n,kind,[gamma])[0]))
        truth_co=np.array([t[-1] for t in truths]);ytr=apply_coeff(truth_co,B);yte=apply_coeff(truth_co,D)
        muA=apply_coeff(truth_co,A);E=ytr-muA
        e_full=np.linalg.norm(E,axis=(-2,-1));e_u=oa.norm_projected(E)
        sse,cs,alphas=fit_all(A.ravel(),XA,ytr.reshape(len(truths),-1));gi=sse.argmin(1)
        ct=np.einsum('tnp,np->t',ytr,A)/np.sum(A*A)
        pred0=ct[:,None,None]*C;err0=np.sum((pred0-yte)**2,axis=(-2,-1))
        errs=[err0];fitinfo=[];predtrain=[]
        for ki,kind in enumerate(KINDS):
          ar=np.arange(len(truths));ixg=gi[ki];c=cs[ki,ixg,ar];a=alphas[ki,ixg,ar]
          pred=c[:,None,None]*C+a[:,None,None]*XC[ki,ixg]
          predA=c[:,None,None]*A+a[:,None,None]*XA[ki,ixg]
          errs.append(np.sum((pred-yte)**2,axis=(-2,-1)))
          norm=np.array([geom.residual_norm('symmetric_transition' if kind=='balanced_fb' else kind,GOLD[g]) for g in ixg])
          # balanced_fb is proportional to symmetric after current removal.
          if kind=='balanced_fb':
            g=GOLD[ixg]/(1+np.sqrt(1-GOLD[ixg]**2));norm*=((1-g*g)/(1+g*g)**2)
          active=a*norm/oa.norm>=1e-10
          fitinfo.append((ixg,c,a,active));predtrain.append(predA)
        errs=np.array(errs);den=np.sum(yte*yte,axis=(-2,-1))
        for t,(kind,gamma,rho,alpha,co) in enumerate(truths):
          tol=1e-10*max(den[t],1e-30);best=errs[:,t].min();tied=np.where(errs[:,t]<=best+tol)[0]
          classes=[]
          for q in tied:
            classes.append('current' if q==0 or not fitinfo[q-1][3][t] else 'symmetric' if q in (3,4) else KINDS[q-1])
          classes=sorted(set(classes));truthclass='symmetric' if kind=='symmetric_transition' else kind
          bb=int(np.argmin(errs[1:3,t]));bl=('forward','backward')[bb]
          unique=abs(errs[1,t]-errs[2,t])>tol
          row=dict(phase=label,states=n,channels=nc,truth=kind,true_gamma=gamma,predictive_strength_ratio=rho,
                   true_gain=alpha,winner_classes='|'.join(classes),correct_family=(classes==[truthclass]),
                   binary_forced_correct=(bl==kind),binary_visible_unique_correct=(unique and bl==kind and bool(fitinfo[bb][3][t])),
                   current_relative_sse=float(err0[t]/den[t]),repeat_error_train_full=float(e_full[t]),
                   repeat_error_train_krylov=float(e_u[t]),norm_reference_train=oa.norm)
          if kind!='current':
            qtrain=float(np.linalg.norm(predtrain[KINDS.index(kind)][t]-muA[t]));row['q_train']=qtrain
            margins=lookup[(kind,gamma)]
            for scope,rv in [('binary',['backward' if kind=='forward' else 'forward'] if kind in ('forward','backward') else []),('multi',[k for k in KINDS[:3] if k!=kind]+['current'])]:
              if not rv:continue
              unit=min(margins[k] for k in rv);m=rho*oa.norm*unit
              row['unit_geometry_'+scope]=unit;row['margin_'+scope+'_lower']=m
              for sp,e in [('full',e_full[t]),('krylov',e_u[t])]:
                row['I_pair_q_'+scope+'_'+sp]=m/max(2*e+qtrain,1e-30)
                row['I_pair_'+scope+'_'+sp]=m/max(2*e,1e-30)
                row['I_naive_'+scope+'_'+sp]=m/max(np.sqrt(2)*e,1e-30)
          rows.append(row)
          for ki,candidate in enumerate(KINDS):
            gg,c,a,active=fitinfo[ki]
            scores.append(dict(phase=label,states=n,channels=nc,truth=kind,true_gamma=gamma,predictive_strength_ratio=rho,
                               candidate=candidate,selected_grid_gamma=float(GOLD[gg[t]]),current_gain=float(c[t]),predictive_gain=float(a[t]),
                               active=bool(active[t]),train_relative_sse=float(sse[ki,gg[t],t]/np.sum(ytr[t]**2)),test_relative_sse=float(errs[ki+1,t]/den[t])))
    d=pd.DataFrame(rows);d.to_csv(OUT/(label+'_trials.csv'),index=False)
    pd.DataFrame(scores).to_csv(OUT/(label+'_model_scores.csv'),index=False)
    pd.DataFrame(geometry).to_csv(OUT/(label+'_geometry.csv'),index=False)
    pd.DataFrame(splits).to_csv(OUT/(label+'_splits.csv'),index=False)
    return d

def analyses(d):
    fb=d[d.truth.isin(['forward','backward'])];pred=d[d.truth!='current']
    for cohort,scope,outcome,frame in [('FB','binary','binary_visible_unique_correct',fb),('FB','multi','correct_family',fb),('all_predictive','multi','correct_family',pred)]:
      cols={'I_pair_q_krylov':'I_pair_q_'+scope+'_krylov','I_pair_q_full':'I_pair_q_'+scope+'_full',
            'I_pair_krylov':'I_pair_'+scope+'_krylov','I_naive_krylov':'I_naive_'+scope+'_krylov',
            'unit_geometry':'unit_geometry_'+scope,'strength_only':'predictive_strength_ratio'}
      yield cohort,scope,outcome,frame,cols

def fit_rules(d):
    rules={}
    for cohort,scope,outcome,frame,cols in analyses(d):
      y=frame[outcome].to_numpy(bool)
      for name,col in cols.items():
        x=frame[col].to_numpy(float);iso=IsotonicRegression(increasing=True,out_of_bounds='clip').fit(x,y)
        thresholds=[]
        for risk in (.1,.2):
          best=None
          for t in np.unique(x):
            choose=x>=t;n=int(choose.sum());wrong=int((~y[choose]).sum())
            if n<20 or wrong/n>risk+1e-12:continue
            candidate={'threshold':float(t),'calibration_selected':n,'calibration_wrong':wrong,'calibration_error':wrong/n}
            if best is None or (n,float(t))>(best['calibration_selected'],best['threshold']):best=candidate
          thresholds.append({'target_empirical_error':risk,'status':'selected' if best else 'abstain_all',**(best or {'threshold':None,'calibration_selected':0,'calibration_wrong':0,'calibration_error':None})})
        rules[cohort+'__'+scope+'__'+name]={'cohort':cohort,'scope':scope,'outcome':outcome,'predictor':name,'column':col,
                                           'isotonic_x':iso.X_thresholds_.tolist(),'isotonic_y':iso.y_thresholds_.tolist(),
                                           'calibration_n':len(x),'calibration_accuracy':float(y.mean()),
                                           'risk_thresholds':thresholds,'calibration_raw_auc':float(roc_auc_score(y,x))}
    payload={'protocol_version':'v2','protocol_sha256':hashlib.sha256((HERE/'session_disjoint_contract_v2.md').read_bytes()).hexdigest(),
             'calibration_sessions':[1,2,3,4,5,6],'validation_sessions_not_yet_computed':[7,8,9,10,11,12],
             'primary_predictor':'I_pair_q_krylov','rule_count':len(rules),'rules':rules}
    p=OUT/'calibration_rules.json';p.write_text(json.dumps(payload,ensure_ascii=False,indent=2))
    return payload,hashlib.sha256(p.read_bytes()).hexdigest()

def evaluate(d,payload,sha):
    frames={(c,s):(f,o) for c,s,o,f,cols in analyses(d)};rows=[];threshold=[];predrows=[]
    for key,r in payload['rules'].items():
      f,out=frames[(r['cohort'],r['scope'])];x=f[r['column']].to_numpy(float);y=f[out].to_numpy(bool)
      yp=np.interp(x,np.asarray(r['isotonic_x']),np.asarray(r['isotonic_y']))
      rows.append(dict(cohort=r['cohort'],scope=r['scope'],predictor=r['predictor'],n=len(x),
                       calibration_accuracy=r['calibration_accuracy'],validation_accuracy=float(y.mean()),
                       calibration_auc=r['calibration_raw_auc'],validation_raw_auc=float(roc_auc_score(y,x)),
                       validation_calibrated_auc=float(roc_auc_score(y,yp)),mean_predicted_success=float(yp.mean()),
                       brier=float(brier_score_loss(y,yp)),constant_brier=float(brier_score_loss(y,np.full(len(y),r['calibration_accuracy'])))))
      for rule in r['risk_thresholds']:
        sel=np.zeros(len(x),bool) if rule['threshold'] is None else x>=rule['threshold']
        threshold.append(dict(cohort=r['cohort'],scope=r['scope'],predictor=r['predictor'],**rule,
                              validation_n=len(x),validation_selected=int(sel.sum()),validation_wrong=int((~y[sel]).sum()),
                              validation_error=float((~y[sel]).mean()) if sel.any() else np.nan,
                              validation_coverage=float(sel.mean())))
      for (_,rr),score,prob,outcome in zip(f.iterrows(),x,yp,y):
        predrows.append(dict(cohort=r['cohort'],scope=r['scope'],predictor=r['predictor'],states=rr.states,channels=rr.channels,
                             truth=rr.truth,true_gamma=rr.true_gamma,predictive_strength_ratio=rr.predictive_strength_ratio,
                             score=score,calibrated_success_probability=prob,correct=bool(outcome)))
    pd.DataFrame(rows).to_csv(OUT/'validation_performance.csv',index=False)
    pd.DataFrame(threshold).to_csv(OUT/'validation_thresholds.csv',index=False)
    pd.DataFrame(predrows).to_csv(OUT/'validation_predictions.csv',index=False)
    assert hashlib.sha256((OUT/'calibration_rules.json').read_bytes()).hexdigest()==sha
    return rows,threshold

def main():
    print('Computing calibration sessions 1--6',flush=True)
    cal=run_half(0,'calibration');rules,sha=fit_rules(cal)
    print('Calibration rules frozen',sha,flush=True)
    print('Computing validation sessions 7--12',flush=True)
    val=run_half(6,'validation');performance,thresholds=evaluate(val,rules,sha)
    summary={'calibration_n':len(cal),'validation_n':len(val),'sessions_disjoint':True,'rules_sha256':sha,
             'primary_predictor':'I_pair_q_krylov',
             'primary_performance':[x for x in performance if x['predictor']=='I_pair_q_krylov'],
             'primary_thresholds':[x for x in thresholds if x['predictor']=='I_pair_q_krylov'],
             'current_controls':{'calibration_correct':int(cal.loc[cal.truth=='current','correct_family'].sum()),
                                 'validation_correct':int(val.loc[val.truth=='current','correct_family'].sum()),'n_each':30},
             'scope':'Retrospective one-participant record-disjoint transfer. Same images and nested computational conditions; no participant-level confidence or FDR claims.'}
    (OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
