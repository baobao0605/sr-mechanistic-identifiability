"""THINGS conditional model margins and actual-repeat error budgets.

Reads frozen round-2 fits; does not refit or overwrite the original study.
Grid distances are UPPER bounds for the continuous-family margin.
Training-only empirical budgets are not confidence bounds. Test-side known-
transform inequalities are retrospective calibration diagnostics.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'prior/sr_round2_20261124'
OUT=Path(__file__).resolve().parent
KINDS=('forward','backward','symmetric_transition','balanced_fb')
GOLD=np.round(np.arange(0,.951,.05),8)
GFINE=np.round(np.arange(0,.9501,.005),8)
STATES=(4,8,16,32,64,100); CHANNELS=(6,12,24,48,72); REPEATS=(1,2,3)
TG=(.25,.5,.8)
KEY=['states','channels','repeats','fold','truth','true_gamma','predictive_strength_ratio']
GROUPS=[np.arange(3*i,3*i+3) for i in range(4)]
FOLDS=[(0,1,2,3),(1,2,3,0),(2,3,0,1),(3,0,1,2)]

def multiplier(n,kind,g):
    g=np.asarray(g,float).reshape(-1,1)
    lam=np.exp(-2j*np.pi*np.arange(n)/n)[None,:]
    if kind=='forward': z=lam
    elif kind=='backward': z=lam.conjugate()
    elif kind=='symmetric_transition': z=lam.real
    elif kind=='balanced_fb':
        gb=g/(1+np.sqrt(1-g*g))
        return .5*(lam/(1-gb*lam)+lam.conjugate()/(1-gb*lam.conjugate()))
    else:raise ValueError(kind)
    return z/(1-g*z)

def apply_coeff(coeff,x):
    # coefficient: (n,) or (trials,n); x is (n,p)
    return np.fft.ifft(np.asarray(coeff)[..., :,None]*np.fft.fft(x,axis=0),axis=-2).real

class Orbit:
    def __init__(self,x):
        self.x=x;self.n=len(x);self.f=np.fft.fft(x,axis=0)
        self.energy=np.sum(abs(self.f)**2,axis=1)
        # All references were exactly state-centered before whitening.
        self.energy[0]=0.;self.f[0]=0.
        self.w=self.energy/self.n
        self.valid=self.energy>np.max(self.energy)*1e-24
        self.norm=float(np.sqrt(np.sum(self.w)))
    def norm_projected(self,e):
        ef=np.fft.fft(e,axis=-2)
        dot=np.sum(ef*self.f.conjugate(),axis=-1)
        sq=np.divide(abs(dot)**2,self.energy,out=np.zeros_like(dot.real),where=self.valid)
        return np.sqrt(np.maximum(0,np.sum(sq,axis=-1)/self.n))
    def project(self,e):
        ef=np.fft.fft(e,axis=-2)
        dot=np.sum(ef*self.f.conjugate(),axis=-1)
        coef=np.divide(dot,self.energy,out=np.zeros_like(dot),where=self.valid)
        return np.fft.ifft(coef[..., :,None]*self.f,axis=-2).real
    def rays(self,kind,gs):
        co=multiplier(self.n,kind,gs)
        mu=np.sum(co*self.w,axis=1).real/np.sum(self.w)
        rem=co-mu[:,None]
        norm=np.sqrt(np.sum(abs(rem)**2*self.w,axis=1))
        valid=norm>1e-12*np.maximum(np.sqrt(np.sum(abs(co)**2*self.w,axis=1)),1e-30)
        u=np.divide(rem*np.sqrt(self.w),norm[:,None],out=np.zeros_like(rem),where=valid[:,None])
        return u,norm,valid

def point_grid(orbit,anchor,design):
    rays={k:orbit.rays(k,GFINE) for k in KINDS[:3]}
    out=[];lookup={}
    for kind in KINDS[:3]:
      for gamma in TG:
        u,norm,valid=orbit.rays(kind,[gamma]);u=u[0];norm=float(norm[0]);valid=bool(valid[0])
        distances={}
        for rival in KINDS[:3]:
          if rival==kind:continue
          rv,rn,ok=rays[rival]
          cosine=(rv@u.conjugate()).real
          cosine=np.where(ok,cosine,-np.inf)
          index=int(np.argmax(cosine));maxcos=float(np.clip(cosine[index],-1,1))
          distance=float(np.sqrt(max(0,1-max(0,maxcos)**2))) if valid else 0.
          distances[rival]=distance
          out.append({**design,'anchor':anchor,'truth':kind,'true_gamma':gamma,'rival':rival,
                      'unit_margin_grid_upper':distance,'grid_best_gamma':float(GFINE[index]),
                      'max_cosine_grid':maxcos,'true_residual_norm':norm,'true_shape_visible':valid,
                      'grid_step':.005,'intervals':'0:0.95'})
        # Same-mechanism gamma discrimination at an a priori distance 0.2.
        rv,rn,ok=rays[kind]
        eligible=(np.abs(GFINE-gamma)>=.2-1e-12)&ok
        cosine=(rv@u.conjugate()).real
        cosine=np.where(eligible,cosine,-np.inf)
        index=int(np.argmax(cosine));maxcos=float(np.clip(cosine[index],-1,1))
        dg=float(np.sqrt(max(0,1-max(0,maxcos)**2))) if valid else 0.
        distances['gamma_0p2']=dg;distances['current']=1. if valid else 0.
        out.append({**design,'anchor':anchor,'truth':kind,'true_gamma':gamma,'rival':'gamma_0p2',
                    'unit_margin_grid_upper':dg,'grid_best_gamma':float(GFINE[index]),
                    'max_cosine_grid':maxcos,'true_residual_norm':norm,'true_shape_visible':valid,
                    'grid_step':.005,'intervals':json.dumps([[0,max(0,gamma-.2)]]+([[gamma+.2,.95]] if gamma+.2<=.95 else []))})
        distances['true_norm']=norm;lookup[(kind,gamma)]=distances
    # Uniform two-source gamma grid geometry; independent amplitudes include 0 in rival.
    family=[]
    if anchor=='train':
      for ka,kb in [('forward','backward'),('forward','symmetric_transition'),('backward','symmetric_transition')]:
        ua,_,va=rays[ka];ub,_,vb=rays[kb]
        co=(ua@ub.conjugate().T).real
        co[~(va[:,None]&vb[None,:])]=-np.inf
        ix=np.unravel_index(np.argmax(co),co.shape)
        c=float(np.clip(co[ix],-1,1));d=float(np.sqrt(max(0,1-max(c,0)**2)))
        family.append({**design,'kind_a':ka,'kind_b':kb,'unit_margin_grid_upper':d,'gamma_a_grid':GFINE[ix[0]],'gamma_b_grid':GFINE[ix[1]],'grid_step':.005})
    return out,lookup,family

def main():
    src=ROOT/'prior/restored/complete/analysis/D05/analysis_input.npz'
    z=np.load(src);data0=z['patterns'];roi=z['voxel_roi']
    old=pd.read_csv(OLD/'calibration/recovery_trials.csv')
    fits=pd.read_csv(OLD/'calibration/model_scores.csv').set_index(KEY+['candidate'])
    groups=old.groupby(['states','channels','repeats','fold'],sort=False)
    rows=[];geo=[];fam=[];cache={};checks=[]
    max_score=0.;max_current=0.;max_projection_diff=0.;max_true_span=0.;checked_models=0
    for n in STATES:
      for nc in CHANNELS:
        ix=np.concatenate([np.where(roi==s)[0][:nc//3] for s in ('V1','lFFA','rFFA')])
        data=data0[:,:n,ix];data=data-data.mean(1,keepdims=True)
        for fold,(ga,gb,gc,gd) in enumerate(FOLDS):
          rr=data[GROUPS[ga]]-data[GROUPS[ga]].mean(0)
          cov=LedoitWolf().fit(rr.reshape(-1,nc));ev,V=np.linalg.eigh(cov.covariance_)
          W=(V/np.sqrt(np.maximum(ev,1e-20)))@V.T;d=data@W
          for nr in REPEATS:
            A,B,C,D=[d[GROUPS[g][:nr]].mean(0) for g in (ga,gb,gc,gd)]
            design=dict(states=n,channels=nc,repeats=nr,fold=fold)
            cachekey=f'n{n}_p{nc}_r{nr}_f{fold}'
            for label,arr in zip('ABCD',(A,B,C,D)):cache[cachekey+'_'+label]=arr
            oa,oc=Orbit(A),Orbit(C)
            gg,la,ff=point_grid(oa,'train',design);geo+=gg;fam+=ff
            gg,lc,ff=point_grid(oc,'test_oracle',design);geo+=gg
            oldg=groups.get_group((n,nc,nr,fold))
            coeff=[]
            for t in oldg.itertuples(index=False):
              co=np.ones(n) if t.truth=='current' else 1+t.true_predictive_gain*multiplier(n,t.truth,[t.true_gamma])[0]
              coeff.append(co)
            coeff=np.array(coeff)
            ytr=apply_coeff(coeff,B);yte=apply_coeff(coeff,D)
            muA=apply_coeff(coeff,A);muC=apply_coeff(coeff,C)
            Etr=ytr-muA;Ete=yte-muC
            norm=lambda x:np.sqrt(np.sum(x*x,axis=(-2,-1)))
            eps_tr=norm(Etr);eps_te=norm(Ete)
            eps_tr_u=oa.norm_projected(Etr);eps_te_u=oc.norm_projected(Ete)
            ct=np.einsum('tnp,np->t',ytr,A)/np.sum(A*A)
            predcur=ct[:,None,None]*C
            den=np.sum(yte*yte,axis=(-2,-1))
            err_cur=norm(yte-predcur)
            max_current=max(max_current,float(np.max(abs(err_cur**2/den-oldg.current_sse.to_numpy()))))
            max_true_span=max(max_true_span,float(np.max(norm(muC-oc.project(muC)))))
            for ti,t in enumerate(oldg.itertuples(index=False)):
              key=(n,nc,nr,fold,t.truth,t.true_gamma,t.predictive_strength_ratio)
              fit_pred={};fit_active={};fit_gam={};fit_est_gam={};fit_train_pred={}
              for kind in KINDS:
                row=fits.loc[key+(kind,)]
                co=row.current_gain+row.predictive_gain*multiplier(n,kind,[row.selected_grid_gamma])[0]
                fit_pred[kind]=apply_coeff(co,C);fit_train_pred[kind]=apply_coeff(co,A)
                fit_active[kind]=bool(row.active_predictive_fit)
                fit_gam[kind]=float(row.selected_grid_gamma);fit_est_gam[kind]=float(row.estimated_gamma)
                err=float(np.sum((fit_pred[kind]-yte[ti])**2)/den[ti])
                max_score=max(max_score,abs(err-float(row.heldout_relative_sse)));checked_models+=1
              preds=[predcur[ti]]+[fit_pred[k] for k in KINDS]
              full=np.array([np.sum((p-yte[ti])**2) for p in preds])
              projected=np.array([float(oc.norm_projected(p-yte[ti]))**2 for p in preds])
              diff=(full-full[0])-(projected-projected[0])
              max_projection_diff=max(max_projection_diff,float(np.max(abs(diff))/max(den[ti],1e-30)))
              tol=1e-10*max(den[ti],1e-30)
              bbest=0 if full[1]<=full[2] else 1;blabel=('forward','backward')[bbest]
              b_unique=abs(full[1]-full[2])>tol
              bk=t.truth=='backward';opp='forward' if bk else 'backward'
              result=t._asdict()
              result.update(binary_forced_choice=blabel,binary_forced_correct=blabel==t.truth,
                            binary_unique_correct=(blabel==t.truth and b_unique),
                            binary_visible_unique_correct=(blabel==t.truth and b_unique and fit_active[blabel]),
                            norm_reference_train=oa.norm,norm_reference_test=oc.norm,
                            repeat_error_train_full=float(eps_tr[ti]),repeat_error_train_krylov=float(eps_tr_u[ti]),
                            actual_error_test_full=float(eps_te[ti]),actual_error_test_krylov=float(eps_te_u[ti]),
                            krylov_rank_train=int(oa.valid.sum()),krylov_rank_test=int(oc.valid.sum()),
                            norm_target_test=float(np.sqrt(den[ti])))
              if t.truth!='current':
                geomA=la[(t.truth,t.true_gamma)];geomC=lc[(t.truth,t.true_gamma)]
                a=t.true_predictive_gain
                ma=a*geomA['true_norm'];mc=a*geomC['true_norm']
                rivals=[k for k in KINDS[:3] if k!=t.truth]
                result.update(true_prediction_norm_train=ma,true_prediction_norm_test=mc,
                              gamma_reportable=np.isfinite(fit_est_gam[t.truth]),
                              gamma_error_conditional=abs(fit_est_gam[t.truth]-t.true_gamma),
                              gamma_correct_0p2_legacy_float=(np.isfinite(fit_est_gam[t.truth]) and abs(fit_est_gam[t.truth]-t.true_gamma)<.2),
                              gamma_correct_0p2=(np.isfinite(fit_est_gam[t.truth]) and abs(fit_est_gam[t.truth]-t.true_gamma)<.2-1e-12),
                              fit_gamma_true_family=fit_est_gam[t.truth])
                for scope,kk in [('binary',[opp] if t.truth in ('forward','backward') else []),('multi',rivals+['current']),('gamma',['gamma_0p2'])]:
                  if not kk:continue
                  da=min(geomA[k] for k in kk);dc=min(geomC[k] for k in kk)
                  result['unit_margin_train_'+scope+'_grid_upper']=da
                  result['margin_train_'+scope+'_grid_upper']=ma*da
                  result['unit_margin_test_'+scope+'_grid_upper']=dc
                  result['margin_test_'+scope+'_grid_upper']=mc*dc
                  # sqrt(2) converts a repeat difference into an empirical
                  # single-repeat RMS scale, without coverage guarantees.
                  for space,ee in [('full',eps_tr[ti]),('krylov',eps_tr_u[ti])]:
                    result['I_train_'+scope+'_'+space+'_grid_upper']=ma*da/max(np.sqrt(2)*ee,1e-30)
                  for rival in kk:
                    result['margin_train_to_'+rival+'_grid_upper']=ma*geomA[rival]
                    result['margin_test_to_'+rival+'_grid_upper']=mc*geomC[rival]
                truepred=fit_pred[t.truth];trainpred=fit_train_pred[t.truth]
                result.update(fit_error_test_true_full=float(np.linalg.norm(truepred-yte[ti])),
                              fit_error_test_true_krylov=float(oc.norm_projected(truepred-yte[ti])),
                              fit_parameter_shift_test_full=float(np.linalg.norm(truepred-muC[ti])),
                              fit_parameter_shift_test_krylov=float(oc.norm_projected(truepred-muC[ti])),
                              fit_parameter_shift_train=float(np.linalg.norm(trainpred-muA[ti])))
                for kind in KINDS[:3]:
                  result['test_error_'+kind+'_full']=float(np.linalg.norm(fit_pred[kind]-yte[ti]))
                  result['test_error_'+kind+'_krylov']=float(oc.norm_projected(fit_pred[kind]-yte[ti]))
                result['test_error_current_full']=float(err_cur[ti])
                result['test_error_current_krylov']=float(oc.norm_projected(predcur[ti]-yte[ti]))
              rows.append(result)
        print('finished',n,nc,flush=True)
    pd.DataFrame(rows).to_csv(OUT/'conditional_margin_trials_grid.csv',index=False)
    pd.DataFrame(geo).to_csv(OUT/'point_geometry_grid.csv',index=False)
    pd.DataFrame(fam).to_csv(OUT/'family_geometry_grid.csv',index=False)
    np.savez_compressed(OUT/'repeat_inputs.npz',**cache)
    audit={'n_trials':len(rows),'n_geometry_rows':len(geo),'n_family_geometry_rows':len(fam),
           'original_model_rows_checked':checked_models,'max_original_relative_sse_difference':max_score,
           'max_current_relative_sse_difference':max_current,
           'max_pairwise_SSE_difference_after_orbit_projection_relative':max_projection_diff,
           'max_true_prediction_outside_orbit':max_true_span,
           'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
           'protocol_sha256':hashlib.sha256((OUT/'analysis_contract_v1.md').read_bytes()).hexdigest()}
    assert max_score<1e-10 and max_current<1e-10 and max_projection_diff<1e-10
    df=pd.DataFrame(rows);pr=df.query('states==100 and channels==72 and repeats==3 and predictive_strength_ratio==1 and truth in ["forward","backward"]')
    audit['primary24']={'n':len(pr),'forced_binary_correct':int(pr.binary_forced_correct.sum()),'visible_unique_binary_correct':int(pr.binary_visible_unique_correct.sum()),'multi_correct':int(pr.correct_family.sum())}
    pr.to_csv(OUT/'primary24_grid.csv',index=False)
    (OUT/'frozen_fit_and_projection_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
    print(json.dumps(audit,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
