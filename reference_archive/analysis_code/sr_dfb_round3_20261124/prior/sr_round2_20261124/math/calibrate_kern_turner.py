"""Real repeated-reference SR discrimination calibration, frozen protocol.json.

Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python calibrate_kern_turner.py
This program never generates people or neural noise. Known linear transforms
are applied to held-out measured means, including their measurement noise.
"""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.linalg import helmert
from sklearn.covariance import LedoitWolf
from sr_geometry import templates,residualize_current,ray_geometry,kernel_design,fit_affine

ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'restored/complete'; O=Path(__file__).resolve().parent
G=np.linspace(0,.95,96)
KINDS=['forward','backward','symmetric_transition','balanced_fb']
GRIDS={k:(G/(1+np.sqrt(1-G*G)) if k=='balanced_fb' else G) for k in KINDS}
POST=['P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4','PO8','O1','Oz','O2']
INPUTS=[]; GEOM=[]; CV=[]; PROFILE=[]; REPEAT=[]; DESIGN=[]; SPLIT=[]

def record(p):
    INPUTS.append({'path':str(p.relative_to(OLD)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})

def means(x,y,n):return np.array([x[y==i].mean(0) for i in range(n)])

def halves(x,y,n):
    ix=[[],[]]
    for j in range(n):
        ind=np.flatnonzero(y==j); q=len(ind)//2
        ix[0].extend(ind[:q]);ix[1].extend(ind[q:])
    return [(x[np.sort(i)],y[np.sort(i)],np.sort(i)) for i in ix]

def fourier(n,ks):
    t=2*np.pi*np.arange(n)/n
    return np.array([f(k*t) for k in ks for f in (np.cos,np.sin)])*np.sqrt(2/n)

def whiten_fit(x,y,n):
    b=means(x,y,n);r=x-b[y]
    cov=LedoitWolf().fit(r).covariance_
    v,U=np.linalg.eigh(cov)
    return (U/np.sqrt(np.maximum(v, max(float(v.max()),1e-300)*1e-12)))@U.T

def get_data():
    out=[]
    kp=OLD/'prior/kern/derived/shared_ssp';x={};seq=np.array([ord(s)-65 for s in 'ABCDEFGEHIBJ'])
    for s in (1,2):
        a=kp/f'localizer{s}_filtered_features.npz';b=kp/f'localizer{s}_trials.csv';record(a);record(b)
        z=np.load(a);d=pd.read_csv(b);use=d.standard.to_numpy(bool)&~d.ssp_artifact_flag.to_numpy(bool)
        for win in ('visual_early','visual_late'):
            x[s,win]=(z[win][use].astype(float),d.label.to_numpy()[use])
    out.append(('Kern',x,10,np.eye(10)[seq],
        {'full11':helmert(12,full=False),'design_k4_2D':fourier(12,[4]),'design_k4_k1_4D':fourier(12,[4,1])}))
    x={}
    for s in (1,2):
        p=OLD/f'analysis/D08/session{s}_localizer.npz';record(p);z=np.load(p)
        channels=list(z['channels']);ix=[channels.index(c) for c in POST]
        for wi,win in ((1,'early'),(2,'late')):
            x[s,win]=((z['features'][:,wi]-z['baseline'])[:,ix],z['labels'])
    out.append(('Turner',x,40,np.eye(40),
        {'full39':helmert(40,full=False),'Fourier1_2D':fourier(40,[1]),'Fourier12_4D':fourier(40,[1,2])}))
    return out

for dataset,data,n,L,spaces in get_data():
    state_n=len(L);P=np.roll(np.eye(state_n),1,axis=1)
    windows=sorted(set(w for s,w in data))
    for win in windows:
        split={s:halves(*data[s,win],n) for s in (1,2)}
        for source,target in ((1,2),(2,1)):
            for trainhalf in (0,1):
                fold=f'{source}to{target}_train{trainhalf}'
                xa,ya,ia=split[source][trainhalf];xc,yc,ic=split[source][1-trainhalf]
                xb,yb,ib=split[target][trainhalf];xd,yd,id_=split[target][1-trainhalf]
                assert len(set(ia)&set(ic))==0 and len(set(ib)&set(id_))==0
                W=whiten_fit(xa,ya,n)
                AA=L@means(xa@W,ya,n);CC=L@means(xc@W,yc,n)
                BB=L@means(xb@W,yb,n);DD=L@means(xd@W,yd,n)
                SPLIT.append({'dataset':dataset,'window':win,'fold':fold,'reference_train_n':len(ya),'reference_test_n':len(yc),
                    'target_train_n':len(yb),'target_test_n':len(yd),'within_partition_train_test_overlap':0})
                for space,C in spaces.items():
                    d={'dataset':dataset,'window':win,'fold':fold,'space':space}
                    A=C@AA;B=C@BB;At=C@CC;Bt=C@DD
                    err=C@P-(C@P@C.T)@C
                    assert np.max(np.abs(err))<1e-12
                    sv=np.linalg.svd(A,compute_uv=False)
                    noise=np.linalg.norm((A-At)/2,2)
                    nr=len(A)
                    # min singular includes structural zeros if rows > channels.
                    smin=float(sv[-1]) if len(sv)==nr else 0.
                    Ta={k:templates(P,AA,GRIDS[k],k,C) for k in KINDS}
                    Tc={k:templates(P,CC,GRIDS[k],k,C) for k in KINDS}
                    Ra={k:residualize_current(Ta[k],A) for k in KINDS}
                    for i,k in enumerate(KINDS):
                        for k2 in KINDS[i+1:]:
                            rg=ray_geometry(Ra[k],Ra[k2]);v=rg['residual_per_unit']
                            if np.isfinite(v).any():
                                idx=np.unravel_index(np.nanargmin(v),v.shape)
                                GEOM.append({**d,'comparison':k+'_vs_'+k2,'minimum_unit_cone_residual_grid':float(v[idx]),
                                    'gamma_1':float(GRIDS[k][idx[0]]),'gamma_2':float(GRIDS[k2][idx[1]]),'undefined':False})
                            else:GEOM.append({**d,'comparison':k+'_vs_'+k2,'minimum_unit_cone_residual_grid':np.nan,
                                             'gamma_1':np.nan,'gamma_2':np.nan,'undefined':True})
                    # Gamma-only information after refitting both gains.
                    for k in ('forward','backward'):
                        rg=ray_geometry(Ra[k],Ra[k]);vv=rg['residual_per_unit']
                        for g0 in (.2,.5,.8):
                            gi=int(round(g0*100));far=np.abs(G-g0)>=.1999999
                            PROFILE.append({**d,'kind':k,'gamma':g0,'min_residual_for_gamma_error_ge_point2':float(np.nanmin(vv[gi,far])),
                                'maximum_ray_difference_over_all_gamma':float(np.nanmax(rg['unit_chord'][gi]))})
                    D=kernel_design(P,AA,3,C,True)
                    ds=np.linalg.svd(D,compute_uv=False)
                    sr=ds/max(float(ds.max()),np.finfo(float).tiny)
                    kr=int(np.sum(sr>1e-10))
                    delta_a=C@(P-P.T)@AA;delta_c=C@(P-P.T)@CC
                    cos=float(np.sum(delta_a*delta_c)/(np.linalg.norm(delta_a)*np.linalg.norm(delta_c)))
                    DESIGN.append({**d,'n_state_contrasts':nr,'n_channels':A.shape[1],'numeric_reference_rank':int(np.linalg.matrix_rank(A)),
                        'full_reference_smin':smin,'empirical_split_operator_noise':noise,'smin_to_empirical_noise':smin/noise,
                        'kernel_lags_1to3_rank_after_current':kr,'kernel_smin_to_smax':float(sr[-1]),
                        'design_reduction_error':float(np.linalg.norm(err))})
                    REPEAT.append({**d,'forward_backward_contrast_repeat_cosine':cos,
                         'contrast_mean_norm':float(np.linalg.norm((delta_a+delta_c)/2)),
                         'contrast_difference_half_norm':float(np.linalg.norm((delta_a-delta_c)/2))})
                    # Independent target-repeat cross energy. The metric uses
                    # source A only, never either target repeat B or D.
                    db=(templates(P,BB,[.5],'forward',C)[0]-templates(P,BB,[.5],'backward',C)[0])
                    dd=(templates(P,DD,[.5],'forward',C)[0]-templates(P,DD,[.5],'backward',C)[0])
                    cross=float(np.sum(db*dd)); naive=float((np.sum(db*db)+np.sum(dd*dd))/2)
                    REPEAT[-1].update({'independent_target_delta_gamma':.5,'target_cross_energy':cross,
                        'target_naive_mean_squared_energy':naive,'target_cross_to_naive_energy':cross/naive,
                        'cross_energy_assumption':'zero_mean_independent_target_repeat_errors_conditional_on_source_metric',
                        'cross_energy_caveat':'same_person_halves_and_fold_reuse_no_population_interval'})
                    current_c=float(np.sum(A*B)/np.sum(A*A))
                    truth=[('current',0.,0.)]+[(k,g,rho) for k in ('forward','backward') for g in (.2,.5,.8) for rho in (.1,.3,1.)]
                    for true_kind,g0,rho in truth:
                        if true_kind=='current':
                            Ytr=B;Yte=Bt;alpha0=0.
                        else:
                            nb=templates(P,BB,[g0],true_kind,C)[0]
                            nd=templates(P,DD,[g0],true_kind,C)[0]
                            rnb=residualize_current(nb[None],B)[0]
                            alpha0=rho*np.linalg.norm(B)/np.linalg.norm(rnb)
                            Ytr=B+alpha0*nb;Yte=Bt+alpha0*nd
                        cc=float(np.sum(A*Ytr)/np.sum(A*A))
                        null_sse=float(np.sum((Yte-cc*At)**2))
                        norm=float(np.sum(Yte*Yte))
                        base={**d,'true_kind':true_kind,'true_gamma':g0,'prediction_residual_strength':rho,'true_alpha':alpha0}
                        CV.append({**base,'candidate':'current','gamma_hat':np.nan,'current_gain':cc,'predictive_gain':0.,
                           'train_sse':float(np.sum((Ytr-cc*A)**2)),'test_sse':null_sse,'test_target_norm2':norm,
                           'test_relative_sse':null_sse/norm,'test_improvement_vs_current':0.})
                        for k in KINDS:
                            fit=fit_affine(A,Ta[k],Ytr);j=int(np.argmin(fit['sse']))
                            pred=fit['current_gain'][j]*At+fit['predictive_gain'][j]*Tc[k][j]
                            loss=float(np.sum((Yte-pred)**2))
                            # Equivalent minimizers are reported instead of fabricating a unique gamma.
                            eps=1e-10*max(float(np.sum(Ytr*Ytr)),1.)
                            tie=GRIDS[k][fit['sse']<=fit['sse'][j]+eps]
                            gamma_point=float(GRIDS[k][j]) if len(tie)<len(G) and fit['predictive_gain'][j]>1e-12 else np.nan
                            CV.append({**base,'candidate':k,'gamma_hat':gamma_point,
                                'gamma_grid_argmin_audit_only':float(GRIDS[k][j]),'gamma_train_equal_min_low':float(tie.min()),
                                'gamma_train_equal_min_high':float(tie.max()),'current_gain':float(fit['current_gain'][j]),
                                'predictive_gain':float(fit['predictive_gain'][j]),'train_sse':float(fit['sse'][j]),
                                'test_sse':loss,'test_target_norm2':norm,'test_relative_sse':loss/norm,
                                'test_improvement_vs_current':(null_sse-loss)/null_sse})
        print('completed',dataset,win,flush=True)

for name,rows in [('geometry',GEOM),('gamma_geometry',PROFILE),('design_diagnostics',DESIGN),('repeat_contrast',REPEAT),('calibration_folds',CV),('split_audit',SPLIT)]:
    pd.DataFrame(rows).to_csv(O/f'{name}.csv',index=False)
df=pd.DataFrame(CV)
keys=['dataset','window','space','true_kind','true_gamma','prediction_residual_strength']
agg=df.groupby(keys+['candidate'],dropna=False)[['test_sse','test_target_norm2']].sum().reset_index()
agg['test_relative_sse']=agg.test_sse/agg.test_target_norm2
agg.to_csv(O/'calibration_aggregate.csv',index=False)
win=[]
for key,g in agg.groupby(keys,dropna=False):
    best=g.loc[g.test_sse.idxmin()]
    tol=1e-10*float(best.test_target_norm2)
    tie=g.loc[g.test_sse<=float(best.test_sse)+tol,'candidate'].tolist()
    classes=sorted(set('symmetric_equivalence' if x in ('symmetric_transition','balanced_fb') else x for x in tie))
    base=dict(zip(keys,key));base.update({'best_candidate':best.candidate,
        'best_set':'|'.join(sorted(tie)),'number_numerically_tied':len(tie),
        'best_class_set':'|'.join(classes),'number_tied_classes':len(classes),
        'unique_numerical_winner':len(classes)==1,'best_relative_sse':float(best.test_relative_sse)})
    current=float(g.loc[g.candidate=='current','test_sse'].iloc[0]);base['best_improvement_vs_current']=(current-float(best.test_sse))/current
    if base['true_kind']!='current':
        tt=g.loc[g.candidate==base['true_kind']].iloc[0];base['true_family_relative_sse']=float(tt.test_relative_sse)
        base['true_family_improvement_vs_current']=(current-float(tt.test_sse))/current
    base['true_family_in_best_set']=base['true_kind'] in classes
    base['unique_correct_family']=len(classes)==1 and base['true_kind'] in classes
    win.append(base)
pd.DataFrame(win).to_csv(O/'calibration_selection.csv',index=False)
(O/'input_manifest.json').write_text(json.dumps(INPUTS,indent=2))
print('COMPLETE',len(CV),'fold-candidate evaluations',len(win),'calibration conditions',flush=True)
