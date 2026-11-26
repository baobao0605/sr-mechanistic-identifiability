"""Ekman localizer analysis, frozen 2026-11-24; no simulated participants.

The primary analyses are of observed, cross-validated confusion probabilities.
They are NOT direct neural basis estimates, contamination eta estimates, or a
reconstruction of a clean pre-learning reference. HPa is a source filename only.
"""
from pathlib import Path
import json, hashlib, warnings
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import ttest_1samp, ttest_rel
from threadpoolctl import threadpool_limits

threadpool_limits(2)
O = Path(__file__).resolve().parent
A = O/'input/analysis_scripts'
F = A/'02_Analyses_and_article_Figures'
assert (O/'analysis_contract.json').exists()
subjects = np.loadtxt(A/'subjects_final.txt',dtype=str)
seqid = dict(np.loadtxt(A/'sequence_ids.txt',dtype=str))
mapping={1:([0,4,6,2],[1,3,5,7]),2:([5,1,3,7],[0,4,6,2]),3:([1,5,7,3],[0,4,6,2]),4:([6,2,4,0],[3,1,5,7]),5:([2,6,0,4],[7,3,1,5]),6:([7,3,5,1],[6,2,4,0]),7:([3,7,1,5],[0,4,6,2]),8:([4,0,2,6],[3,1,5,7])}
inputs={'HP_Figure6':F/'Figure_06/data/N35_HPa_localizer_LR_loo_probs.npy',
        'V1_Figure6':F/'Figure_06/data/N35_V1_localizer_LR_loo_probs.npy',
        'HP_Figure5_sensitivity':F/'Figure_05/data/N35_HP_localizer_LR_loo_probs.npy'}
data={k:np.load(p) for k,p in inputs.items()}
rng=np.random.default_rng(2026112402)

def stats(v):
    v=np.asarray(v,float); n=len(v)
    bs=v[rng.integers(0,n,(20000,n))].mean(1)
    ci=np.quantile(bs,[.025,.975]); t,p=ttest_1samp(v,0)
    null=(rng.choice([-1,1],(100000,n))*v).mean(1)
    return dict(n=n,mean=float(v.mean()),ci_low=float(ci[0]),ci_high=float(ci[1]),
                t=float(t),p_t=float(p),p_signflip=float((1+np.sum(abs(null)>=abs(v.mean())))/(len(null)+1)))

def holm(p):
    p=np.asarray(p); ii=np.argsort(p); q=np.empty(len(p))
    q[ii]=np.minimum(1,np.maximum.accumulate(p[ii]*(len(p)-np.arange(len(p)))))
    return q

def design(idx,include_row=True):
    i,j=np.where(~np.eye(8,dtype=bool))
    ang=np.minimum(abs(i-j),8-abs(i-j))
    pos=np.full(8,-1);pos[idx]=np.arange(4)
    learned=(pos[i]>=0)&(pos[j]>=0)
    lag=pos[j]-pos[i]
    # Full ring-distance categories absorb any isotropic spatial falloff.
    z=[np.ones(len(i))]
    if include_row:z += [(i==k).astype(float) for k in range(1,8)]
    z += [(j==k).astype(float) for k in range(1,8)]
    z += [(ang==k).astype(float) for k in range(1,4)]
    z += [learned.astype(float)]
    Z=np.column_stack(z)
    T=-abs(lag)*learned
    U=np.sign(lag)*learned
    return i,j,Z,T.astype(float),U.astype(float),lag,learned,ang

rows=[];aud=[];pred=[];curves=[]
for name,arr in data.items():
    for s,sub in enumerate(subjects):
        idx,ctrl=mapping[int(seqid[sub])]
        i,j,Z,T,U,lag,learned,ang=design(idx)
        y=arr[s,i,j]
        X=np.column_stack([Z,T,U]);coef=np.linalg.lstsq(X,y,rcond=None)[0]
        R=np.eye(len(y))-Z@np.linalg.pinv(Z)
        unique=np.column_stack([R@T,R@U]); ss=np.linalg.svd(unique,compute_uv=False)
        assert np.linalg.matrix_rank(X)==np.linalg.matrix_rank(Z)+2
        aud.append(dict(source=name,subject=sub,n_cells=len(y),n_nuisance_columns=Z.shape[1],
            nuisance_rank=np.linalg.matrix_rank(Z),full_rank=np.linalg.matrix_rank(X),
            unique_temporal_norm=np.linalg.norm(R@T),unique_direction_norm=np.linalg.norm(R@U),
            unique_min_sv=ss[-1],row_probability_sum_max_error=np.max(abs(arr[s].sum(1)-1))))
        for metric,value in [('temporal_proximity_unique',coef[-2]),('direction_unique',coef[-1])]:
            rows.append(dict(source=name,subject=sub,metric=metric,value=value))
        mat=arr[s][np.ix_(idx,idx)]
        for k in [-3,-2,-1,0,1,2,3]:
            curves.append(dict(source=name,subject=sub,lag=k,value=np.diag(mat,k=k).mean()))
        # Predict all held cue cells using only other cue rows. No estimand treats
        # those correlated folds as independent participants.
        i,j,Z,T,U,lag,learned,ang=design(idx,include_row=False)
        y=arr[s,i,j]
        for held in range(8):
            tr=i!=held;te=~tr;zt=Z[tr];zz=Z[te]
            q=np.eye(sum(tr))-zt@np.linalg.pinv(zt);ry=q@y[tr]
            for model in ['spatial_baseline','forward','backward','symmetric','temporal_linear']:
                best=None
                for g in ([.5] if model in ['spatial_baseline','temporal_linear'] else np.linspace(.025,.975,39)):
                    x=np.zeros(len(y))
                    if model=='temporal_linear':x=T
                    elif model=='forward':x=np.where(learned&(lag>0),g**np.maximum(abs(lag)-1,0),0)
                    elif model=='backward':x=np.where(learned&(lag<0),g**np.maximum(abs(lag)-1,0),0)
                    elif model=='symmetric':x=np.where(learned,g**np.maximum(abs(lag)-1,0),0)
                    xr=q@x[tr];den=xr@xr
                    amp=max(0,float(xr@ry/den)) if den>1e-20 else 0
                    b=np.linalg.lstsq(zt,y[tr]-amp*x[tr],rcond=None)[0]
                    loss=np.sum((y[tr]-zt@b-amp*x[tr])**2)
                    if best is None or loss<best[0]:best=(loss,g,amp,b,x)
                loss,g,amp,b,x=best;yh=zz@b+amp*x[te]
                pred.append(dict(source=name,subject=sub,held_cue=held,model=model,training_SSE=loss,
                    test_SSE=np.sum((y[te]-yh)**2),gamma=g,associative_gain=amp))
    print('DONE MATRIX',name,flush=True)

df=pd.DataFrame(rows);df.to_csv(O/'localizer_unique_subject.csv',index=False)
st=[]
for (name,metric),d in df.groupby(['source','metric']):st.append(dict(source=name,metric=metric,**stats(d.value)))
st=pd.DataFrame(st);st['primary']=~st.source.str.contains('sensitivity');st['holm_primary4']=np.nan
sel=st.primary;st.loc[sel,'holm_primary4']=holm(st.loc[sel,'p_signflip'])
st.to_csv(O/'localizer_unique_group.csv',index=False)
pd.DataFrame(aud).to_csv(O/'localizer_design_audit.csv',index=False)
pd.DataFrame(curves).to_csv(O/'localizer_lag_curves.csv',index=False)
pr=pd.DataFrame(pred);pr.to_csv(O/'localizer_heldcue_predictions.csv',index=False)
su=pr.groupby(['source','subject','model']).test_SSE.sum().unstack('model')
out=[]
for (source,sub),s in su.iterrows():
    for m in ['forward','backward','symmetric','temporal_linear']:
        out.append(dict(source=source,subject=sub,model=m,improvement=1-s[m]/s.spatial_baseline))
out=pd.DataFrame(out);out.to_csv(O/'localizer_heldcue_subject.csv',index=False)
pst=[dict(source=r,model=m,**stats(d.improvement)) for (r,m),d in out.groupby(['source','model'])]
pst=pd.DataFrame(pst);pst['primary_source']=~pst.source.str.contains('sensitivity');pst['holm_eight_predictive_comparisons']=np.nan
sel=pst.primary_source;pst.loc[sel,'holm_eight_predictive_comparisons']=holm(pst.loc[sel,'p_signflip'])
pst.to_csv(O/'localizer_heldcue_group.csv',index=False)

# Source-code execution audit: rerun only numerical cells, save original cell
# hashes/execution counters/output. This makes a stale notebook discrepancy visible.
book=json.loads((F/'Figure_06/Figure_06_tuning.ipynb').read_text())
notes=[]
for k in [4,5,8,9,12,13,14]:
    c=book['cells'][k];text=''.join(c.get('source',[]))
    notes.append(dict(cell=k,execution_count=c.get('execution_count'),source=text,
        source_sha256=hashlib.sha256(text.encode()).hexdigest(),stored_outputs=c.get('outputs',[])))
(O/'figure6_notebook_execution_audit.json').write_text(json.dumps(notes,ensure_ascii=False,indent=2))

def original_profiles(arr):
    out=[]
    for s,sub in enumerate(subjects):
        idx,_=mapping[int(seqid[sub])];mat=arr[s][np.ix_(idx,idx)]
        out.append([np.diag(mat,k=k).mean() for k in [-3,-2,-1,0,1,2,3]])
    return np.array(out)

def bic_sse(y,yh,k):
    sse=np.sum((y-yh)**2);n=len(y)
    return n*(np.log(2*np.pi*sse/n)+1)+k*np.log(n)

def f_sr(x,a,g,c):
    out=np.zeros(7);out[3:]=a*g**np.arange(4)+c;return out

modrows=[];pairs=[];decomp=[]
orig_sp=np.array([9.89,9.89,14,0,14,9.89,9.89])
correct_sp=np.array([7*np.sqrt(2),7*np.sqrt(2),(28+7*np.sqrt(2))/3,0,(28+7*np.sqrt(2))/3,7*np.sqrt(2),7*np.sqrt(2)])
temp=3-abs(np.arange(-3,4))
for name in data:
    yp=original_profiles(data[name])
    for s,sub in enumerate(subjects):
        y=yp[s]
        for m,x in [('source_temporal',temp),('source_spatial',orig_sp),('correct_spatial',correct_sp),('fitted_null',np.arange(7)==3)]:
            X=np.column_stack([x,np.ones(7)]);yh=X@np.linalg.lstsq(X,y,rcond=None)[0]
            modrows.append(dict(source=name,subject=sub,model=m,BIC=bic_sse(y,yh,2),SSE=np.sum((y-yh)**2),gamma=np.nan))
        yh=np.full(7,(y[0]+y[-1])/2);yh[3]=y[3]
        modrows.append(dict(source=name,subject=sub,model='source_overwritten_null',BIC=bic_sse(y,yh,1),SSE=np.sum((y-yh)**2),gamma=np.nan))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            cf,_=curve_fit(f_sr,np.arange(7),y,maxfev=10000,bounds=([-np.inf]*3,[np.inf]*3))
        yh=f_sr(np.arange(7),*cf)
        modrows.append(dict(source=name,subject=sub,model='source_SR',BIC=bic_sse(y,yh,2),SSE=np.sum((y-yh)**2),gamma=cf[1]))
        # Fair in-sample comparison: separate current response for every family,
        # shared baseline, finite-chain positive prediction, bounded discount.
        best={}
        for model in ['current','forward','backward','symmetric','spatial']:
            for g in ([.5] if model in ['current','spatial'] else np.linspace(0,1,401)):
                lag=np.arange(-3,4)
                if model=='current':x=np.zeros(7)
                elif model=='spatial':x=-correct_sp
                elif model=='forward':x=np.where(lag>0,g**np.maximum(abs(lag)-1,0),0)
                elif model=='backward':x=np.where(lag<0,g**np.maximum(abs(lag)-1,0),0)
                else:x=np.where(lag!=0,g**np.maximum(abs(lag)-1,0),0)
                Z=np.column_stack([np.ones(7),lag==0]);q=np.eye(7)-Z@np.linalg.pinv(Z);xr=q@x;yr=q@y
                amp=max(0,xr@yr/(xr@xr)) if xr@xr>1e-20 else 0
                b=np.linalg.lstsq(Z,y-amp*x,rcond=None)[0];pred=Z@b+amp*x;sse=np.sum((y-pred)**2)
                if model not in best or sse<best[model][0]:best[model]=(sse,g,amp,pred)
            sse,g,amp,yh=best[model]
            modrows.append(dict(source=name,subject=sub,model='fair_'+model,BIC=bic_sse(y,yh,2 if model=='current' else 3 if model=='spatial' else 4),SSE=sse,gamma=g))
        # Diagnostic: the temporal model can win partly because its central peak
        # is tied to the offdiagonal slope. Remove the diagonal and directly test
        # the observed offdiagonal proximity slope, before spatial adjustment.
        mask=np.arange(7)!=3;X=np.column_stack([np.ones(6),-abs(np.arange(-3,4))[mask]])
        b=np.linalg.lstsq(X,y[mask],rcond=None)[0]
        decomp.append(dict(source=name,subject=sub,offdiag_proximity_slope=b[1],current_minus_offdiag=y[3]-y[mask].mean()))
    print('DONE SOURCE AUDIT',name,flush=True)
md=pd.DataFrame(modrows);md.to_csv(O/'figure6_model_audit_subject.csv',index=False)
for source,dd in md.groupby('source'):
    w=dd.pivot(index='subject',columns='model',values='BIC')
    for a,b in [('source_temporal','source_spatial'),('source_temporal','source_overwritten_null'),('source_temporal','source_SR'),('source_overwritten_null','source_temporal'),('source_overwritten_null','source_spatial'),('source_overwritten_null','source_SR'),('source_temporal','fitted_null'),('source_temporal','correct_spatial'),('fair_symmetric','fair_current'),('fair_forward','fair_current'),('fair_forward','fair_backward')]:
        pairs.append(dict(source=source,a=a,b=b,**stats(w[a]-w[b])))
pd.DataFrame(pairs).to_csv(O/'figure6_model_audit_group.csv',index=False)
dc=pd.DataFrame(decomp);dc.to_csv(O/'localizer_diagonal_sensitivity_subject.csv',index=False)
pd.DataFrame([dict(source=r,metric=m,**stats(d[m])) for r,d in dc.groupby('source') for m in ['offdiag_proximity_slope','current_minus_offdiag']]).to_csv(O/'localizer_diagonal_sensitivity_group.csv',index=False)

manifest={k:{'path':str(p.relative_to(O)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'shape':list(data[k].shape)} for k,p in inputs.items()}
manifest['HP_filename_discrepancy']={'max_absolute_difference':float(np.max(abs(data['HP_Figure6']-data['HP_Figure5_sensitivity']))),'meaning_of_HPa':'unresolved; not assumed anterior or whole hippocampus'}
manifest['scope']='Observed participant matrices; no individual run/repetition axis. Inference uses 35 participant units, folds are not independent samples.'
(O/'localizer_source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
print(st.to_string(index=False))
print(pst.to_string(index=False))
