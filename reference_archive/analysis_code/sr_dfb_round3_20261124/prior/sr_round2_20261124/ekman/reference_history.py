"""Actual-event HRF mixing and V1 readout stability, Ekman S01/S02.

Analytical convolutions on actual EVs, not synthetic participants. HRF mixing is
conditional on an explicit linear time-invariant response. It cannot estimate
learned contamination in an unknown neural basis, or directly correct nonlinear
classifier probabilities. Repeated reference partitions are sensitivity settings,
not independent observations or confidence intervals.
"""
from pathlib import Path
from itertools import combinations
import json,io,hashlib
import numpy as np
import pandas as pd
from scipy.stats import gamma
from scipy.linalg import helmert
from threadpoolctl import threadpool_limits
threadpool_limits(2)

O=Path(__file__).resolve().parent
OLD=O.parents[1]/'restored/complete/branches/D02'
# O is .../sr_round2_20261124/ekman; parents[1] is workspace.
S=O/'input';A=S/'analysis_scripts'
raw=np.load(OLD/'measured_patterns.npz');glm=np.load(OLD/'GLM_measured_patterns.npz')
mapping={'S01':([3,7,1,5],[0,4,6,2]),'S02':([2,6,0,4],[7,3,1,5])}
H8=helmert(8,full=False)
J=np.eye(4)-np.ones((4,4))/4
W=np.sign(np.arange(4)[None,:]-np.arange(4)[:,None]);JWJ=J@W@J
splits=list(combinations(range(8),4))

def direction(mat):
    return {'triangular_FminusB':float(mat[np.triu_indices(4,1)].mean()-mat[np.tril_indices(4,-1)].mean()),
            'rowcol_orth_direction':float(np.sum(mat*JWJ)/np.sum(JWJ**2))}

def highpass(n):
    ind=np.arange(n);dx=ind[None,:]-ind[:,None]
    w=np.exp(-.5*(dx/(100/3))**2);s0=w.sum(1);s1=(w*dx).sum(1);s2=(w*dx**2).sum(1)
    h=(s2[:,None]*w-s1[:,None]*w*dx)/(s0*s2-s1*s1)[:,None]
    return np.eye(n)-h

def select_voxels(X,labels,top=25):
    means=np.array([X[labels==k].mean(0) for k in range(8)])
    res=X-means[labels];var=np.sum(res**2,axis=0)/(len(X)-8)
    # Each state has equal repetitions; common denominator scaling does not
    # change per-location ranking, so ratio is sufficient for selection.
    t=(8/7)*(means-means.mean(0))/np.sqrt(np.maximum(var,1e-20))
    picks=[np.argsort(t[k])[::-1][:top] for k in range(8)]
    counts=np.bincount(np.concatenate(picks),minlength=X.shape[1])
    picks=[x[counts[x]==1] for x in picks];m=min(map(len,picks))
    if m==0:raise RuntimeError('No unique spatial groups')
    return [x[:m] for x in picks]

mixrows=[];mixdetails=[];historyrows=[];designrows=[];sensitivity=[];selectionrows=[]
cached={}
for sub,(idx,ctrl) in mapping.items():
    evdir=A/f'subject_data/{sub}/evs'
    db=np.load(OLD/f'{sub}_localizer_roi_timeseries.npz');n=db['timeseries'].shape[1]
    des=np.loadtxt(io.StringIO((S/f'data/{sub}/func/localizer_hp100_fwhm5_v3.feat/design.mat').read_text().split('/Matrix')[1]))
    times=(np.arange(n)+.5)*1.5;Fhp=highpass(n)
    events=[]
    for k in range(8):
        for rep,(on,dur,amp) in enumerate(np.loadtxt(evdir/f'localizer_dot{k+1}.txt')):
            events.append(dict(onset=on,duration=dur,amplitude=amp,label=k,rep=rep))
    # This order matches the archived measured_patterns.npz exactly.
    onsets=np.array([e['onset'] for e in events]);labels=np.array([e['label'] for e in events]);rep=np.array([e['rep'] for e in events])
    nulls=[dict(onset=t,duration=d,amplitude=a,label=8,rep=j) for j,(t,d,a) in enumerate(np.loadtxt(evdir/'localizer_null.txt'))]
    ordered=sorted(events+nulls,key=lambda x:x['onset'])
    predecessor=np.zeros((8,9),int)
    for i,e in enumerate(ordered):
        if e['label']<8:
            prv=ordered[i-1]['label'] if i else -1
            if prv>=0:predecessor[e['label'],prv]+=1
            historyrows.append(dict(subject=sub,**e,previous_label=prv,previous_onset=ordered[i-1]['onset'] if i else np.nan))
    pd.DataFrame(predecessor,columns=[f'previous_{k}' for k in range(9)]).to_csv(O/f'{sub}_actual_predecessor_counts.csv',index_label='current_label')

    # Canonical source HRF; sensitivity means 4,6,8s with fixed coefficient of
    # variation .5. Source file uses mean6, SD3, a=4, scale=1.5.
    for mean in [4.,6.,8.]:
        L=np.column_stack([e['amplitude']*(gamma.cdf(times-e['onset'],a=4,scale=mean/4)-gamma.cdf(times-e['onset']-e['duration'],a=4,scale=mean/4)) for e in events])
        G=np.column_stack([L[:,labels==k].sum(1) for k in range(8)])
        for filtering in ['unfiltered_LTI','source_highpass100s']:
            GG=G if filtering=='unfiltered_LTI' else Fhp@G
            for window,start,end in [('3to6s',2,4),('3to13p5s',2,9),('9to13p5s',6,9)]:
                B=np.zeros((64,n))
                for e,on in enumerate(onsets):B[e,int(np.floor(on/1.5))+start:int(np.floor(on/1.5))+end]=1/(end-start)
                # Each row is an actual event's averaged exposure to all 8 labels.
                M=B@GG
                C=np.array([M[labels==k].mean(0) for k in range(8)])
                K=H8@C@H8.T;sv=np.linalg.svd(K,compute_uv=False)
                direct=np.mean(np.diag(C));N=C/direct
                row=dict(subject=sub,HRF_mean_s=mean,HRF_SD_s=mean/2,filtering=filtering,window=window,
                    diagonal_mean=direct,offdiag_mean=float(C[~np.eye(8,dtype=bool)].mean()),
                    full_rank=np.linalg.matrix_rank(C),full_condition=np.linalg.cond(C),
                    common_component_leakage_norm=np.linalg.norm((np.eye(8)-np.ones((8,8))/8)@C@np.ones(8)),
                    contrast_rank=np.linalg.matrix_rank(K),contrast_condition=sv[0]/sv[-1],contrast_smin=sv[-1],
                    max_abs_normalized_offdiag=float(np.max(abs(N[~np.eye(8,dtype=bool)]))),
                    **direction(N[np.ix_(idx,idx)]))
                mixrows.append(row)
                for i in range(8):
                    for j in range(8):mixdetails.append(dict(subject=sub,HRF_mean_s=mean,filtering=filtering,window=window,current=i,represented=j,mixing=C[i,j],relative_to_current=N[i,j]))
                cached[f'{sub}_{mean}_{filtering}_{window}']=C
        if mean==6.:
            Lhp=Fhp@L
            c=np.array([np.corrcoef((Fhp@G)[:,k],des[:,k])[0,1] for k in range(8)])
            # Original 8 condition columns replaced by 64 event columns; retain
            # original instruction/motion nuisance and an intercept.
            Z=np.column_stack([Lhp,des[:,8:],np.ones(n)])
            pinv=np.linalg.pinv(Z)
            Y=db['timeseries'][db['V1']].T.astype(float);mu=Y.mean(0);Y=100*(Y-mu)/np.maximum(abs(mu),1)
            trialbeta=(pinv@Y)[:64]
            cached[f'{sub}_V1_trialbeta_LSA']=trialbeta
            designrows.append(dict(subject=sub,n_volumes=n,n_columns=Z.shape[1],rank=np.linalg.matrix_rank(Z),condition=np.linalg.cond(Z),
                min_correlation_source8design=float(c.min()),max_correlation_source8design=float(c.max()),
                trial_betas='64 LSA estimates; correlated estimation errors, not independent raw observations'))
    np.savez_compressed(O/f'{sub}_reference_LSA.npz',trialbeta=trialbeta,label=labels,rep=rep,onset=onsets)
    for mode,X in [('window3to13p5s',raw[f'{sub}_V1_none_localizer']),('LSA_history_adjusted',trialbeta)]:
        for fold,trainreps in enumerate(splits):
            tr=np.isin(rep,trainreps);te=~tr
            groups=select_voxels(X[tr],labels[tr]);nv=len(groups[0])
            # Readout is built ONLY from localizer training repetition subsets.
            loc=np.array([[X[te&(labels==k)][:,g].mean() for g in groups] for k in range(8)])
            lm=loc[np.ix_(idx,idx)]-loc[np.ix_(idx,ctrl)].mean(1,keepdims=True)
            locdir=direction(lm)
            for run,M in enumerate(glm[f'{sub}_V1_main_betas'],1):
                out=np.array([[M[k,g].mean() for g in groups] for k in range(4)])
                mat=out[:,idx]-out[:,ctrl].mean(1,keepdims=True)
                for metric,value in direction(mat).items():
                    sensitivity.append(dict(subject=sub,mode=mode,partition=fold,train_repetitions='-'.join(map(str,trainreps)),
                        main_run=run,metric=metric,value=value,localizer_heldout_value=locdir[metric],voxels_per_location=nv))
            full=select_voxels(X,labels)
            overlap=np.mean([len(np.intersect1d(groups[k],full[k]))/len(groups[k]) for k in range(8)])
            selectionrows.append(dict(subject=sub,mode=mode,partition=fold,voxels_per_location=nv,mean_fraction_matching_full_reference=overlap))
    print('DONE REFERENCE',sub,flush=True)

pd.DataFrame(mixrows).to_csv(O/'actual_history_mixing_summary.csv',index=False)
pd.DataFrame(mixdetails).to_csv(O/'actual_history_mixing_cells.csv',index=False)
pd.DataFrame(historyrows).to_csv(O/'actual_localizer_event_history.csv',index=False)
pd.DataFrame(designrows).to_csv(O/'reference_LSA_design_audit.csv',index=False)
pd.DataFrame(selectionrows).to_csv(O/'reference_selection_stability.csv',index=False)
s=pd.DataFrame(sensitivity);s.to_csv(O/'reference_main_direction_sensitivity.csv',index=False)
g=s.groupby(['subject','mode','partition','metric']).value.mean().reset_index()
summ=[]
for (sub,mode,metric),a in g.groupby(['subject','mode','metric']):
    v=a.value.to_numpy()
    summ.append(dict(subject=sub,mode=mode,metric=metric,n_readout_partitions=len(v),mean=v.mean(),median=np.median(v),
        minimum=v.min(),maximum=v.max(),q025=np.quantile(v,.025),q975=np.quantile(v,.975),positive_fraction=np.mean(v>0),
        interpretation='Sensitivity across overlapping localizer partitions; quantiles are NOT confidence intervals'))
pd.DataFrame(summ).to_csv(O/'reference_main_direction_summary.csv',index=False)
np.savez_compressed(O/'actual_history_mixing_matrices.npz',**cached)
facts={sub:{'ordered_pair_min':int(np.min(np.loadtxt(O/f'{sub}_actual_predecessor_counts.csv',delimiter=',',skiprows=1)[:,1:9][~np.eye(8,dtype=bool)])),
            'ordered_pair_max':int(np.max(np.loadtxt(O/f'{sub}_actual_predecessor_counts.csv',delimiter=',',skiprows=1)[:,1:9])),
            'all_current_labels_follow_blank_once':True,'full_first_order_counterbalanced':False} for sub in mapping}
facts['limitations']=['N=2 selected participants only; no population inference','LTI mixing is exact for the specified HRF and supplied EVs; empirical learned-reference contamination remains unknown','V1 source-style selection on trial features is a sensitivity pipeline, not exact original FEAT selection','Main beta cache uses supplemental SPM HRF GLM lacking original motion regressors; sign changes do not prove published direction errors','Localizer LSA removes modeled inter-block HRF overlap; it does not remove genuine learned associations or fix an unknown neural basis']
(O/'reference_history_audit.json').write_text(json.dumps(facts,ensure_ascii=False,indent=2))
print(pd.DataFrame(mixrows).query("HRF_mean_s==6 and window=='3to13p5s'").to_string(index=False))
print(pd.DataFrame(summ).to_string(index=False))
