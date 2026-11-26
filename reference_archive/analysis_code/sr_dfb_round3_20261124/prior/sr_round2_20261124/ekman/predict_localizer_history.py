"""Held chronological block-group prediction from actual localizer history.

Actual 64 localizer blocks per participant, not simulations. Test units are whole
rounds of eight locations, with neighbouring rounds excluded from fitting. The
eight folds and two participants are descriptive; no voxel/trial pseudoreplication
is used for population inference. Hyperparameters use only training rounds.
"""
from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.stats import gamma
from threadpoolctl import threadpool_limits
threadpool_limits(2)
O=Path(__file__).resolve().parent;OLD=O.parents[1]/'restored/complete/branches/D02'
events=pd.read_csv(O/'actual_localizer_event_history.csv')
rows=[];checks=[];alphas=[0.,.1,1.,10.,100.]

def fit_predict(X,Y,tr,te,alpha):
    # Train-only centering. Current/history stimulus regressors have common 0/1
    # units; nuisance time is standardized by its known experiment range.
    xm=X[tr].mean(0);ym=Y[tr].mean(0)
    xt=X[tr]-xm;yt=Y[tr]-ym
    beta=np.linalg.pinv(xt.T@xt+alpha*np.eye(X.shape[1]))@xt.T@yt
    return (X[te]-xm)@beta+ym

for sub,ev in events.groupby('subject',sort=True):
    ev=ev.sort_values(['label','rep']).reset_index(drop=True)
    label=ev.label.to_numpy(int);onset=ev.onset.to_numpy();previous=ev.previous_label.to_numpy(int)
    null=np.loadtxt(O/f'input/analysis_scripts/subject_data/{sub}/evs/localizer_null.txt')[:,0]
    rounds=np.searchsorted(null,onset,side='right')-1
    assert all(np.array_equal(np.sort(label[rounds==g]),np.arange(8)) for g in range(8))
    nextlab=np.full(64,8,int)
    for n in range(64):
        cand=np.where(onset>onset[n])[0]
        if len(cand):
            nx=cand[np.argmin(onset[cand])]
            if rounds[nx]==rounds[n]:nextlab[n]=label[nx]
    hot=np.eye(8)[label]
    prev=np.where((previous<8)[:,None],np.eye(9)[np.maximum(previous,0)][:,:8],0)
    nxt=np.eye(9)[nextlab][:,:8]
    nuisance=np.column_stack([(previous==8).astype(float),(onset-onset.mean())/(onset.max()-onset.min())])
    db=np.load(OLD/f'{sub}_localizer_roi_timeseries.npz');n=db['timeseries'].shape[1]
    times=(np.arange(n)+.5)*1.5
    ind=np.arange(n);dx=ind[None,:]-ind[:,None];w=np.exp(-.5*(dx/(100/3))**2)
    s0=w.sum(1);s1=(w*dx).sum(1);s2=(w*dx**2).sum(1)
    high=np.eye(n)-(s2[:,None]*w-s1[:,None]*w*dx)/(s0*s2-s1*s1)[:,None]
    physical={}
    for mean in [4.,6.,8.]:
        L=np.column_stack([gamma.cdf(times-t,a=4,scale=mean/4)-gamma.cdf(times-t-13.5,a=4,scale=mean/4) for t in onset])
        G=high@np.column_stack([L[:,label==k].sum(1) for k in range(8)])
        physical[mean]=G
    for window,start,end in [('primary3to13p5s',2,9),('late9to13p5s',6,9),('code6to13p5s_sensitivity',4,9)]:
        avg=np.zeros((64,n))
        for k,t in enumerate(onset):avg[k,int(np.floor(t/1.5))+start:int(np.floor(t/1.5))+end]=1/(end-start)
        designs={'current':[(None,np.column_stack([hot,nuisance]))],
                 'current_previous':[(None,np.column_stack([hot,prev,nuisance]))],
                 'current_next_control':[(None,np.column_stack([hot,nxt,nuisance]))],
                 'physical_LTI':[(m,np.column_stack([avg@physical[m],nuisance])) for m in [4.,6.,8.]]}
        for roi in ['V1','HP']:
            Y=db['timeseries'][db[roi]].T.astype(float);mu=Y.mean(0)
            Y=100*(Y-mu)/np.maximum(abs(mu),1);Y=avg@Y
            for held in range(8):
                te=rounds==held;trainrounds=[g for g in range(8) if abs(g-held)>1];tr=np.isin(rounds,trainrounds)
                for model,variants in designs.items():
                    best=None
                    for param,X in variants:
                        for alpha in alphas:
                            loss=0
                            for inner in trainrounds:
                                va=rounds==inner;train=np.isin(rounds,[g for g in trainrounds if g!=inner])
                                pred=fit_predict(X,Y,train,va,alpha)
                                loss+=np.sum((Y[va]-pred)**2)
                            if best is None or loss<best[0]:best=(loss,param,X,alpha)
                    loss,param,X,alpha=best;pred=fit_predict(X,Y,tr,te,alpha)
                    sse=np.sum((Y[te]-pred)**2);base=np.sum((Y[te]-Y[tr].mean(0))**2)
                    rows.append(dict(subject=sub,roi=roi,window=window,test_round=held,model=model,train_rounds='-'.join(map(str,trainrounds)),
                        train_trials=int(tr.sum()),test_trials=int(te.sum()),inner_SSE=loss,test_SSE=sse,meanbaseline_SSE=base,alpha=alpha,HRF_mean=param))
                    checks.append(dict(subject=sub,window=window,test_round=held,model=model,rank_train=np.linalg.matrix_rank(X[tr]-X[tr].mean(0)),n_columns=X.shape[1]))
    print('DONE HISTORY PREDICT',sub,flush=True)

d=pd.DataFrame(rows);d.to_csv(O/'localizer_history_prediction_folds.csv',index=False)
pd.DataFrame(checks).to_csv(O/'localizer_history_prediction_designs.csv',index=False)
summary=[]
for (sub,roi,win),g in d.groupby(['subject','roi','window']):
    w=g.pivot(index='test_round',columns='model',values='test_SSE')
    for model in ['current_previous','current_next_control','physical_LTI']:
        a=1-w[model]/w.current
        summary.append(dict(subject=sub,roi=roi,window=win,model=model,
            pooled_improvement=1-w[model].sum()/w.current.sum(),mean_fold_improvement=a.mean(),positive_folds=int((a>0).sum()),n_folds=len(a),
            inference='Two selected participants; overlapping folds are descriptive, not independent samples'))
s=pd.DataFrame(summary);s.to_csv(O/'localizer_history_prediction_summary.csv',index=False)
(O/'localizer_history_prediction_contract.json').write_text(json.dumps({'stage':'Targeted follow-up after actual event imbalance was established; no outcome-based model or fold selection',
 'test':'8 chronological rounds, each with all 8 locations. Exclude adjacent chronological rounds from outer training. Inner leave-one-round tuning among remaining training rounds.',
 'models':['current state + after-null + elapsed time','same + previous actually presented location','same + next location negative/control (author highpass is noncausal)','same-size physical current+history design with HRF mean4/6/8 chosen using only training'],
 'ridge_grid':alphas,'units':'percent-signal-change raw anatomical ROI windows; all voxels retained, no main-task outcome selection',
 'scope':'Prediction of measured localizer windows and sensitivity of reference construction. This is not a test that cognitive SR is present or absent.'},ensure_ascii=False,indent=2))
print(s.to_string(index=False))
