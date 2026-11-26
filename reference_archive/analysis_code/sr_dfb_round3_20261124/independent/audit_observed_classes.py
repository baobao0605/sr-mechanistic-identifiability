"""Align old and target-refit endpoints without overwriting original labels.

Old saved candidate labels are preserved. An additional observation-class
endpoint maps zero predictive gain and current-equivalent shapes to current,
using the same rule as the target-refit ablation.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from build_inputs import kernel
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
KEYS=['dataset','window','fold','space','true_kind','true_gamma','prediction_residual_strength']

def main():
    old=pd.read_csv(ROOT/'prior/sr_round2_20261124/math/calibration_folds.csv')
    f=pd.read_csv(HERE/'margin_folds.csv');profile=pd.read_csv(HERE/'profiled_decision_folds.csv')
    cache=np.load(HERE/'repeated_inputs.npz');meta=pd.read_csv(HERE/'input_folds.csv');eq={};rows=[]
    for m in meta.itertuples(index=False):
        A,S=cache[m.fold_id+'_A'],cache[m.fold_id+'_S'];ss=float(np.sum(S*S))
        for kind in ('forward','backward','symmetric_transition'):
            X=kernel(A,0.,kind)@S;R=X-S*(float(np.sum(S*X))/ss)
            eq[m.fold_id,kind]=bool(np.linalg.norm(R)<=1e-12*max(np.linalg.norm(X),1e-300))
    for r in f.itertuples(index=False):
        s=old.copy()
        for k in KEYS:s=s[s[k].eq(getattr(r,k))]
        assert len(s)==5
        row={k:getattr(r,k) for k in KEYS};row['fold_id']=r.fold_id
        tf=s[s.candidate.eq(r.true_kind)].iloc[0]
        row['true_family_fitted_predictive_gain']=float(tf.predictive_gain)
        row['true_family_shape_current_equivalent']=False if r.true_kind=='current' else eq[r.fold_id,r.true_kind]
        row['true_family_visible']=r.true_kind!='current' and tf.predictive_gain>1e-12 and not row['true_family_shape_current_equivalent']
        for mode,candidates in [('binary',['forward','backward']),('multi',['forward','backward','symmetric_transition','current'])]:
            part=s[s.candidate.isin(candidates)];best=part.test_sse.min()
            tol=1e-10*float(part.test_target_norm2.iloc[0]);ties=part[part.test_sse<=best+tol]
            raw=[];obs=[]
            for q in ties.itertuples(index=False):
                label='symmetric' if q.candidate=='symmetric_transition' else q.candidate;raw.append(label)
                current=q.candidate=='current' or q.predictive_gain<=1e-12 or eq.get((r.fold_id,q.candidate),False)
                obs.append('current' if current else label)
            raw='|'.join(sorted(set(raw)));obs='|'.join(sorted(set(obs)))
            assert raw==getattr(r,mode+'_winners')
            row['old_'+mode+'_label_winners']=raw;row['old_'+mode+'_observed_winners']=obs
            row[mode+'_label_unique_correct']=bool(getattr(r,mode+'_correct'))
            row[mode+'_visible_unique_correct']=bool(r.true_kind!='current' and obs==r.true_kind)
            if mode=='multi' and r.true_kind=='current':row[mode+'_visible_unique_correct']=obs=='current'
        rows.append(row)
    a=pd.DataFrame(rows);a.to_csv(HERE/'observed_class_audit.csv',index=False)
    d=a[a.true_kind.ne('current')];p=profile[profile.true_kind.ne('current')]
    summary={'n_direction_conditions':len(d),'boundary_rule':'Predictive gain <= 1e-12 or current-equivalent shape maps to current; same rule as target-refit profile. Original label columns preserved.','whole_sample':{},'paired_gates':{}}
    for mode in ('binary','multi'):
        summary['whole_sample'][mode]={'old_label_unique_correct':int(d[mode+'_label_unique_correct'].sum()),
           'old_observed_unique_correct':int(d[mode+'_visible_unique_correct'].sum()),
           'profile_observed_unique_correct':int(p[mode+'_correct'].sum())}
        keep=p[p[mode+'_oracle_Q_U_bound_pass'].eq(True)][KEYS]
        z=d.merge(keep,on=KEYS,validate='one_to_one');obs=z['old_'+mode+'_observed_winners'].str.split('|')
        summary['paired_gates'][mode]={'conditions':len(z),'old_label_unique_correct':int(z[mode+'_label_unique_correct'].sum()),
           'old_observed_unique_correct':int(z[mode+'_visible_unique_correct'].sum()),
           'old_observed_tied':int((obs.map(len)>1).sum()),
           'old_observed_unique_wrong':int(((obs.map(len)==1)&~z[mode+'_visible_unique_correct']).sum()),
           'old_label_correct_but_current':int((z[mode+'_label_unique_correct']&~z[mode+'_visible_unique_correct']).sum()),
           'old_observed_current_only':int(z['old_'+mode+'_observed_winners'].eq('current').sum()),
           'old_observed_symmetric':int(z['old_'+mode+'_observed_winners'].eq('symmetric').sum()),
           'old_observed_opposite_direction':int((z['old_'+mode+'_observed_winners'].isin(['forward','backward'])&~z[mode+'_visible_unique_correct']).sum()),
           'profile_observed_unique_correct':len(z)}
    assert summary['whole_sample']['binary']['old_observed_unique_correct']==371
    assert summary['paired_gates']['binary']['old_observed_unique_correct']==42
    assert summary['whole_sample']['multi']['old_observed_unique_correct']==293
    assert summary['paired_gates']['multi']['old_observed_unique_correct']==39
    aggkeys=[k for k in KEYS if k!='fold'];agg=pd.read_csv(HERE/'margin_aggregate.csv')
    visible=a.groupby(aggkeys).true_family_visible.any().rename('at_least_one_fold_visible').reset_index()
    agg=agg.merge(visible,on=aggkeys,validate='one_to_one');agg=agg[agg.true_kind.ne('current')]
    for mode in ('binary','multi'):agg[mode+'_visible_unique_correct']=agg[mode+'_correct']&agg.at_least_one_fold_visible
    agg[aggkeys+['at_least_one_fold_visible','binary_correct','multi_correct','binary_visible_unique_correct','multi_visible_unique_correct']].to_csv(HERE/'observed_class_aggregate_audit.csv',index=False)
    summary['aggregate_conditions']={'n':len(agg),'binary_label_unique_correct':int(agg.binary_correct.sum()),'binary_visible_unique_correct':int(agg.binary_visible_unique_correct.sum()),'multi_label_unique_correct':int(agg.multi_correct.sum()),'multi_visible_unique_correct':int(agg.multi_visible_unique_correct.sum())}
    (HERE/'observed_class_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
