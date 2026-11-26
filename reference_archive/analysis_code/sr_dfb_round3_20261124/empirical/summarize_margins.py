"""Join continuous brackets to frozen recovery and summarize declared tests."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
ROOT=Path(__file__).resolve().parents[1];OUT=Path(__file__).resolve().parent
KEY=['states','channels','repeats','fold','anchor','truth','true_gamma','rival']

def main():
    d=pd.read_csv(OUT/'conditional_margin_trials_grid.csv')
    # A02: the rival domain includes |gamma - gamma_true| >= 0.2.
    # Grid endpoints such as 0.7 - 0.5 must not become successes merely
    # because their binary floating-point subtraction is just below 0.2.
    if 'gamma_correct_0p2_legacy_float' not in d:
        d['gamma_correct_0p2_legacy_float']=d['gamma_correct_0p2']
    d['gamma_correct_0p2']=d['gamma_error_conditional'].notna() & (d['gamma_error_conditional']<.2-1e-12)
    g=pd.read_csv(OUT/'point_geometry_continuous.csv').set_index(KEY)
    rows=[]
    for r in d.to_dict('records'):
      if r['truth']!='current':
        rivals=[k for k in ('forward','backward','symmetric_transition') if k!=r['truth']]
        for anchor in ('train','test'):
          name='train' if anchor=='train' else 'test_oracle'
          scale=r['true_prediction_norm_'+anchor]
          points={}
          for rv in rivals+['gamma_0p2']:
            q=g.loc[(r['states'],r['channels'],r['repeats'],r['fold'],name,r['truth'],r['true_gamma'],rv)]
            points[rv]=(q.lower,q.upper)
            r['margin_'+anchor+'_to_'+rv+'_lower']=scale*q.lower
            r['margin_'+anchor+'_to_'+rv+'_upper']=scale*q.upper
          points['current']=(1.,1.)
          for scope,kk in [('binary',[('backward' if r['truth']=='forward' else 'forward')] if r['truth'] in ('forward','backward') else []),('multi',rivals+['current']),('gamma',['gamma_0p2'])]:
            if not kk:continue
            lo=min(points[k][0] for k in kk);hi=min(points[k][1] for k in kk)
            r['unit_margin_'+anchor+'_'+scope+'_lower']=lo
            r['unit_margin_'+anchor+'_'+scope+'_upper']=hi
            r['margin_'+anchor+'_'+scope+'_lower']=lo*scale
            r['margin_'+anchor+'_'+scope+'_upper']=hi*scale
            if anchor=='train':
              for sp in ('full','krylov'):
                e=r['repeat_error_train_'+sp]/np.sqrt(2)
                r['I_train_'+scope+'_'+sp+'_lower']=scale*lo/max(2*e,1e-30)
                r['I_train_'+scope+'_'+sp+'_upper']=scale*hi/max(2*e,1e-30)
                epair=r['repeat_error_train_'+sp]
                qtrain=r['fit_parameter_shift_train']
                for budget_name,budget in [('pair',2*epair),('pair_q',2*epair+qtrain)]:
                  r['I_train_'+budget_name+'_'+scope+'_'+sp+'_lower']=scale*lo/max(budget,1e-30)
                  r['I_train_'+budget_name+'_'+scope+'_'+sp+'_upper']=scale*hi/max(budget,1e-30)
            else:
              for sp in ('full','krylov'):
                e=r['actual_error_test_'+sp];rf=r['fit_error_test_true_'+sp];h=r['fit_parameter_shift_test_'+sp]
                budget=e+rf;tol=1e-8*max(1.,budget,scale*hi)
                r['oracle_'+scope+'_'+sp+'_budget']=budget
                r['oracle_'+scope+'_'+sp+'_index_lower']=scale*lo/max(budget,1e-30)
                r['oracle_'+scope+'_'+sp+'_guarantee']=bool(scale*lo>budget+tol)
                r['oracle_'+scope+'_'+sp+'_looser_guarantee']=bool(scale*lo>2*e+h+tol)
        # Same subspace for all hypotheses: deterministic triangle bound.
        r['min_test_wrong_residual_full']=min(r['test_error_'+k+'_full'] for k in rivals+['current'])
        r['min_test_wrong_residual_krylov']=min(r['test_error_'+k+'_krylov'] for k in rivals+['current'])
      rows.append(r)
    d=pd.DataFrame(rows)
    d.to_csv(OUT/'conditional_margin_trials.csv',index=False)
    primary=d.query('states==100 and channels==72 and repeats==3 and predictive_strength_ratio==1 and truth in ["forward","backward"]')
    primary.to_csv(OUT/'primary24.csv',index=False)
    fb=d.query('truth in ["forward","backward"]').copy();pred=d.query('truth!="current"').copy()
    stats=[];risk=[];oracles=[]
    for cohort,sub in [('FB',fb),('all_predictive',pred)]:
      specs=[('multi','correct_family'),('gamma','gamma_correct_0p2')]
      if cohort=='FB':specs.insert(0,('binary','binary_visible_unique_correct'))
      for scope,outcome in specs:
        for partition,ss in [('pooled',sub)]+[(f'fold_{f}',sub[sub.fold==f]) for f in range(4)]:
          y=ss[outcome].astype(bool)
          predictors={'unit_geometry':ss['unit_margin_train_'+scope+'_lower'],
                      'strength_only':ss.predictive_strength_ratio,
                      'strength_times_geometry':ss.predictive_strength_ratio*ss['unit_margin_train_'+scope+'_lower'],
                      'I_naive_full':ss['I_train_'+scope+'_full_lower'],
                      'I_naive_krylov':ss['I_train_'+scope+'_krylov_lower'],
                      'I_pair_full':ss['I_train_pair_'+scope+'_full_lower'],
                      'I_pair_krylov':ss['I_train_pair_'+scope+'_krylov_lower'],
                      'I_pair_q_full':ss['I_train_pair_q_'+scope+'_full_lower'],
                      'I_pair_q_krylov':ss['I_train_pair_q_'+scope+'_krylov_lower']}
          for label,x in predictors.items():
            valid=x.notna();xx=x[valid];yy=y[valid]
            passing=xx>1+1e-8
            stats.append(dict(cohort=cohort,scope=scope,partition=partition,predictor=label,n=len(xx),
                              correct=int(yy.sum()),accuracy=float(yy.mean()),
                              descriptive_auc=float(roc_auc_score(yy,xx)) if yy.nunique()==2 else np.nan,
                              threshold_1_selected=int(passing.sum()) if label.startswith('I_') else np.nan,
                              threshold_1_correct=int(yy[passing].sum()) if label.startswith('I_') else np.nan))
            if partition=='pooled':
              for cover in (.05,.1,.2,.3,.5,.75,1.):
                cutoff=float(xx.quantile(1-cover));sel=xx>=cutoff
                risk.append(dict(cohort=cohort,scope=scope,predictor=label,requested_coverage=cover,
                                 cutoff=cutoff,n=len(xx),selected=int(sel.sum()),coverage=float(sel.mean()),
                                 wrong=int((~yy[sel]).sum()),error_rate=float((~yy[sel]).mean())))
        for sp in ('full','krylov'):
          selected=sub['oracle_'+scope+'_'+sp+'_guarantee'].fillna(False).astype(bool)
          correct=sub[outcome].astype(bool)
          oracles.append(dict(cohort=cohort,scope=scope,space=sp,n=len(sub),
                              selected=int(selected.sum()),correct=int(correct[selected].sum()),
                              wrong=int((~correct[selected]).sum()),
                              looser_selected=int(sub['oracle_'+scope+'_'+sp+'_looser_guarantee'].fillna(False).sum()),
                              max_index=float(sub['oracle_'+scope+'_'+sp+'_index_lower'].max())))
          assert not ((~correct)&selected).any(),(cohort,scope,sp)
    pd.DataFrame(stats).to_csv(OUT/'predictor_performance.csv',index=False)
    pd.DataFrame(risk).to_csv(OUT/'risk_coverage_descriptive.csv',index=False)
    pd.DataFrame(oracles).to_csv(OUT/'oracle_inequality_results.csv',index=False)
    selected=fb[fb.I_train_binary_krylov_lower>1+1e-8]
    selected.to_csv(OUT/'empirical_threshold_failures.csv',index=False)
    # Validate the entire-family geometric lower bound against actual frozen
    # wrong fitted residuals, not only the selected cases.
    checks=[]
    for sp in ('full','krylov'):
      residual=pred['min_test_wrong_residual_'+sp]
      lower=pred['margin_test_multi_lower']-pred['actual_error_test_'+sp]
      violation=float((lower-residual).max())
      checks.append({'space':sp,'n':len(pred),'max_lower_minus_observed_wrong_residual':violation})
      assert violation<1e-7
    perf=pd.DataFrame(stats)
    summary={'n_trials':len(d),'n_FB':len(fb),'n_all_predictive':len(pred),
             'primary24':{'binary_forced_correct':int(primary.binary_forced_correct.sum()),
                          'binary_visible_unique_correct':int(primary.binary_visible_unique_correct.sum()),
                          'multi_unique_correct':int(primary.correct_family.sum()),
                          'primary_I_pair_q_binary_krylov_above_1':int((primary.I_train_pair_q_binary_krylov_lower>1).sum()),
                          'primary_I_pair_q_multi_krylov_above_1':int((primary.I_train_pair_q_multi_krylov_lower>1).sum()),
                          'primary_I_pair_q_binary_krylov_min':float(primary.I_train_pair_q_binary_krylov_lower.min()),
                          'primary_I_pair_q_binary_krylov_max':float(primary.I_train_pair_q_binary_krylov_lower.max()),
                          'primary_I_pair_q_multi_krylov_min':float(primary.I_train_pair_q_multi_krylov_lower.min()),
                          'primary_I_pair_q_multi_krylov_max':float(primary.I_train_pair_q_multi_krylov_lower.max()),
                          'naive_mismatched_I_binary_krylov_min':float(primary.I_train_binary_krylov_lower.min()),
                          'naive_mismatched_I_binary_krylov_max':float(primary.I_train_binary_krylov_lower.max()),
                          'unit_binary_margin_min':float(primary.unit_margin_train_binary_lower.min()),
                          'unit_binary_margin_max':float(primary.unit_margin_train_binary_lower.max())},
             'FB_pooled_performance':perf.query('cohort=="FB" and partition=="pooled"').to_dict('records'),
             'oracle_results':oracles,'triangle_inequality_checks':checks,
             'naive_mismatched_budget_control':{'selected':len(selected),'correct':int(selected.binary_visible_unique_correct.sum()),
                                               'n_states':sorted(selected.states.unique().tolist()),
                                               'interpretation':'This obsolete index substitutes pair difference / sqrt(2) for the actual pair deviation and omits fitted-parameter error. Its failures do not test the deterministic separation theorem.'},
             'primary_index':'I_pair_q, amended before the disjoint 6+6 experiment; naive single-repeat index retained as an explicitly mismatched comparison',
             'meaning':'All outcome counts are conditional computations on one participant, not independent participants. Empirical repeat radius has no coverage guarantee. Oracle inequalities use known imposed transform and heldout data.'}
    (OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
