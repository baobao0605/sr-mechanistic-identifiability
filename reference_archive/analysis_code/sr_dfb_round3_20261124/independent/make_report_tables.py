"""Assemble paired decision-rule and discount-resolution reporting tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
def main():
    f=pd.read_csv(HERE/'margin_folds.csv');p=pd.read_csv(HERE/'profiled_decision_folds.csv')
    keys=['dataset','window','fold','space','true_kind','true_gamma','prediction_residual_strength','fold_id']
    cols=['binary_correct','multi_correct','binary_winners','multi_winners','correct_fit_error_test',
          'correct_fit_test_residual_U','oracle_binary_U_bound_pass','oracle_multi_U_bound_pass',
          'oracle_binary_U_tight_pass','oracle_multi_U_tight_pass','true_alpha']
    pair=p.merge(f[keys+cols],on=keys,validate='one_to_one',suffixes=('_profile','_frozen'))
    audit=pd.read_csv(HERE/'observed_class_audit.csv')
    auditcols=['true_family_fitted_predictive_gain','true_family_visible','old_binary_observed_winners','old_multi_observed_winners','binary_visible_unique_correct','multi_visible_unique_correct']
    pair=pair.merge(audit[keys+auditcols],on=keys,validate='one_to_one')
    rows=[]
    for mode in ('binary','multi'):
        d=pair[pair.true_kind.ne('current') & pair[mode+'_oracle_Q_U_bound_pass'].fillna(False).astype(bool)].copy()
        d['gate_slack']=d[mode+'_margin_lower']-2*d.oracle_error_Q_U
        for rule in ('frozen','profile'):
            label=d[mode+'_winners_'+rule].str.split('|')
            d[rule+'_tied']=label.map(len)>1
            d[rule+'_unique_wrong']=(label.map(len)==1)&~d[mode+'_correct_'+rule]
            d[rule+'_tied_including_truth']=[len(v)>1 and t in v for v,t in zip(label,d.true_kind)]
            d[rule+'_tied_excluding_truth']=[len(v)>1 and t not in v for v,t in zip(label,d.true_kind)]
        d['scope']='Target-level reprofiled nearest family; old frozen training-fit correctness shown only as decision-rule ablation'
        d.to_csv(HERE/('profiled_'+mode+'_threshold_cases.csv'),index=False)
        for ds,g in [('all',d)]+[(s,d[d.dataset.eq(s)]) for s in ('Kern','Turner')]:
            rows.append({'mode':mode,'dataset':ds,'conditions':len(g),
                         'frozen_label_unique_correct':int(g[mode+'_correct_frozen'].sum()),
                         'frozen_label_not_uniquely_recovered':int((~g[mode+'_correct_frozen']).sum()),
                         'frozen_label_tied':int(g.frozen_tied.sum()),
                         'frozen_label_tied_including_truth':int(g.frozen_tied_including_truth.sum()),
                         'frozen_label_tied_excluding_truth':int(g.frozen_tied_excluding_truth.sum()),
                         'frozen_label_unique_wrong':int(g.frozen_unique_wrong.sum()),
                         'frozen_observed_unique_correct':int(g[mode+'_visible_unique_correct'].sum()),
                         'frozen_observed_current_only':int(g['old_'+mode+'_observed_winners'].eq('current').sum()),
                         'frozen_observed_symmetric':int(g['old_'+mode+'_observed_winners'].eq('symmetric').sum()),
                         'frozen_observed_tied':int(g['old_'+mode+'_observed_winners'].str.contains('|',regex=False).sum()),
                         'profile_observed_unique_correct':int(g[mode+'_correct_profile'].sum()),
                         'profile_not_uniquely_recovered':int((~g[mode+'_correct_profile']).sum()),
                         'profile_tied':int(g.profile_tied.sum()),
                         'profile_unique_wrong':int(g.profile_unique_wrong.sum()),
                         'minimum_margin_minus_2error':float(g.gate_slack.min())})
    pd.DataFrame(rows).to_csv(HERE/'paired_gate_summary.csv',index=False)
    (HERE/'paired_gate_summary.json').write_text(json.dumps(rows,indent=2))
    g=pd.read_csv(HERE/'continuous_geometry.csv');g=g[g.rival.eq('gamma_far')]
    rows=[]
    for (ds,space),dd in f[f.true_kind.ne('current')].groupby(['dataset','space']):
        gg=g[g.dataset.eq(ds)&g.space.eq(space)]
        rows.append({'dataset':ds,'space':space,'n_direction_fold_conditions':len(dd),
                     'unit_margin_lower_min':float(gg.lower.min()),'unit_margin_lower_median':float(gg.lower.median()),
                     'unit_margin_lower_max':float(gg.lower.max()),'unit_margin_upper_max':float(gg.upper.max()),
                     'old_gamma_reportable':int(dd.gamma_reportable.sum()),'old_gamma_error_lt_point2':int(dd.gamma_error_lt_point2.sum()),
                     'single_plane_exact_zero':space in ('design_k4_2D','Fourier1_2D')})
    pd.DataFrame(rows).to_csv(HERE/'discount_resolution_summary.csv',index=False)
    rows=[]
    for ds,dd in [('all',p[p.true_kind.ne('current')])]+[(ds,p[p.true_kind.ne('current')&p.dataset.eq(ds)]) for ds in ('Kern','Turner')]:
        rows.append({'dataset':ds,'n':len(dd),'median_U_error_over_full':float((dd.oracle_error_U/dd.oracle_error_full).median()),
                     'median_Q_U_error_over_full':float((dd.oracle_error_Q_U/dd.oracle_error_full).median())})
    pd.DataFrame(rows).to_csv(HERE/'projection_ratio_summary.csv',index=False)

if __name__=='__main__':main()
