"""Numerical AUC uncertainty induced by continuous-margin brackets.

These are worst-case ordering bounds, not statistical confidence intervals.
They prevent adaptive-optimizer tolerances from being read as meaningful
ordering when true distances are saturated near one.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

OUT=Path(__file__).resolve().parent

def auc_brackets(y,lo,hi):
    y=np.asarray(y,bool);lo=np.asarray(lo,float);hi=np.asarray(hi,float)
    ok=np.isfinite(lo)&np.isfinite(hi);y=y[ok];lo=lo[ok];hi=hi[ok]
    assert np.all(lo<=hi+1e-12)
    pl,pu=lo[y],hi[y];nl,nu=lo[~y],hi[~y]
    total=len(pl)*len(nl)
    certain_win=np.searchsorted(np.sort(nu),pl,side='left').sum()
    certain_loss=(len(nl)-np.searchsorted(np.sort(nl),pu,side='right')).sum()
    negpoint=nl[nl==nu]
    vals,counts=np.unique(negpoint,return_counts=True)
    count_map=dict(zip(vals.tolist(),counts.tolist()))
    exact_ties=sum(count_map.get(float(v),0) for v in pl[pl==pu])
    return dict(n=len(y),positive=int(y.sum()),
                auc_at_lower=float(roc_auc_score(y,lo)),
                auc_at_upper=float(roc_auc_score(y,hi)),
                auc_ordering_lower=float((certain_win+.5*exact_ties)/total),
                auc_ordering_upper=float((total-certain_loss-.5*exact_ties)/total),
                unresolved_pair_fraction=float(1-(certain_win+certain_loss+exact_ties)/total),
                max_score_bracket_width=float(np.max(hi-lo)))

def main():
    d=pd.read_csv(OUT/'conditional_margin_trials.csv')
    fb=d[d.truth.isin(['forward','backward'])]
    rows=[]
    for scope,out in [('binary','binary_visible_unique_correct'),('multi','correct_family'),('gamma','gamma_correct_0p2')]:
        specs={'unit_geometry':'unit_margin_train_'+scope,
               'I_pair_q_full':'I_train_pair_q_'+scope+'_full',
               'I_pair_q_krylov':'I_train_pair_q_'+scope+'_krylov'}
        for name,col in specs.items():
            rows.append(dict(scope=scope,predictor=name,**auc_brackets(fb[out],fb[col+'_lower'],fb[col+'_upper'])))
        x=fb.predictive_strength_ratio
        rows.append(dict(scope=scope,predictor='strength_only',**auc_brackets(fb[out],x,x)))
    pd.DataFrame(rows).to_csv(OUT/'auc_precision_bounds.csv',index=False)
    stratified=[]
    for gamma,sub in fb.groupby('true_gamma'):
        for name,col in [('unit_geometry','unit_margin_train_gamma'),
                         ('I_pair_q_krylov','I_train_pair_q_gamma_krylov')]:
            stratified.append(dict(true_gamma=gamma,predictor=name,
                accuracy=float(sub.gamma_correct_0p2.mean()),
                **auc_brackets(sub.gamma_correct_0p2,sub[col+'_lower'],sub[col+'_upper'])))
        x=sub.predictive_strength_ratio
        stratified.append(dict(true_gamma=gamma,predictor='strength_only',
            accuracy=float(sub.gamma_correct_0p2.mean()),**auc_brackets(sub.gamma_correct_0p2,x,x)))
    pd.DataFrame(stratified).to_csv(OUT/'gamma_stratified_descriptive.csv',index=False)
    quant=[]
    for scope in ('binary','multi','gamma'):
        for bound in ('lower','upper'):
            x=fb['unit_margin_train_'+scope+'_'+bound]
            quant.append(dict(scope=scope,bound=bound,min=float(x.min()),q10=float(x.quantile(.1)),
                              median=float(x.median()),q90=float(x.quantile(.9)),max=float(x.max())))
    p24=d.query('states==100 and channels==72 and repeats==3 and predictive_strength_ratio==1 and truth in ["forward","backward"]')
    diagnosis={
      'meaning':'Worst-case AUC ordering bounds from analytic margin brackets; not sampling confidence intervals. Correlations and repeated conditions are not independent observations.',
      'unit_margin_quantiles_FB':quant,
      'auc_brackets':rows,
      'primary24':{
          'binary_pair_q_index_min':float(p24.I_train_pair_q_binary_krylov_lower.min()),
          'binary_pair_q_index_max':float(p24.I_train_pair_q_binary_krylov_lower.max()),
          'multi_pair_q_index_min':float(p24.I_train_pair_q_multi_krylov_lower.min()),
          'multi_pair_q_index_max':float(p24.I_train_pair_q_multi_krylov_lower.max()),
          'binary_distance_lower_min':float(p24.unit_margin_train_binary_lower.min()),
          'binary_distance_upper_max':float(p24.unit_margin_train_binary_upper.max())}}
    (OUT/'auc_precision_audit.json').write_text(json.dumps(diagnosis,indent=2))
    print(pd.DataFrame(rows).to_string(index=False))
    print(json.dumps(diagnosis['primary24'],indent=2))

if __name__=='__main__':main()
