"""Bound descriptive AUROC uncertainty from continuous margin brackets.

These are worst-case pairwise ranking bounds, not sampling confidence
intervals. No inference about independent participants is performed.
"""
from pathlib import Path
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
def bounds(y,lo,hi,point=False):
    y=np.asarray(y,bool);lo=np.asarray(lo,float);hi=np.asarray(hi,float)
    if point:
        d=lo[y,None]-lo[None,~y]
        a=float((d>0).mean()+.5*(d==0).mean());return a,a
    pl=lo[y,None];ph=hi[y,None];nl=lo[None,~y];nh=hi[None,~y]
    # Touching/overlapping intervals may produce a tie; assigning their
    # possible contribution the wider [0,1] remains a conservative bound.
    lower=float((pl>nh+1e-12).mean());upper=float((ph>=nl-1e-12).mean())
    return lower,upper
def main():
    d=pd.read_csv(HERE/'margin_aggregate.csv');d=d[d.true_kind.ne('current')]
    auc=pd.read_csv(HERE/'predictive_auroc.csv');rows=[]
    for group,x in [('all',d)]+[(s,d[d.dataset.eq(s)]) for s in ('Kern','Turner')]:
      for mode in ('binary','multi'):
        lo=x['train_'+mode+'_margin_lower'].to_numpy();hi=x['train_'+mode+'_margin_upper'].to_numpy()
        metrics={
          'unit_geometry':(lo/x.train_current_residual_strength,hi/x.train_current_residual_strength),
          'relative_margin':(lo/x.train_reference_norm,hi/x.train_reference_norm),
          'transport_ratio':(lo/(2*x.train_transport_error_U+x.correct_fit_error_train),hi/(2*x.train_transport_error_U+x.correct_fit_error_train)),
          'repeat_ratio':(lo/(2*x.reference_repeat_error_U+x.correct_fit_error_train),hi/(2*x.reference_repeat_error_U+x.correct_fit_error_train)),
          'prediction_strength':(x.prediction_residual_strength,x.prediction_residual_strength)}
        for name,(l,h) in metrics.items():
            low,high=bounds(x[mode+'_correct'],l,h,point=name=='prediction_strength')
            a=float(auc.loc[auc.group.eq(group)&auc['mode'].eq(mode)&auc.metric.eq(name),'auroc'].iloc[0])
            assert low-1e-10<=a<=high+1e-10
            rows.append({'group':group,'mode':mode,'metric':name,'n_conditions':len(x),'observed_score_auc':a,
                         'auc_lower_from_geometry_brackets':low,'auc_upper_from_geometry_brackets':high,
                         'scope':'Conservative numerical-ranking interval; not a sampling confidence interval'})
    pd.DataFrame(rows).to_csv(HERE/'auc_precision_bounds.csv',index=False)
    print(pd.DataFrame(rows).query("group != 'all'").to_string(index=False))
if __name__=='__main__':main()
