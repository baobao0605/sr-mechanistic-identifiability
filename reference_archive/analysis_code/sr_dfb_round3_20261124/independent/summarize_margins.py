"""Aggregate the four reused folds and describe predictive utility, without p-values."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
OLD=ROOT/'prior/sr_round2_20261124/math'
KEYS=['dataset','window','space','true_kind','true_gamma','prediction_residual_strength']
SCORE_DECIMALS=12 # Merge arithmetic ties, not genuine geometry differences.

def rss(x):return float(np.sqrt(np.sum(np.asarray(x,float)**2)))
def auc(y,s):
    ok=np.isfinite(np.asarray(s,float))
    y=np.asarray(y,bool)[ok];s=np.round(np.asarray(s,float)[ok],SCORE_DECIMALS)
    return float(roc_auc_score(y,s)) if len(np.unique(y))==2 else float('nan')

def main():
    d=pd.read_csv(HERE/'margin_folds.csv');old=pd.read_csv(OLD/'calibration_selection.csv');ag=pd.read_csv(OLD/'calibration_aggregate.csv')
    out=[]
    for key,g in d.groupby(KEYS,dropna=False):
        r=dict(zip(KEYS,key));
        o=old.copy();a=ag.copy()
        for k,v in r.items():o=o[o[k].eq(v)];a=a[a[k].eq(v)]
        assert len(o)==1 and len(g)==4
        r['multi_correct']=bool(o.unique_correct_family.iloc[0]);r['multi_winners']=o.best_class_set.iloc[0]
        sub=a[a.candidate.isin(['forward','backward'])]
        eps=1e-10*float(sub.test_target_norm2.iloc[0]);best=float(sub.test_sse.min())
        names=sub.loc[sub.test_sse<=best+eps,'candidate'].tolist()
        r['binary_correct']=len(names)==1 and names[0]==r['true_kind'];r['binary_winners']='|'.join(sorted(names))
        for name in ['train_current_residual_strength','test_current_residual_strength','train_transport_error_full','train_transport_error_U','reference_repeat_error_full','reference_repeat_error_U','oracle_test_error_full','oracle_test_error_U','correct_fit_error_train','correct_fit_error_test','correct_fit_test_residual_U','train_reference_norm']:
            r[name]=rss(g[name])
        r['test_error_U_fraction']=r['oracle_test_error_U']/r['oracle_test_error_full']
        r['gamma_reportable_folds']=int(g.gamma_reportable.sum())
        r['gamma_close_folds']=int(g.gamma_error_lt_point2.sum())
        if r['true_kind']!='current':
            for phase in ('train','test'):
                for rival in ('binary','symmetric','current','gamma_far'):
                    for bound in ('lower','upper'):
                        r[f'{phase}_{rival}_margin_{bound}']=rss(g[f'{phase}_{rival}_margin_{bound}'])
                for bound in ('lower','upper'):
                    r[f'{phase}_multi_margin_{bound}']=min(r[f'{phase}_{k}_margin_{bound}'] for k in ('binary','symmetric','current'))
            for mode in ('binary','multi'):
                q=r['correct_fit_error_test'];dm=r[f'test_{mode}_margin_lower']
                r[f'oracle_{mode}_full_bound_pass']=dm>2*r['oracle_test_error_full']+q+1e-10
                r[f'oracle_{mode}_U_bound_pass']=dm>2*r['oracle_test_error_U']+q+1e-10
                r[f'oracle_{mode}_U_tight_pass']=dm>r['oracle_test_error_U']+r['correct_fit_test_residual_U']+1e-10
                for s in ('transport','repeat'):
                    e=r['train_transport_error_U' if s=='transport' else 'reference_repeat_error_U']
                    r[f'pretest_{mode}_{s}_ratio']=r[f'train_{mode}_margin_lower']/max(2*e+r['correct_fit_error_train'],1e-300)
                r[f'pretest_{mode}_unit_geometry']=r[f'train_{mode}_margin_lower']/max(r['train_current_residual_strength'],1e-300)
                r[f'pretest_{mode}_relative_margin']=r[f'train_{mode}_margin_lower']/r['train_reference_norm']
        out.append(r)
    z=pd.DataFrame(out);z.to_csv(HERE/'margin_aggregate.csv',index=False)
    x=z[z.true_kind!='current'];des=[]
    for group,sub in [('all',x)]+[(ds,x[x.dataset.eq(ds)]) for ds in ('Kern','Turner')]:
        for mode in ('binary','multi'):
            metrics={'unit_geometry':f'pretest_{mode}_unit_geometry','relative_margin':f'pretest_{mode}_relative_margin','transport_ratio':f'pretest_{mode}_transport_ratio','repeat_ratio':f'pretest_{mode}_repeat_ratio','prediction_strength':'prediction_residual_strength'}
            for metric,col in metrics.items():
                des.append({'group':group,'mode':mode,'metric':metric,'n_conditions':len(sub),'n_correct':int(sub[f'{mode}_correct'].sum()),'auroc':auc(sub[f'{mode}_correct'],sub[col]),'within_strength_macro_auroc':float(np.nanmean([auc(g[f'{mode}_correct'],g[col]) for _,g in sub.groupby('prediction_residual_strength')]))})
    pd.DataFrame(des).to_csv(HERE/'predictive_auroc.csv',index=False)
    # Choose one high-score acceptance threshold from one entire dataset,
    # then evaluate it on the other. No thresholds use held-out labels.
    cross=[]
    for train_ds,test_ds in [('Kern','Turner'),('Turner','Kern')]:
        tr=x[x.dataset.eq(train_ds)];te=x[x.dataset.eq(test_ds)]
        for mode in ('binary','multi'):
            for metric,col in [('unit_geometry',f'pretest_{mode}_unit_geometry'),('relative_margin',f'pretest_{mode}_relative_margin'),('transport_ratio',f'pretest_{mode}_transport_ratio'),('repeat_ratio',f'pretest_{mode}_repeat_ratio'),('prediction_strength','prediction_residual_strength')]:
                options=[]
                training_score=np.round(tr[col],SCORE_DECIMALS)
                validation_score=np.round(te[col],SCORE_DECIMALS)
                for threshold in np.unique(training_score):
                    sel=training_score>=threshold;n=int(sel.sum())
                    if n>=10 and float(tr.loc[sel,f'{mode}_correct'].mean())>=.9:options.append((n,float(threshold)))
                if options:
                    n,th=max(options,key=lambda v:(v[0],-v[1]));a=training_score>=th;b=validation_score>=th
                    cross.append({'train_dataset':train_ds,'test_dataset':test_ds,'mode':mode,'metric':metric,'status':'training_threshold_found','threshold':th,'train_selected':int(a.sum()),'train_correct':int(tr.loc[a,f'{mode}_correct'].sum()),'test_selected':int(b.sum()),'test_correct':int(te.loc[b,f'{mode}_correct'].sum()),'test_precision':float(te.loc[b,f'{mode}_correct'].mean()) if b.any() else np.nan})
                else:cross.append({'train_dataset':train_ds,'test_dataset':test_ds,'mode':mode,'metric':metric,'status':'no_training_threshold_with_90pct_and_10_conditions','threshold':np.nan,'train_selected':0,'train_correct':0,'test_selected':0,'test_correct':0,'test_precision':np.nan})
    pd.DataFrame(cross).to_csv(HERE/'cross_dataset_thresholds.csv',index=False)
    by=[]
    for (ds,win,space),g in x.groupby(['dataset','window','space']):
        r={'dataset':ds,'window':win,'space':space,'n':len(g),'binary_correct':int(g.binary_correct.sum()),'multi_correct':int(g.multi_correct.sum()),'median_projected_error_fraction':float(g.test_error_U_fraction.median())}
        for mode in ('binary','multi'):
            for suffix in ('full_bound','U_bound','U_tight'):
                a=g[f'oracle_{mode}_{suffix}_pass'].astype(bool)
                r[f'oracle_{mode}_{suffix}_selected']=int(a.sum());r[f'oracle_{mode}_{suffix}_wrong']=int((a&~g[f'{mode}_correct']).sum())
            for scale in ('transport','repeat'):
                a=g[f'pretest_{mode}_{scale}_ratio']>1;r[f'pretest_{mode}_{scale}_selected']=int(a.sum());r[f'pretest_{mode}_{scale}_wrong']=int((a&~g[f'{mode}_correct']).sum())
        by.append(r)
    pd.DataFrame(by).to_csv(HERE/'group_summary.csv',index=False)
    summary={'conditions_total':len(z),'direction_conditions':len(x),'binary_correct':int(x.binary_correct.sum()),'multi_correct':int(x.multi_correct.sum()),'oracle':{},'cross_dataset_thresholds':cross,'auroc':des,'score_rounding_decimals':SCORE_DECIMALS,'scope':'Descriptive across reused calibration conditions; no population p values; threshold selection uses only source dataset. Arithmetic score ties are merged at 12 decimal places.'}
    for mode in ('binary','multi'):
        for suf in ('full_bound','U_bound','U_tight'):
            a=x[f'oracle_{mode}_{suf}_pass'].astype(bool);assert x.loc[a,f'{mode}_correct'].all()
            summary['oracle'][f'{mode}_{suf}']={'selected':int(a.sum()),'wrong':int((a&~x[f'{mode}_correct']).sum())}
    (HERE/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=True))
    print(json.dumps({k:v for k,v in summary.items() if k not in ('cross_dataset_thresholds','auroc')},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
