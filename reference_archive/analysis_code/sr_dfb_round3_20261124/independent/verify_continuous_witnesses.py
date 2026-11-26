"""Independent original-matrix check of continuous-margin upper witnesses."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from build_inputs import kernel

HERE=Path(__file__).resolve().parent

def residual(B,X):
    out=X-float(np.sum(B*X)/np.sum(B*B))*B
    if np.linalg.norm(out)<=1e-12*max(np.linalg.norm(X),np.finfo(float).tiny):return np.zeros_like(X)
    return out

def main():
    z=np.load(HERE/'repeated_inputs.npz');g=pd.read_csv(HERE/'continuous_geometry.csv');rows=[]
    for _,r in g.iterrows():
        A=z[f'{r.fold_id}_A'];B=z[f'{r.fold_id}_'+('R' if r.phase=='train' else 'S')]
        true=residual(B,kernel(A,r.true_gamma,r.true_kind)@B);true=true/np.linalg.norm(true)
        if not np.isfinite(r.gamma_rival):d=1.
        else:
            kind=r.true_kind if r.rival=='gamma_far' else r.rival
            wrong=residual(B,kernel(A,r.gamma_rival,kind)@B)
            n2=float(np.sum(wrong*wrong))
            a=max(0.,float(np.sum(true*wrong)/n2)) if n2>0 else 0.
            d=float(np.linalg.norm(true-a*wrong))
        diff=abs(d-float(r.upper));lower_ok=float(r.lower)<=d+2e-8
        rows.append({'fold_id':r.fold_id,'phase':r.phase,'true_kind':r.true_kind,'true_gamma':r.true_gamma,'rival':r.rival,'direct_witness_distance':d,'reported_upper':float(r.upper),'reported_lower':float(r.lower),'upper_witness_abs_error':diff,'lower_below_witness':lower_ok})
    d=pd.DataFrame(rows);d.to_csv(HERE/'continuous_witness_checks.csv',index=False)
    summary={'checks':len(d),'upper_witness_max_abs_error':float(d.upper_witness_abs_error.max()),'lower_below_witness_all':bool(d.lower_below_witness.all()),'status':'PASS'}
    assert summary['upper_witness_max_abs_error']<5e-8 and summary['lower_below_witness_all']
    (HERE/'continuous_witness_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
