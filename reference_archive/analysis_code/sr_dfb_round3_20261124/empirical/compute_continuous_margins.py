"""Continuous gamma brackets for every THINGS conditional shape.

Uses analytic interval bounds supplied by the independent theory branch.
The intervals are evaluated in float64, not outward-rounded arithmetic.
The cache references are actual training-whitened repeat averages.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,sys,time,hashlib
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'theory'))
from continuous_geometry import SpectralSRGeometry
from analyze_margin_inputs import Orbit,TG,KINDS

def work(design):
    n,nc,nr,fold=design
    dd=dict(states=n,channels=nc,repeats=nr,fold=fold)
    key=f'n{n}_p{nc}_r{nr}_f{fold}'
    points=[];families=[]
    with np.load(OUT/'repeat_inputs.npz') as z:
      for letter,anchor in [('A','train'),('C','test_oracle')]:
        x=z[key+'_'+letter];orbit=Orbit(x)
        geom=SpectralSRGeometry(np.exp(-2j*np.pi*np.arange(n)/n),orbit.w)
        for kind in KINDS[:3]:
          for gamma in TG:
            for rival in [k for k in KINDS[:3] if k!=kind]+['gamma_0p2']:
              ints=[(0.,.95)]
              if rival=='gamma_0p2':
                ints=[]
                if gamma-.2>=0:ints.append((0.,gamma-.2))
                if gamma+.2<=.95:ints.append((gamma+.2,.95))
              r=geom.point_margin(kind,gamma,kind if rival=='gamma_0p2' else rival,
                                  intervals=ints,atol=1e-4,max_evals=4097)
              points.append({**dd,'anchor':anchor,'truth':kind,'true_gamma':gamma,'rival':rival,
                             **r,'true_residual_norm':geom.residual_norm(kind,gamma),
                             'intervals':json.dumps(ints),'requested_atol':1e-4,'max_evaluations':4097})
        if anchor=='train':
          for ka,kb in [('forward','backward'),('forward','symmetric_transition'),('backward','symmetric_transition')]:
            r=geom.family_margin(ka,kb,atol=5e-4,max_evals=20001)
            witness=r.pop('gamma_witness',None)
            families.append({**dd,'kind_a':ka,'kind_b':kb,**r,
                             'gamma_a_witness':None if witness is None else witness[0],
                             'gamma_b_witness':None if witness is None else witness[1],
                             'requested_atol':5e-4,'max_evaluations':20001})
    return points,families

def main():
    old=pd.read_csv(OUT/'conditional_margin_trials_grid.csv')
    designs=list(old[['states','channels','repeats','fold']].drop_duplicates().itertuples(index=False,name=None))
    start=time.time();points=[];families=[]
    checkpoint=OUT/'continuous_checkpoint.jsonl'
    with checkpoint.open('w') as stream,ProcessPoolExecutor(max_workers=4) as pool:
      jobs={pool.submit(work,d):d for d in designs}
      for i,f in enumerate(as_completed(jobs),1):
        p,g=f.result();points+=p;families+=g
        stream.write(json.dumps({'design':jobs[f],'points':p,'families':g})+'\n');stream.flush()
        if i%20==0:print(i,'/',len(designs),'seconds',round(time.time()-start,2),flush=True)
    p=pd.DataFrame(points).sort_values(['states','channels','repeats','fold','anchor','truth','true_gamma','rival'])
    f=pd.DataFrame(families).sort_values(['states','channels','repeats','fold','kind_a','kind_b'])
    p.to_csv(OUT/'point_geometry_continuous.csv',index=False);f.to_csv(OUT/'family_geometry_continuous.csv',index=False)
    pg=pd.read_csv(OUT/'point_geometry_grid.csv')
    keys=['states','channels','repeats','fold','anchor','truth','true_gamma','rival']
    j=p.merge(pg,on=keys,suffixes=('','_grid'),validate='one_to_one')
    lower_violation=(j.lower-j.unit_margin_grid_upper).max()
    residual_norm_difference=abs(j.true_residual_norm-j.true_residual_norm_grid).max()
    assert lower_violation<1e-7 and residual_norm_difference<1e-8
    fg=pd.read_csv(OUT/'family_geometry_grid.csv')
    jf=f.merge(fg,on=['states','channels','repeats','fold','kind_a','kind_b'],validate='one_to_one')
    assert (jf.lower-jf.unit_margin_grid_upper).max()<1e-7
    report={'status':'complete','designs':len(designs),'points':len(p),'families':len(f),
            'seconds':time.time()-start,'point_status':p.status.value_counts().to_dict(),
            'family_status':f.status.value_counts().to_dict(),
            'point_max_gap':float((p.upper-p.lower).max()),
            'family_max_ray_gap':float((f.upper-f.lower).max()),
            'family_max_chord_gap':float((f.chord_upper-f.chord_lower).max()),
            'lower_minus_grid_max':float(lower_violation),
            'residual_norm_max_difference':float(residual_norm_difference),
            'numeric_scope':'analytic interval bounds evaluated in float64; not outward-rounded machine certificates',
            'theory_implementation_sha256':hashlib.sha256((ROOT/'theory/continuous_geometry.py').read_bytes()).hexdigest()}
    (OUT/'continuous_geometry_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
