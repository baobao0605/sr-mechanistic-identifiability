"""Continuous compact-gamma family separation on real reference geometry.

Primary: gamma in [0,.95], same domain as the previous fit grid but searched
continuously.  Known reference, common current nuisance, positive residual
strength r>=rho_min*||B||.  Report both two-sided r-floor and one-sided floor
versus a rival allowed to have zero prediction. They answer different questions.

All bounds use the companion analytic-enclosure solver in float64; no claim of
outward-rounded machine certificates or calibrated probabilistic noise balls.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import sys,json,time
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'theory'))
from continuous_geometry import SpectralSRGeometry

def main():
 meta=json.loads((OUT/'spectral_reference_cases.json').read_text())
 weights=np.load(OUT/'spectral_reference_cases.npz')['weights']
 rows=[];scaled=[];start=time.time()
 for idx,m in enumerate(meta):
  n=m['states'];lam=np.exp(-2j*np.pi*np.arange(n)/n)
  geom=SpectralSRGeometry(lam,weights[m['case_id'],:n])
  # B/S equals F/S by complex conjugation with the real symmetric template.
  result={}
  for a,b in [('F','B'),('F','S')]:
   r=geom.family_margin(a,b,intervals_a=[(0.,.95)],intervals_b=[(0.,.95)],atol=2e-4,max_evals=40001)
   result[a+b]=r
   rows.append({**m,'pair':a+'__'+b,'gamma_low':0.,'gamma_high':.95,**{k:v for k,v in r.items() if k!='gamma_witness'},'gamma_a_witness':r['gamma_witness'][0] if r.get('gamma_witness') else np.nan,'gamma_b_witness':r['gamma_witness'][1] if r.get('gamma_witness') else np.nan})
  # Merged S is one observable class; current has distance rho_min.
  for rho in (.1,.3,1.):
   for two_sided in (False,True):
    lowkey='chord_lower' if two_sided else 'lower'
    upkey='chord_upper' if two_sided else 'upper'
    lo=rho*result['FB'][lowkey];hi=rho*result['FB'][upkey]
    ml=rho*min(1.,result['FB'][lowkey],result['FS'][lowkey])
    mu=rho*min(1.,result['FB'][upkey],result['FS'][upkey])
    scaled.append({**m,'rho_min':rho,'rival_strength_floor':'same_positive_floor' if two_sided else 'zero_allowed','binary_distance_over_reference_lower':lo,'binary_distance_over_reference_upper':hi,'multiway_distance_over_reference_lower':ml,'multiway_distance_over_reference_upper':mu,'binary_max_error_over_reference_sufficient':lo/2,'multiway_max_error_over_reference_sufficient':ml/2,'scope':'known reference; free current profiled; error norm in identical metric; no guarantee these error budgets hold in measurements'})
  if (idx+1)%30==0:print('completed',idx+1,'of',len(meta),'seconds',round(time.time()-start,2),flush=True)
 df=pd.DataFrame(rows);sf=pd.DataFrame(scaled)
 df.to_csv(OUT/'continuous_uniform_geometry.csv',index=False)
 sf.to_csv(OUT/'restricted_theta_separation.csv',index=False)
 primary=sf.query('states==100 and channels==72 and repeats==3')
 primary.to_csv(OUT/'primary_restricted_theta.csv',index=False)
 report={'reference_cases':len(meta),'continuous_family_optimizations':len(df),'restricted_theta_rows':len(sf),'nonconverged':int((df.status!='analytic_bracket_tolerance').sum()),'maximum_ray_width':float((df.upper-df.lower).max()),'maximum_chord_width':float((df.chord_upper-df.chord_lower).max()),'seconds':time.time()-start,'primary_ranges':primary.groupby(['rho_min','rival_strength_floor'])[['binary_distance_over_reference_lower','binary_distance_over_reference_upper','multiway_distance_over_reference_lower','multiway_distance_over_reference_upper']].agg(['min','max']).to_json(),'prior_grid_geometry_is_upper_bound':True,'noise_statement':'Thresholds are allowable deterministic error budgets, not estimated coverage/confidence intervals.'}
 (OUT/'summary.json').write_text(json.dumps(report,indent=2))
 print(json.dumps({k:v for k,v in report.items() if k!='primary_ranges'},indent=2))

if __name__=='__main__':main()
