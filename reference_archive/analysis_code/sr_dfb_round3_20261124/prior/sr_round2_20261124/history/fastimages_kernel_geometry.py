from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.linalg import svd
from threadpoolctl import threadpool_limits
threadpool_limits(2)
P=Path(__file__).resolve().parent
z=np.load(P/'fastimages_kernel_input.npz',allow_pickle=True)
E=z['localizer_epochs']; lab=z['localizer_labels']; blocks=z['localizer_blocks']; rows=[]; time_rows=[]
for fold in range(8):
 tr=blocks!=fold;te=blocks==fold
 Htr=np.array([E[tr&(lab==j)].mean(0) for j in range(5)])
 Hte=np.array([E[te&(lab==j)].mean(0) for j in range(5)])
 res=E[tr]-Htr[lab[tr]];sc=np.sqrt(np.mean(res*res,axis=(0,1)));sc=np.maximum(sc,np.median(sc)*1e-6)
 A=(Htr-Htr.mean(0))/sc;B=(Hte-Hte.mean(0))/sc
 # Time by combined stimulus/sensor axes: rank-one is a common temporal kernel.
 A=A.transpose(1,0,2).reshape(101,-1); B=B.transpose(1,0,2).reshape(101,-1)
 u,s,v=svd(A,full_matrices=False,check_finite=False)
 for q in [0,1,2,4,8,16,32,101]:
  pred=(u[:,:q]*s[:q])@v[:q] if q else np.zeros_like(A)
  rows.append({'fold':fold,'rank':q,'SSE':float(np.sum((B-pred)**2)),'test_energy':float(np.sum(B*B)),'training_energy_fraction':float(np.sum(s[:q]**2)/np.sum(s*s))})
 reference=A[10:25].mean(0);den=reference@reference
 for t in range(101):time_rows.append({'fold':fold,'latency_ms':t*10,'training_pattern_projection':float(A[t]@reference/den),'heldout_pattern_projection':float(B[t]@reference/den)})
r=pd.DataFrame(rows);r.to_csv(P/'fastimages_temporal_rank_folds.csv',index=False)
s=r.groupby('rank',as_index=False)[['SSE','test_energy']].sum();s['R2_to_zero_contrast']=1-s.SSE/s.test_energy;s.to_csv(P/'fastimages_temporal_rank_summary.csv',index=False)
tm=pd.DataFrame(time_rows);tm.to_csv(P/'fastimages_temporal_projection.csv',index=False)
print(s.to_string(index=False))
g=pd.read_csv(P/'fastimages_kernel_discrimination_geometry.csv');print(g.groupby('gamma').agg({'remaining_distance_fraction':['mean','min','max'],'removed_information_fraction':'mean'}))
(P/'fastimages_geometry_interpretation.json').write_text(json.dumps({'separability':'Rank-one across time and stimulus-sensor axes is the shared scalar kernel h(t) B assumption. Cross-block mean-pattern prediction evaluates it without selecting held-out temporal axes. Best low-rank performance is descriptive across stated ranks, not independently tuned model selection.','projection':'Projected coefficient versus independently estimated current template is not a literal probability or discounted future occupancy. Its sign and non-monotonicity are valid descriptions of the observed neural readout only.','nuisance':'Distance removed by richer physical-history nuisance is a linear observation-space calculation based on measured localizer kernels. It is neither proof that nuisance occurs with that strength nor a statistical confidence bound.'},indent=2))
