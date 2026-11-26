from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.linalg import svd
from threadpoolctl import threadpool_limits
threadpool_limits(2)
P=Path(__file__).resolve().parent
z=np.load(P/'fastimages_kernel_input.npz',allow_pickle=True);E=z['localizer_epochs'];lab=z['localizer_labels'];groups=z['localizer_blocks'];qs=[0,1,2,4,8,16,32,101]

def prepare(tr,te):
 H=np.array([E[tr&(lab==j)].mean(0) for j in range(5)]);V=np.array([E[te&(lab==j)].mean(0) for j in range(5)])
 R=E[tr]-H[lab[tr]];sc=np.sqrt(np.mean(R*R,axis=(0,1)));sc=np.maximum(sc,np.median(sc)*1e-6)
 A=((H-H.mean(0))/sc).transpose(1,0,2).reshape(101,-1);B=((V-V.mean(0))/sc).transpose(1,0,2).reshape(101,-1)
 u,s,v=svd(A,full_matrices=False,check_finite=False)
 return [(u[:,:q]*s[:q])@v[:q] if q else np.zeros_like(A) for q in qs],B
rows=[]
for fold in range(8):
 nums=np.zeros(len(qs));dens=nums.copy()
 for inner in range(8):
  if inner==fold:continue
  pred,B=prepare((groups!=fold)&(groups!=inner),groups==inner)
  for i,A in enumerate(pred):nums[i]+=np.sum(A*B);dens[i]+=np.sum(A*A)
 gains=np.divide(nums,dens,out=np.zeros_like(nums),where=dens>0);gains=np.maximum(gains,0)
 pred,B=prepare(groups!=fold,groups==fold)
 for i,(q,A) in enumerate(zip(qs,pred)):
  rows.append({'fold':fold,'rank':q,'gain':gains[i],'SSE':float(np.sum((B-gains[i]*A)**2)),'test_energy':float(np.sum(B*B))})
r=pd.DataFrame(rows);r.to_csv(P/'fastimages_temporal_rank_nested_folds.csv',index=False)
s=r.groupby('rank',as_index=False)[['SSE','test_energy']].sum();s['R2_to_zero_contrast']=1-s.SSE/s.test_energy;s.to_csv(P/'fastimages_temporal_rank_nested_summary.csv',index=False)
print(s.to_string(index=False))
(P/'temporal_rank_amendment.json').write_text(json.dumps({'reason':'Unshrunk training class-mean temporal kernels overfit noisy heldout block means at all tested ranks. Added nested nonnegative scalar gain calibration, selected only from the seven outer-training blocks, before comparing any outer-heldout means.','scope':'This controls a real bias-variance problem. Neither rank-one failure nor higher-rank advantage is presumed. All stated ranks reported, one person only.'},indent=2))
