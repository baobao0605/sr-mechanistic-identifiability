"""Exact cyclic SR geometry compressed from actual THINGS references.

Only reference-training records set the whitening and spectral weights.
The compression uses Parseval, not a learned state reduction. No virtual
subjects/no artificial measurement noise. Eigenvalue 1 is removed by centering.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib
import numpy as np
from sklearn.covariance import LedoitWolf

ROOT=Path(__file__).resolve().parents[1]
OUT=Path(__file__).resolve().parent
SOURCE=ROOT/'prior/restored/complete/analysis/D05/analysis_input.npz'
STATES=(4,8,16,32,64,100)
CHANNELS=(6,12,24,48,72)
REPEATS=(1,2,3)
FOLDS=((0,1,2,3),(1,2,3,0),(2,3,0,1),(3,0,1,2))

def main():
 z=np.load(SOURCE);patterns=z['patterns'];roi=z['voxel_roi']
 groups=[np.arange(3*i,3*i+3) for i in range(4)]
 meta=[];weights=[];checks=[]
 for n in STATES:
  for nc in CHANNELS:
   ix=np.concatenate([np.flatnonzero(roi==r)[:nc//3] for r in ('V1','lFFA','rFFA')])
   data=patterns[:,:n,ix];data=data-data.mean(1,keepdims=True)
   for fold,(ga,gb,gc,gd) in enumerate(FOLDS):
    resid=data[groups[ga]]-data[groups[ga]].mean(0)
    cov=LedoitWolf().fit(resid.reshape(-1,nc)).covariance_
    ev,V=np.linalg.eigh(cov);W=(V/np.sqrt(np.maximum(ev,1e-20)))@V.T
    dw=data@W
    for nr in REPEATS:
     A=dw[groups[ga][:nr]].mean(0)
     F=np.fft.fft(A,axis=0,norm='ortho')
     w=np.sum(np.abs(F)**2,axis=1); w[0]=0
     padded=np.zeros(100);padded[:n]=w;weights.append(padded)
     meta.append({'case_id':len(meta),'states':n,'channels':nc,'repeats':nr,'fold':fold,'reference_norm':float(np.linalg.norm(A)),'metric_group':ga,'scope':'training-reference geometry; not clean latent B'})
     # Independent matrix identity checks on every case, with two unequal g.
     P=np.roll(np.eye(n),-1,axis=1)
     X=P@np.linalg.solve(np.eye(n)-.31*P,A)
     Z=P.T@np.linalg.solve(np.eye(n)-.78*P.T,A)
     xf=X-(np.sum(A*X)/np.sum(A*A))*A
     zf=Z-(np.sum(A*Z)/np.sum(A*A))*A
     lam=np.exp(-2j*np.pi*np.arange(n)/n)
     f=lam/(1-.31*lam);b=lam.conj()/(1-.78*lam.conj())
     f-=np.sum(w*f)/np.sum(w);b-=np.sum(w*b)/np.sum(w)
     numeric=float(np.sum(xf*zf));spectral=float(np.real(np.sum(w*f.conj()*b)))
     err=abs(numeric-spectral)/max(np.linalg.norm(xf)*np.linalg.norm(zf),1e-30)
     checks.append(err)
 assert max(checks)<1e-11
 np.savez_compressed(OUT/'spectral_reference_cases.npz',weights=np.array(weights))
 (OUT/'spectral_reference_cases.json').write_text(json.dumps(meta,indent=2))
 (OUT/'compression_checks.json').write_text(json.dumps({'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'cases':len(meta),'maximum_normalized_inner_product_difference':max(checks),'construction':'Parseval exact inner products of circulant functions; no estimated noise floor or latent B claim'},indent=2))
 print(json.dumps({'cases':len(meta),'maximum_normalized_difference':max(checks)}))

if __name__=='__main__': main()
