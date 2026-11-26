from pathlib import Path
import h5py,numpy as np,json,hashlib
P=Path(__file__).resolve().parent
ROOT=P.parent
C=P/'inputs/natural_scene_counts.npz'
z=np.load(C);sel=z['areas']=='VISp'; ids=z['unit_ids'][sel]
cache=P/'allen_visp_25ms_counts.npz'
if not cache.exists():
 with h5py.File(ROOT/'raw/D14/sub-707296975_ses-721123822.nwb','r') as f:
  unitrows={a:i for i,a in enumerate(f['units/id'][:])};ends=f['units/spike_times_index'][:];beg=np.r_[0,ends[:-1]]; edges=np.linspace(0,.25,11)
  counts=np.zeros((len(z['state']),10,len(ids)),np.uint16)
  for j,uid in enumerate(ids):
   u=unitrows[uid];ts=f['units/spike_times'][beg[u]:ends[u]];ix=np.searchsorted(ts,z['start'][:,None]+edges,side='left');counts[:,:,j]=np.diff(ix,axis=1)
  np.savez_compressed(cache,counts=counts,unit_ids=ids,edges=edges)
# Same guarded outer groups and nested regularization as the final event-history model,
# but fit each actual 25ms window's predictive coefficients separately.
s=(P/'allen_event_history.py').read_text()
s=s.replace("YY=z['response'].astype(float); rois={'VISp':z['areas']=='VISp','all_visual':np.ones(len(z['areas']),bool)}", "bins=np.load(OUT/'allen_visp_25ms_counts.npz'); YY=bins['counts'].astype(float).reshape(len(y),-1); rois={f'bin_{j*25:03d}_{(j+1)*25:03d}ms':np.repeat(np.arange(10)==j,110) for j in range(10)}")
s=s.replace("'allen_history_", "'allen_timebins_")
s=s.replace("['SSE_standardized','baseline_SSE_standardized','n_test']", "['SSE_standardized','baseline_SSE_standardized','n_test']")
s=s.replace("'lag_units':'one actual event, median spacing seconds '","'timebin_scope':'VISp 110 units only; ten nonoverlapping 25ms windows from 0 to250ms; actual next onset always >250ms; hyperparameters tuned separately per bin using training folds.','lag_units':'one actual event, median spacing seconds '")
(P/'allen_timebins_model.py').write_text(s)
exec(compile(s,str(P/'allen_timebins_model.py'),'exec'))
