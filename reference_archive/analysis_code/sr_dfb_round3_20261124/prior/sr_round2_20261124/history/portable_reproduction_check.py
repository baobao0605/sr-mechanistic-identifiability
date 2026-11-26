"""Meaningful portability gate: rerun three primary analyses from an unrelated tree."""
from pathlib import Path
import tempfile,shutil,subprocess,sys,json,os
import numpy as np,pandas as pd
P=Path(__file__).resolve().parent
T=Path(tempfile.mkdtemp(prefix='sr_history_portable_'))
try:
 for name in ['allen_event_history.py','allen_timebins.py','allen_visp_25ms_counts.npz','fastimages_kernel_competition.py','fastimages_kernel_input.npz','fastimages_extraction_audit.json']:
  shutil.copy2(P/name,T/name)
 shutil.copytree(P/'inputs',T/'inputs')
 # No raw directory is copied. Numerical source caches are sufficient.
 env=dict(os.environ,OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
 runs=[]
 for name in ['allen_event_history.py','allen_timebins.py','fastimages_kernel_competition.py']:
  cp=subprocess.run([sys.executable,str(T/name)],cwd=T,env=env,text=True,capture_output=True)
  if cp.returncode: raise RuntimeError(name+'\n'+cp.stdout+'\n'+cp.stderr)
  runs.append({'script':name,'exit_code':cp.returncode});print('portable run complete',name,flush=True)
 checks=[]
 for filename in ['allen_history_summary.csv','allen_timebins_summary.csv','fastimages_kernel_summary.csv']:
  a=pd.read_csv(P/filename);b=pd.read_csv(T/filename);assert a.shape==b.shape
  nonnum=list(a.select_dtypes(exclude='number').columns);a=a.sort_values(nonnum).reset_index(drop=True);b=b.sort_values(nonnum).reset_index(drop=True)
  for c in nonnum:assert a[c].equals(b[c])
  c=list(a.select_dtypes('number').columns);err=float(np.max(np.abs(a[c].to_numpy()-b[c].to_numpy())))
  assert np.allclose(a[c],b[c],rtol=1e-12,atol=1e-10)
  checks.append({'file':filename,'max_abs_difference':err,'identical_row_labels':True})
 result={'status':'passed','no_raw_files_copied':True,'absolute_original_data_paths_required':False,'runs':runs,'numerical_checks':checks}
 (P/'portable_reproduction_audit.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(result,indent=2),flush=True)
finally:shutil.rmtree(T)
