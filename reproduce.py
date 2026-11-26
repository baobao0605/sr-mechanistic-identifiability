#!/usr/bin/env python3
"""Run release checks in a temporary copy, leaving source data unchanged."""
import argparse,hashlib,json,os,platform,shutil,subprocess,sys,tempfile
from datetime import datetime,timezone
from importlib.metadata import version
from pathlib import Path
ROOT=Path(__file__).resolve().parent
SHA='1c4cdfeb625a195530f8a08b685a03f7f654a8d033b9db37e3152308f87f60ff'
STEPS={
 'math_checks.py':['source_data/math_checks.json'],
 'math_additions.py':['source_data/math_additions.json','source_data/partial_discount_bounds.csv'],
 'empirical_check_cyclic_lag_coordinates.py':['empirical/finite_kernel_moment_scope.json'],
 'empirical_check_dfb_tables.py':['empirical/dFB_table_checks.json'],
 'check_fig07_tables.py':['empirical/fig07_table_recheck.json'],
}
OPTIONAL='empirical_check_fig07_portable.py'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--require-measured-input',action='store_true')
    p.add_argument('--measured-input',type=Path)
    p.add_argument('--out',type=Path,default=ROOT/'reproduction_runs/latest')
    a=p.parse_args()
    src=a.measured_input or ROOT/'source_data/fig07_empirical_arrays.npz'
    if (a.measured_input or a.require_measured_input) and not src.is_file():p.error('Measured input missing; see source_data/README_DATA.md.')
    if src.is_file() and hashlib.sha256(src.read_bytes()).hexdigest()!=SHA:p.error('Measured input SHA-256 does not match.')
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
    steps=dict(STEPS)
    if src.is_file():steps[OPTIONAL]=['empirical/fig07_portable_checks.json']
    report=[]
    with tempfile.TemporaryDirectory(prefix='sr-reproduce-') as t:
        work=Path(t)
        for folder in ('code','source_data','empirical'):
            shutil.copytree(ROOT/folder,work/folder,ignore=shutil.ignore_patterns('fig07_empirical_arrays.npz','__pycache__'))
        if src.is_file():shutil.copy2(src,work/'source_data/fig07_empirical_arrays.npz')
        for script, outputs in steps.items():
            # Remove archived output before running so stale records cannot appear as fresh results.
            for rel in outputs:(work/rel).unlink(missing_ok=True)
            job=subprocess.run([sys.executable,str(work/'code'/script)],cwd=work,capture_output=True,text=True,env={**os.environ,'PYTHONUTF8':'1'})
            (out/(Path(script).stem+'.log')).write_text(job.stdout+'\n'+job.stderr,encoding='utf-8')
            ok=job.returncode==0 and all((work/rel).is_file() for rel in outputs)
            status='PASS' if ok else 'FAIL'
            print(status,script,flush=True);report.append({'script':script,'status':status,'returncode':job.returncode})
            if ok:
                for rel in outputs:
                    dest=out/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(work/rel,dest)
        if not src.is_file():
            report.append({'script':OPTIONAL,'status':'SKIPPED_OPTIONAL_INPUT_MISSING'})
            print('SKIP optional measured-array refit (see data guide).')
    rec={'timestamp_utc':datetime.now(timezone.utc).isoformat(),'python':platform.python_version(),'packages':{k:version(k) for k in ('numpy','pandas','scipy','sympy','PyMuPDF')},'steps':report,'scope':'Symbolic checks, archived table aggregation and optional measured-array downstream refits. No raw-acquisition rerun or natural SR certification.'}
    (out/'report.json').write_text(json.dumps(rec,indent=2)+'\n')
    return int(any(x['status']=='FAIL' for x in report))
if __name__=='__main__':raise SystemExit(main())
