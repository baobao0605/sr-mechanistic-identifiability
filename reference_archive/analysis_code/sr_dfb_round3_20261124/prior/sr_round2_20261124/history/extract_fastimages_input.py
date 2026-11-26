"""Recreate numeric MEG input cache from author's continuous FIF and audited events.
Usage: python extract_fastimages_input.py /absolute/path/to/raw.fif
"""
from pathlib import Path
import argparse,json,hashlib,re
import numpy as np,pandas as pd,mne
p=argparse.ArgumentParser();p.add_argument('fif',type=Path);p.add_argument('--output',type=Path,default=Path(__file__).resolve().parent);args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
e=pd.read_csv(Path(__file__).resolve().parent/'inputs/all_actual_visual_events.csv')
raw=mne.io.read_raw_fif(args.fif,preload=False,verbose=False);sf=raw.info['sfreq'];pick=mne.pick_types(raw.info,meg='mag',exclude=[]);assert sf==100 and len(pick)==102
pat=re.compile(r'^(localizer|sequence) (Face|House|Cat|Shoe|Chair)( distractor)? onset$')
ons=np.array([t-raw.first_time for t,d in zip(raw.annotations.onset,raw.annotations.description) if pat.match(d)])
assert len(ons)==len(e);dt=float(np.max(np.abs(ons-e.onset.to_numpy())));assert dt<1/sf
D=raw.get_data(picks=pick).T;loc=e[(e.stage=='localizer')&~e.distractor];E=[]
for row in loc.itertuples():
 i=int(round(row.onset*sf));E.append(D[i:i+101]-D[i-20:i].mean(0))
ys=[];sample_ids=[];onsets=[];labels=[];blocks=[];soas=[]
for sid,dd in e[e.stage=='sequence'].groupby('sequence_id',sort=True):
 dd=dd.sort_values('sequence_position');o=np.rint(dd.onset.to_numpy()*sf).astype(int);ix=np.arange(o[0],o[-1]+101)
 ys.append(D[ix]-D[o[0]-20:o[0]].mean(0));sample_ids.append(ix);onsets.append(o);labels.append(dd.label.to_numpy());blocks.append(int(dd.block.iloc[0]));soas.append(float(np.median(np.diff(dd.onset))*1000))
np.savez_compressed(out/'fastimages_kernel_input.npz',localizer_epochs=np.array(E),localizer_labels=loc.label.to_numpy(),localizer_blocks=loc.block.to_numpy(),trial_y_flat=np.concatenate(ys),trial_sample_ids_flat=np.concatenate(sample_ids),trial_indptr=np.r_[0,np.cumsum([len(x) for x in ys])],trial_onset_samples=np.array(onsets),trial_labels=np.array(labels),trial_blocks=np.array(blocks),trial_SOA_ms=np.array(soas))
meta={'n_people':1,'n_blocks':8,'n_sequences':len(ys),'n_ordinary_localizer':len(loc),'n_magnetometers':len(pick),'sfreq':sf,'source_first_time_seconds':raw.first_time,'source_highpass_Hz':raw.info['highpass'],'source_lowpass_Hz':raw.info['lowpass'],'source_sha256':hashlib.sha256(args.fif.read_bytes()).hexdigest(),'input_format':'Numeric arrays only; variable sequence lengths stored as flat arrays plus indptr.','max_event_annotation_difference_seconds':dt}
(out/'fastimages_extraction_audit.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta,indent=2))
