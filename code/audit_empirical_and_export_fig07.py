"""Bounded source audit and source-data export for manuscript Figure 7.

This checks the supplied AB release. It does not regenerate upstream EEG
preprocessing, estimate population performance, or establish iid sampling.
No observations, participants or substitute source values are simulated.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
POST = ['P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4','PO8','O1','Oz','O2']

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main(ab):
    sd=ROOT/'source_data'; sd.mkdir(exist_ok=True)
    dest=ROOT/'empirical'; dest.mkdir(exist_ok=True)
    inputs=[]
    def csv(rel):
        p=ab/rel; inputs.append({'relative_path':rel,'sha256':sha(p)})
        return pd.read_csv(p)
    losses=csv('B_reference/source_data/all_control_reference_losses.csv')
    summary=csv('B_reference/source_data/all_control_reference_summary.csv')
    block=losses[(losses.dataset=='Turner')&(losses.context=='local')&(losses.estimator=='cv_fourier')].copy()
    block=block.sort_values(['session','block','window']).reset_index(drop=True)
    assert len(block)==8
    block['block_label']=['R%d B%d'%(r.session,r.block+1) for r in block.itertuples()]
    block['window_start_ms']=block.window.map({'early':100,'late':200})
    block['window_end_ms']=block.window.map({'early':150,'late':250})
    block['gain_percent']=100*block.plane_gain_vs_zero
    block['source']='B_reference/source_data/all_control_reference_losses.csv'
    cols=['key','session','block','block_label','window','window_start_ms','window_end_ms','control_count','state_count_min','state_count_max','family','lambda','plane_mse','plane_zero_mse','gain_percent','outer_mse','zero_mse','source']
    block[cols].to_csv(sd/'fig07_reference_blocks.csv',index=False)

    # Re-run the release's independently implemented refitter on raw supplied
    # derived arrays; this is a fresh execution of that independent audit,
    # not a claim to a third independent implementation.
    auditfile=ab/'reviews/final_repro_audit/audit_real_data.py'
    spec=importlib.util.spec_from_file_location('ab_independent_b_audit',auditfile)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    zpath=ab/'B_reference/source_data/all_control_and_deployment_B.npz'
    inputs.append({'relative_path':str(zpath.relative_to(ab)),'sha256':sha(zpath)})
    saved=np.load(zpath)
    checks=[]; splitrows=[]; raw={}; hier=[]; portable={}
    for session in [1,2]:
        p=ab/f'prior/sr_errorbounds_round4_20261124/inputs/round3/prior/restored/complete/analysis/D08/session{session}_localizer.npz'
        inputs.append({'relative_path':str(p.relative_to(ab)),'sha256':sha(p)})
        raw[session]=np.load(p)
        zz=raw[session]; ix=[list(zz['channels']).index(c) for c in POST]
        portable[f'R{session}__measurements']=np.stack([(zz['features'][:,wi]-zz['baseline'])[:,ix] for wi in [1,2]],axis=1)
        for name in ['labels','time','block']:
            portable[f'R{session}__{name}']=zz[name]
        for window,d in mod.load_turner(ab,session).items():
            for physical_block in sorted(np.unique(d['blocks'])):
                mask=d['blocks']==physical_block;t=d['times'];mid=(t[mask].min()+t[mask].max())/2
                train=mask&(t<=mid-.5);test=mask&(t>=mid+.5)
                pools=[d['x'][train&(d['labels']==j)] for j in range(40)]
                target=mod.center(np.array([d['x'][test&(d['labels']==j)].mean(0) for j in range(40)]))
                fits,counts=mod.independently_fit(pools)
                key=f'Turner_s{session}_b{physical_block}_{window}_all_local'
                overlap=len(set(d['event_ids'][train])&set(d['event_ids'][test]))
                gap=float(t[test].min()-t[train].max())
                assert overlap==0 and gap>=1
                splitrows.append(dict(key=key,session=session,block=int(physical_block),window=window,
                    controls=int(train.sum()),test_events=int(test.sum()),event_overlap=overlap,onset_gap_s=gap))
                for method,(estimate,params) in fits.items():
                    portable['reference__'+key+'__'+method]=saved[key+'__Bhat_'+method]
                    got=mod.losses(estimate,target)
                    reported=losses[(losses.key==key)&(losses.estimator==method)].iloc[0]
                    matrix_err=float(np.max(np.abs(estimate-saved[key+'__Bhat_'+method])))
                    metric_err=max(abs(v-float(reported[k])) for k,v in got.items())
                    checks.append(dict(key=key,estimator=method,max_matrix_error=matrix_err,max_metric_error=metric_err))
                    assert matrix_err<1e-10 and metric_err<1e-10
        z=raw[session]
        for physical_block in sorted(np.unique(z['block'])):
            hier.append(dict(participant='one supplied Turner participant (prior release: ID13)',session=session,
                physical_block=int(physical_block),block_label=f'R{session} B{int(physical_block)+1}',
                eligible_events=int((z['block']==physical_block).sum()),states=len(np.unique(z['labels'])),
                channels=len(POST),source='supplied session%d_localizer.npz'%session))
    pd.DataFrame(checks).to_csv(dest/'fig07_b_recomputation.csv',index=False)
    methods=['zero','mean','cv_scalar','cv_svd','cv_fourier']
    losses[(losses.dataset=='Turner')&(losses.context=='local')&(losses.estimator.isin(methods))].to_csv(sd/'fig07_all_reference_estimators.csv',index=False)
    pd.DataFrame(splitrows).to_csv(sd/'fig07_temporal_splits.csv',index=False)
    pd.DataFrame(hier).to_csv(sd/'fig07_acquisition_hierarchy.csv',index=False)
    agg=[]
    for window,d in block.groupby('window'):
        gain=1-d.plane_mse.mean()/d.plane_zero_mse.mean()
        reported=summary[(summary.dataset=='Turner')&(summary.context=='local')&(summary.estimator=='cv_fourier')&(summary.window==window)].iloc[0]
        # Source summary column naming is preserved in audit input; deriving
        # the aggregate from its four constituent loss rows avoids ambiguity.
        agg.append(dict(window=window,physical_blocks=len(d),mean_plane_mse=d.plane_mse.mean(),
            mean_zero_plane_mse=d.plane_zero_mse.mean(),gain_percent=100*gain,
            positive_blocks=int((d.plane_gain_vs_zero>0).sum()),
            source='ratio of mean MSEs in fig07_reference_blocks.csv'))
    pd.DataFrame(agg).to_csv(sd/'fig07_reference_aggregate.csv',index=False)

    positives=csv('theory_joint/J01_all_positive_audit.csv')
    p=positives[(positives.variant=='B_trained_A')&(positives.truth_family=='canonical_SR')].copy()
    assert len(p)==3 and (p.true_kind=='backward').all()
    p['design']=np.where(p.context.str.contains('fixed'),'Native order, k = 1','Designed order, k = 10')
    p['source']='theory_joint/J01_all_positive_audit.csv'
    p.to_csv(sd/'fig07_conditional_sr.csv',index=False)

    readpath=ab/'joint_integration/source_data/joint_readouts.npz'
    inputs.append({'relative_path':str(readpath.relative_to(ab)),'sha256':sha(readpath)})
    read=np.load(readpath); reconstruction=[]
    z=raw[2];ids=np.arange(len(z['labels']));counts=np.bincount(z['labels'],minlength=40)
    nr=int(counts.min());index=np.array([np.flatnonzero(z['labels']==s)[:nr] for s in range(40)]).T
    channels=[list(z['channels']).index(c) for c in POST]
    x=(z['features'][:,1]-z['baseline'])[:,channels][index]
    spans=z['time'][index].max(1)-z['time'][index].min(1)
    stats=dict(participants=1,records_in_positive_context=1,physical_blocks_in_positive_record=len(np.unique(z['block'])),
               complete_state_maps=nr,training_state_maps=4,training_events=4*x.shape[1],
               inference_events=(nr-4)*x.shape[1],heldout_regrouped_pairs=(nr-4)//2,
               maps_spanning_both_blocks=int(np.sum(np.ptp(z['block'][index],axis=1)>0)),
               median_state_map_span_s=float(np.median(spans)),states=x.shape[1],channels=x.shape[2])
    for context in p.context.unique():
        pref=context+'__B_trained_A__'
        C=read[pref+'Ceff'];W=read[pref+'W'];H=read[pref+'H'];R=read[pref+'R']
        for name in ['Ceff','W','H','R','all_readouts','training_event_ids']:
            portable[pref+name]=read[pref+name]
        got=np.einsum('aij,rij->ra',H,np.einsum('ij,rjk,kl->ril',C,x,W))/np.linalg.norm(R)
        err=float(np.max(abs(got-read[pref+'all_readouts'])))
        eids=index[:4].ravel();reported=read[pref+'training_event_ids'].ravel()
        assert np.array_equal(eids,reported)
        assert not (set(eids)&set(index[4:].ravel()))
        assert err<1e-10
        reconstruction.append(dict(context=context,variant='B_trained_A',maximum_readout_error=err,
            training_event_identity=True,training_inference_event_overlap=0))
    pd.DataFrame(reconstruction).to_csv(dest/'fig07_j01_readout_recomputation.csv',index=False)
    pd.DataFrame([dict(metric=k,value=v,source='reconstructed from supplied session2_localizer.npz') for k,v in stats.items()]).to_csv(sd/'fig07_regrouping_hierarchy.csv',index=False)

    phys=csv('A_measurement/integration_review/physical_block_decomposition.csv')
    selected=phys[(phys.context=='Turner_early_s2_designed_quarter_k10')&(phys.gamma==.95)&(phys.truth=='backward')&(phys.correction=='joint_family')].copy()
    assert len(selected)==2 and (selected.decision=='abstain').all()
    selected.to_csv(sd/'fig07_within_physical_block.csv',index=False)
    tr=csv('joint_integration/source_data/transport_stress_results.csv')
    selectedtr=tr[(tr.transport=='cross_record')&(tr.correction=='joint_family_nominal')].copy()
    assert len(selectedtr)==1200 and (selectedtr.decision=='abstain').all()
    trsummary=selectedtr.groupby('dataset').agg(computational_conditions=('decision','size'),
        declared=('declared','sum'),abstentions=('decision',lambda x:(x=='abstain').sum())).reset_index()
    trsummary['source']='joint_integration/source_data/transport_stress_results.csv'
    trsummary.to_csv(sd/'fig07_transport_summary.csv',index=False)
    selectedtr[['context','variant','dataset','window','session','transport','n_pairs','truth_family','true_kind','gamma','rho','correction','decision']].to_csv(sd/'fig07_transport_conditions.csv',index=False)
    verdict=pd.DataFrame([
        dict(endpoint='Within each physical block',denominator=len(selected),declarations=int((selected.decision!='abstain').sum()),abstentions=int((selected.decision=='abstain').sum()),unit='two block-specific reanalyses of one record; not participants',source='fig07_within_physical_block.csv'),
        dict(endpoint='Cross-record stress',denominator=len(selectedtr),declarations=int(selectedtr.declared.sum()),abstentions=int((selectedtr.decision=='abstain').sum()),unit='computed conditions reusing recordings; not participants',source='fig07_transport_conditions.csv')])
    verdict.to_csv(sd/'fig07_decision_summary.csv',index=False)
    # The optional fMRI example is retained in the supplement, not Figure 7.
    losses[losses.dataset=='D02'].to_csv(sd/'empirical_ekman_reference_all.csv',index=False)
    np.savez_compressed(sd/'fig07_empirical_arrays.npz',**portable)
    result=dict(date='2026-11-25',input_root=str(ab),input_hashes=inputs,
        current_recomputation=dict(reference_contexts=8,reference_estimator_checks=len(checks),
            maximum_matrix_error=max(r['max_matrix_error'] for r in checks),
            maximum_metric_error=max(r['max_metric_error'] for r in checks),
            joint_contexts=len(reconstruction),maximum_joint_readout_error=max(r['maximum_readout_error'] for r in reconstruction),
            event_overlap=0,raw_eeg_preprocessing_rerun=False),
        acquisition_hierarchy=stats,transport_conditions=len(selectedtr),
        transport_all_abstain=bool((selectedtr.decision=='abstain').all()),
        qualifications=['Previously derived measured EEG arrays, not raw EEG refit.',
          'The B refit executes the prior independent implementation; not a third independent algorithm.',
          'No iid, Gaussianity, population performance or endogenous SR claim follows from numerical reproduction.',
          'Early 100–150 ms and late 200–250 ms window bounds come from the prior statistical appendix; arrays retain named windows only.',
          'Participant ID13 is inherited from the prior A report; no participant identifier is embedded in NPZ arrays.'])
    (dest/'empirical_audit_results.json').write_text(json.dumps(result,indent=2,ensure_ascii=False))
    print(json.dumps(result['current_recomputation'],indent=2))
    print(pd.DataFrame(agg).to_string(index=False))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--ab',type=Path,required=True);args=ap.parse_args();main(args.ab)
