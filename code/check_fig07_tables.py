"""Audit archived Figure 7 endpoints without participant-level input arrays."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
sd = ROOT / 'source_data'
b = pd.read_csv(sd/'fig07_reference_blocks.csv')
g = pd.read_csv(sd/'fig07_reference_aggregate.csv')
assert len(b) == 8
values = {}
for w, rows in b.groupby('window'):
    value = 100*(1-rows.plane_mse.mean()/rows.plane_zero_mse.mean())
    assert abs(value-g[g.window==w].gain_percent.iloc[0]) < 1e-10
    values[w] = float(value)
p = pd.read_csv(sd/'fig07_conditional_sr.csv')
assert len(p) == 3 and (p.ratio > 1).all()
assert np.max(abs(p.boundary_rational/p.critical-p.ratio)) < 1e-12
tr = pd.read_csv(sd/'fig07_transport_conditions.csv')
summary = pd.read_csv(sd/'fig07_transport_summary.csv')
assert len(tr)==1200 and (tr.decision=='abstain').all()
for dataset, rows in tr.groupby('dataset'):
    row = summary[summary.dataset==dataset].iloc[0]
    assert len(rows)==row.computational_conditions==row.abstentions
blocks = pd.read_csv(sd/'fig07_within_physical_block.csv')
assert len(blocks)==2 and (blocks.decision=='abstain').all()
result = {'passed':True,'gain_percent':values,'conditional_positive_rows':len(p),'physical_block_abstentions':len(blocks),'cross_record_abstentions':len(tr),'scope':'Archived endpoint aggregation, not measured-array refit or raw preprocessing.'}
(ROOT/'empirical/fig07_table_recheck.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
