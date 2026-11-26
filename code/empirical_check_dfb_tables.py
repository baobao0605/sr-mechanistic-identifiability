"""Recompute published dFB count summaries from archived condition-level tables.

This is a deterministic source-table audit. It neither refits neural models nor
assigns probability coverage to the observed repeat residuals.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def truth(series):
    return series.astype(str).str.lower().eq('true')


def main():
    sd = ROOT / 'source_data'
    p = pd.read_csv(sd / 'dFB_profiled_decision_folds.csv')
    a = pd.read_csv(sd / 'dFB_observed_class_audit.csv')
    keys = ['dataset', 'window', 'fold', 'space', 'true_kind', 'true_gamma',
            'prediction_residual_strength', 'fold_id']
    merged = p.merge(a, on=keys, validate='one_to_one')
    d = merged[merged.true_kind.isin(['forward', 'backward'])].copy()
    assert len(p) == len(a) == 912
    assert len(d) == 864
    counts = []
    for label, suffix in [('full', 'full'), ('common_prediction', 'U'),
                          ('common_prediction_current_removed', 'Q_U')]:
        row = {'space': label, 'total': len(d)}
        for rule, out in [('binary', 'binary'), ('multi', 'multiclass')]:
            passed = truth(d[f'{rule}_oracle_{suffix}_bound_pass'])
            row[out] = int(passed.sum())
            assert truth(d.loc[passed, f'{rule}_correct']).all()
        counts.append(row)
    expected = [(0, 0), (6, 6), (80, 77)]
    assert [(r['binary'], r['multiclass']) for r in counts] == expected
    saved = pd.read_csv(sd / 'dFB_error_space_counts.csv').to_dict('records')
    assert saved == [{k: r[k] for k in ['space', 'binary', 'multiclass', 'total']}
                     for r in counts]
    gate = d[truth(d.multi_oracle_Q_U_bound_pass)]
    old = gate.old_multi_observed_winners.value_counts().to_dict()
    correct = int(truth(gate.multi_visible_unique_correct).sum())
    assert len(gate) == 77 and correct == 39
    assert old.get('current', 0) == 37 and old.get('symmetric', 0) == 1
    lower_margin = d.binary_margin_lower.to_numpy(float)
    error = d.oracle_error_Q_U.to_numpy(float)
    numeric_bound = lower_margin > 2 * error
    assert np.array_equal(numeric_bound, truth(d.binary_oracle_Q_U_bound_pass).to_numpy())
    lower_margin = d.multi_margin_lower.to_numpy(float)
    assert np.array_equal(lower_margin > 2 * error,
                          truth(d.multi_oracle_Q_U_bound_pass).to_numpy())
    result = {'passed': True, 'all_conditions': len(p),
              'direction_conditions': len(d), 'counts': counts,
              'multiclass_gate': {'n': len(gate), 'old_visible_correct': correct,
                                 'old_current': old.get('current', 0),
                                 'old_symmetric': old.get('symmetric', 0),
                                 'matching_rule_correct': int(truth(gate.multi_correct).sum())},
              'scope': 'Reaggregation of archived tables and verification of saved sufficient inequalities.'}
    out = ROOT / 'empirical/dFB_table_checks.json'
    out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
