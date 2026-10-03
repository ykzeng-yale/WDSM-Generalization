#!/usr/bin/env python3
"""Complete-grid reporting only; no RDS decoding, fitting, simulation or RNG."""
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import formal_metrics


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def need(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def endpoints(x):
    if isinstance(x, dict):
        return [x.get('lower'), x.get('upper')]
    return x


def main(argv):
    need(len(argv) == 3, 'Usage: summarize.py TRUTH_DIR CASE_ROOT NEW_SUMMARY_DIR')
    truth_dir, case_root, output = map(Path, argv)
    need(not output.exists() and output.parent.is_dir(), 'Summary output must be fresh')
    source = Path(__file__).resolve().parent
    manifest_hash = sha(source / 'source_inventory.json')
    inventory = json.loads((source / 'source_inventory.json').read_text())
    package_root = source.parents[1]
    source_pins = [(package_root/x['path'],x) for x in inventory['files']]
    for path,pin in source_pins:
        need(sha(path) == pin['sha256'] and path.stat().st_size == pin['bytes'], 'Reporting/common source mismatch')
    truth_receipt_path = truth_dir / 'result.json'
    truth_receipt_hash = sha(truth_receipt_path)
    tr = json.loads(truth_receipt_path.read_text())
    need(tr['schema'] == 'principal_iid_public_truth_v1' and tr['status'] == 'COMPLETE_TRUTH' and
         tr['source_input_end_verified'] is True, 'Complete truth output required')
    need(tr['source_inventory_sha256'] == manifest_hash, 'Truth/source mismatch')
    truth_names = {'model.rds', 'mixture_means.csv', 'truth_orders.csv',
                   'truth_components.csv', 'truth_targets.csv', 'truth.rds'}
    truth_files = tr['files']
    need(isinstance(truth_files,list) and len(truth_files) == 6 and
         {x['path'] for x in truth_files} == truth_names and
         truth_files == tr['artifact_manifest'], 'Exact six-file truth inventory required')
    selected_truth_inputs = {'result.json': truth_receipt_hash}
    truth_pins = [dict(path=str(truth_receipt_path),sha256=truth_receipt_hash,
                       bytes=truth_receipt_path.stat().st_size)]
    for item in truth_files:
        path = truth_dir / item['path']
        need(sha(path) == item['sha256'] and path.stat().st_size == item['bytes'],
             'Selected truth/model artifact hash/size mismatch')
        selected_truth_inputs[item['path']] = item['sha256']
        truth_pins.append(dict(path=str(path),sha256=item['sha256'],bytes=item['bytes']))
    tt = truth_dir / 'truth_targets.csv'
    truths = {}
    with tt.open(newline='') as f:
        for row in csv.DictReader(f):
            key = (row['overlap'], row['estimand'])
            need(key not in truths and row['admitted'] == 'TRUE', 'Nonadmitted/duplicate target')
            truths[key] = {k: float(row[k]) for k in ('target', 'sensitivity_lower', 'sensitivity_upper')}
            truths[key]['truth_method'] = row['truth_method']
    need(set(truths) == {(o,e) for o in ('GoodOverlap','PoorOverlap') for e in formal_metrics.ESTIMANDS}, 'Target universe mismatch')
    cases = sorted(case_root.glob('*/result.json'))
    need(len(cases) == 4000, 'Exactly 4000 completed case results required; no partial metrics')
    seen, records, inputs = set(), [], list(truth_pins)
    for path in cases:
        result = json.loads(path.read_text())
        need(result['schema'] == 'principal_iid_public_case_v1' and result['status'] in
             ('COMPLETE_CASE','COMPLETE_WITH_STATISTICAL_FAILURES'), f'Incomplete/integrity case: {path}')
        need(result['source_input_end_verified'] is True and result['post_count_RNG_unchanged'] is True,
             f'End-state mismatch: {path}')
        need(result['source_inventory_sha256'] == manifest_hash, 'Case source differs')
        case_inputs = result['inputs']
        need(isinstance(case_inputs,list) and len(case_inputs) == 7 and
             all(set(x) == {'path','sha256'} for x in case_inputs) and
             len({x['path'] for x in case_inputs}) == 7 and
             {x['path']:x['sha256'] for x in case_inputs} == selected_truth_inputs,
             'Case did not use the selected complete model/truth artifacts and receipt')
        g, n, rep, overlap = result['group_index'], result['n'], result['replicate'], result['overlap']
        need((g,overlap,n) in formal_metrics.GROUPS and isinstance(rep,int) and 1 <= rep <= 1000, 'Invalid case key')
        case = (g,rep)
        need(case not in seen and result['case_id'] == f'g{g}_{overlap}_n{n}_rep{rep:04d}', 'Duplicate/key mismatch')
        seen.add(case)
        need(result['data_seed'] == 104000000+g*100000+rep and result['count_seed'] == 204000000+g*100000+rep,
             'Predetermined seed mismatch')
        need(result['M'] == 3 and result['B'] == 200, 'M/B mismatch')
        rows = result['rows']
        expected = {f'{f}_{e}' for e in formal_metrics.ESTIMANDS for f in formal_metrics.FAMILIES}
        need(set(rows) == expected, 'Missing method row')
        artifact_map = {x['path']:x for x in result['artifact_manifest']}
        need(len(artifact_map) == len(result['artifact_manifest']), 'Duplicate artifact declaration')
        pairing = set()
        for key, row in rows.items():
            detail = path.parent / (key+'.json')
            pin = artifact_map[key+'.json']
            need(sha(detail) == pin['sha256'] and detail.stat().st_size == pin['bytes'], 'Method JSON hash/size mismatch')
            need(json.loads(detail.read_text()) == row, 'Method/result row mismatch')
            inputs.append(dict(path=str(detail),sha256=pin['sha256']))
            f,e = row['family'],row['estimand']
            expected_p = {'PS':(26,18),'DSM':(42,26),'X6':(68,40)}[f][e=='PATT']
            need(key == f'{f}_{e}' and row['p'] == expected_p and row['case_id'] == result['case_id'] and
                 all(row[k] == result[k] for k in ('group_index','overlap','n','replicate','M','B')), 'Full-p/method binding mismatch')
            need(row['target'] == truths[(overlap,e)]['target'], 'Case target mismatch')
            pairing.add((row['sample_sha256'], row['counts_sha256']))
            need(all(isinstance(row[k],bool) for k in ('analytic_available','refit_available','contribution_ok')), 'Invalid availability status')
            point_ok = row['point_status'] == 'ok'
            need(row['point_status'] in ('ok','failed'), 'Nonterminal point status')
            if point_ok:
                need(finite(row['estimate']) and finite(row['raw_estimate']), 'Invalid point success')
            r = dict(row, point_available=point_ok, contribution_available=row['contribution_ok'])
            need(not row['contribution_ok'] or row['analytic_available'],
                 'Contribution availability inconsistent')
            failed = {}
            for field in ('count_fit_failures','count_draw_failures'):
                # jsonlite auto_unbox preserves a singleton integer as a scalar.
                value = row[field]
                ix = [value] if type(value) is int else value
                need(isinstance(ix,list) and all(type(i) is int and 1 <= i <= 200 for i in ix)
                     and ix == sorted(set(ix)), 'Lost/invalid count failure positions')
                failed[field] = ix
                r[field] = ix
            need(not row['refit_available'] or
                 not (failed['count_fit_failures'] or failed['count_draw_failures']),
                 'All-B refit success conceals failed slots')
            need(set(failed['count_fit_failures']).issubset(failed['count_draw_failures']),
                 'Failed count fit lacks its scalar failure')
            for kind in ('analytic','refit'):
                ci = endpoints(row[kind+'_interval'])
                need(isinstance(ci,list) and len(ci)==2, 'Invalid interval shape')
                r[kind+'_lower'], r[kind+'_upper'] = ci
                if row[kind+'_available']:
                    need(point_ok and finite(row[kind+'_variance']) and row[kind+'_variance'] >= 0 and
                         all(finite(v) for v in ci) and ci[0] <= ci[1], 'Invalid interval success')
            records.append(r)
        need(len(pairing) == 1, 'Methods do not retain one shared sample/count matrix')
        all_success = all(x['analytic_available'] and x['contribution_ok'] and x['refit_available']
                          for x in rows.values())
        need(result['status'] == ('COMPLETE_CASE' if all_success else
             'COMPLETE_WITH_STATISTICAL_FAILURES'), 'Case availability and status disagree')
        inputs.append(dict(path=str(path),sha256=sha(path)))
    need(seen == {(g,r) for g,_,_ in formal_metrics.GROUPS for r in range(1,1001)}, 'Incomplete requested universe')
    reports = formal_metrics.summarize(records, truths)
    need(sha(source/'source_inventory.json') == manifest_hash, 'Source inventory changed during reporting')
    need(all(sha(p) == pin['sha256'] for p,pin in source_pins) and
         all(sha(Path(x['path'])) == x['sha256'] and ('bytes' not in x or Path(x['path']).stat().st_size == x['bytes']) for x in inputs), 'Source/input changed during reporting')
    output.mkdir()
    for name, rows in reports.items():
        fields = list(dict.fromkeys(k for row in rows for k in row))
        with (output/(name+'.csv')).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields); writer.writeheader()
            for row in rows:
                writer.writerow({k:json.dumps(v,sort_keys=True,allow_nan=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})
    (output/'summary.json').write_text(json.dumps(reports,indent=2,allow_nan=False)+'\n')
    receipt=dict(status='COMPLETE_SUMMARY',cases=4000,point_rows=24000,inference_cells=48,
        source_inventory_sha256=manifest_hash,inputs=inputs,
        direct_checks='complete keys, selected six-file truth/model artifact and receipt hashes, exact case inputs, method JSON bytes, declared source/input/RNG end state, all-slot status/point pairing',
        array_scope='Case RDS arrays are not reread; their creation-time hashes remain in case results. Selected truth/model artifacts are hashed without decoding',
        empirical_scope='iid synthetic repeated-sampling comparison; not a proof of population assumptions or finite-sample coverage')
    (output/'receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print('COMPLETE_SUMMARY')


if __name__ == '__main__':
    main(sys.argv[1:])
