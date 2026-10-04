"""One frozen HC129 complete-record preparation; no model or effect fits."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import re

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-selected', type=Path, required=True,
                    help='source_selected.csv.gz produced by decode_source.py')
parser.add_argument('--dictionary', type=Path, required=True,
                    help='source_selected_dictionary.json produced by decode_source.py')
parser.add_argument('--output-dir', type=Path, required=True,
                    help='New private directory for complete records, ledgers and receipt')
args = parser.parse_args()
if not __debug__:
    parser.error('Run without -O: source identity and format assertions are required.')
SPEC = Path(__file__).resolve().with_name('FROZEN_SPEC.json')
SPEC_SHA = '49126a4f29a7883eb610677e89f72288ea8e7256f23149c948ca49125f98225f'
SOURCE = args.source_selected.resolve()
DICTIONARY = args.dictionary.resolve()
OUT = args.output_dir.resolve()
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(SPEC) == SPEC_SHA
spec = json.loads(SPEC.read_text())
# CSV content identity is independent of gzip headers and compression version.
with gzip.open(SOURCE, 'rb') as selected_csv:
    decoded_hash = hashlib.sha256()
    for chunk in iter(lambda: selected_csv.read(1024 * 1024), b''):
        decoded_hash.update(chunk)
source_decoded_sha256 = decoded_hash.hexdigest()
assert source_decoded_sha256 == spec['source']['decoded_sha256']
assert sha(DICTIONARY) == spec['dictionary']['sha256']
OUT.mkdir(parents=True, exist_ok=False)
started = datetime.now(timezone.utc).isoformat()
(OUT / 'FROZEN_SPEC.json').write_bytes(SPEC.read_bytes())
fields = spec['fields_order']
dictionary = json.loads(DICTIONARY.read_text())
labels = {}
missing = {}
for f in fields:
    labels[f] = {int(code): label for code, label in re.findall(
        r"^\s*(-?\d+)\s*=\s*'([^']*)'", dictionary[f]['official_format'], re.M)}
    missing[f] = {c: label for c, label in labels[f].items()
                  if c < 0 or (f == 'MARRY42X' and c == 6)}

with gzip.open(SOURCE, 'rt', newline='') as handle:
    rows = list(csv.DictReader(handle))
assert len({r['DUPERSID'] for r in rows}) == len(rows) == 36855

def val(row, f):
    x = float(row[f])
    if not math.isfinite(x):
        raise ValueError(f'Nonfinite source value {row["source_row"]}/{f}')
    return x

def group(row):
    h, r = val(row, 'HISPANX'), val(row, 'RACEX')
    if h == 1:
        return 'Hispanic'
    if h == 2:
        return {1: 'White', 2: 'Black', 4: 'Asian'}.get(r, 'Other_non_Hispanic')
    return 'Unclassified'

def invalid_reason(row, f):
    x = val(row, f)
    if x in missing[f]:
        return missing[f][x]
    if spec['fields'][f]['type'] == 'categorical':
        if x not in spec['fields'][f]['valid_levels']:
            raise ValueError(f'Unknown categorical code {row["source_row"]}/{f}={x}')
    else:
        rule = spec['continuous_validity'][f]
        if not rule['min'] <= x <= rule['max'] or (rule.get('integer') and x != int(x)):
            raise ValueError(f'Out-of-codebook value {row["source_row"]}/{f}={x}')
    return None

# Validate every source value, before cohort filtering. Missing codes are retained.
invalid = {r['DUPERSID']: {f: reason for f in fields
           if (reason := invalid_reason(r, f)) is not None} for r in rows}
assert all(group(r) == r['race_group'] for r in rows)

def write_csv(name, data, columns=None, compressed=False):
    path = OUT / name
    if columns is None:
        columns = list(data[0]) if data else []
    if compressed:
        raw = path.open('wb')
        stream = io.TextIOWrapper(gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0),
                                  encoding='utf-8', newline='')
    else:
        stream = path.open('w', newline='')
    with stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(data)
    if compressed:
        raw.close()

flow = []
reasons = []
field_summary = []
baseline = []
cohort_records = {}
schemas = {}
base_ledgers = []
meta = ['source_row', 'DUPERSID', 'PANEL', 'race_group', 'Z', 'Y', 'W',
        'TOTEXP09', 'SAQWT09F', 'PERWT09F', 'HISPANX', 'RACEX', 'RACETHNX',
        'SAQELIG', 'SFFLAG42', 'VARSTR', 'VARPSU']
for comparator in spec['pair_order']:
    pair = 'White-' + comparator
    eligible = rows
    def add_flow(step, before, after):
        flow.append({'pair': pair, 'step': step, 'n_before': len(before),
                     'n_excluded_here': len(before)-len(after), 'n_after': len(after),
                     'White_after': sum(group(r)=='White' for r in after),
                     'comparator_after': sum(group(r)==comparator for r in after)})
    filters = [('pair_race', lambda r: group(r) in ['White', comparator]),
               ('adult_AGE42X', lambda r: val(r, 'AGE42X') >= 18),
               ('positive_SAQWT09F', lambda r: val(r, 'SAQWT09F') > 0),
               ('valid_TOTEXP09', lambda r: val(r, 'TOTEXP09') >= 0)]
    for step, predicate in filters:
        after = [r for r in eligible if predicate(r)]
        add_flow(step, eligible, after)
        eligible = after
    base = eligible
    for row in rows:
        if group(row) not in ['White', comparator]:
            continue
        failed_base = [name for name, predicate in filters[1:] if not predicate(row)]
        base_ledgers.append({'pair':pair, 'source_row':row['source_row'], 'DUPERSID':row['DUPERSID'],
            'race_group':group(row), 'AGE42X':row['AGE42X'], 'SAQWT09F':row['SAQWT09F'],
            'PERWT09F':row['PERWT09F'], 'TOTEXP09':row['TOTEXP09'],
            'base_eligible':int(not failed_base), 'failed_base_rules':'|'.join(failed_base),
            'failed_covariates':'|'.join(invalid[row['DUPERSID']]),
            'retained':int(not failed_base and not invalid[row['DUPERSID']])})
    for f in fields:
        failed = [r for r in base if f in invalid[r['DUPERSID']]]
        counts = Counter((group(r), r[f], invalid[r['DUPERSID']][f]) for r in failed)
        for (race, code, reason), count in sorted(counts.items()):
            field_summary.append({'pair':pair, 'field':f, 'race_group':race,
                'raw_code':code, 'source_reason':reason, 'n_excluded_any_order':count})
        for row in failed:
            reasons.append({'pair':pair, 'source_row':row['source_row'], 'DUPERSID':row['DUPERSID'],
                'race_group':group(row), 'W':row['SAQWT09F'], 'SAQWT09F':row['SAQWT09F'],
                'PERWT09F':row['PERWT09F'], 'field':f, 'raw_code':row[f],
                'source_reason':invalid[row['DUPERSID']][f]})
        after = [r for r in eligible if f not in invalid[r['DUPERSID']]]
        add_flow('complete_record:' + f, eligible, after)
        eligible = after
    complete = eligible
    assert all(not invalid[r['DUPERSID']] for r in complete)
    assert [int(r['source_row']) for r in complete] == sorted(int(r['source_row']) for r in complete)
    cohort_records[pair] = complete
    # Same source-defined predictors for all methods; dummy code categorical variables.
    schema_fields = []
    columns = []
    for f in fields:
        values = sorted({val(r, f) for r in complete})
        if spec['fields'][f]['type'] == 'continuous':
            active = len(values) > 1
            names = [f] if active else []
            schema_fields.append({'source_field':f, 'type':'continuous', 'columns':names,
                                 'constant_value':values[0] if not active else None,
                                 'transformation':'identity; no scaling or fitted transformation'})
        else:
            levels = [int(x) for x in values]
            reference = levels[0]
            names = [f + '__' + str(x) for x in levels if x != reference]
            schema_fields.append({'source_field':f, 'type':'categorical',
                'declared_source_levels':spec['fields'][f]['valid_levels'],
                'observed_levels':levels, 'reference_level':reference,
                'zero_count_source_levels':[x for x in spec['fields'][f]['valid_levels'] if x not in levels],
                'columns':names, 'source_labels':labels[f]})
        columns.extend(names)
    schemas[pair] = {'target':spec['target'], 'n':len(complete), 'fields':schema_fields,
        'columns':columns, 'n_columns_excluding_intercept':len(columns),
        'intercept':'Not stored; add one intercept in a regression if required. Do not include it in matching distance.',
        'same_rows_for_all_methods':True, 'row_key':['source_row','DUPERSID'],
        'categorical_model_input':'Use numeric_design file or factor coding below; never fit raw category codes as continuous covariates.',
        'rank_status':'No model fitting or numerical rank-based selection performed.'}
    model_rows, design_rows = [], []
    for row in complete:
        common = {k:row[k] for k in meta if k not in ['Z','Y','W']}
        common.update(Z=int(group(row)=='White'),Y=row['TOTEXP09'],W=row['SAQWT09F'])
        model = dict(common)
        model.update({f:row[f] for f in fields})
        model.update({f+'_label':labels[f][int(val(row,f))] for f in fields if spec['fields'][f]['type']=='categorical'})
        model_rows.append(model)
        design = dict(common)
        for sf in schema_fields:
            f = sf['source_field']
            if sf['type']=='continuous':
                if sf['columns']:design[f]=row[f]
            else:
                for level in sf['observed_levels']:
                    if level != sf['reference_level']:design[f+'__'+str(level)]=int(val(row,f)==level)
        assert list(k for k in design if k not in common) == columns
        design_rows.append(design)
    prefix='white_'+comparator.lower()
    model_cols=meta+fields+[f+'_label' for f in fields if spec['fields'][f]['type']=='categorical']
    write_csv(prefix+'_complete_records.csv.gz',model_rows,model_cols,True)
    write_csv(prefix+'_numeric_design.csv.gz',design_rows,meta+columns,True)
    for race in ['Overall','White',comparator]:
        rr=complete if race=='Overall' else [r for r in complete if group(r)==race]
        ww=[val(r,'SAQWT09F') for r in rr]
        sw=math.fsum(ww);sw2=math.fsum(w*w for w in ww)
        baseline.append({'pair':pair,'race_group':race,'field':'__population__','type':'sample',
            'level':'','label':'','n':len(rr),'sum_W':sw,'weighted_mean':'','weighted_sd':'',
            'level_n':'','level_sum_W':'','weighted_proportion':'','weight_ESS':sw*sw/sw2})
        for f in fields:
            if spec['fields'][f]['type']=='continuous':
                xx=[val(r,f) for r in rr]
                mean=math.fsum(w*x for w,x in zip(ww,xx))/sw
                sd=math.sqrt(math.fsum(w*(x-mean)**2 for w,x in zip(ww,xx))/sw)
                baseline.append({'pair':pair,'race_group':race,'field':f,'type':'continuous',
                    'level':'','label':dictionary[f]['label'],'n':len(rr),'sum_W':sw,
                    'weighted_mean':mean,'weighted_sd':sd,'level_n':'','level_sum_W':'',
                    'weighted_proportion':'','weight_ESS':''})
            else:
                for level in spec['fields'][f]['valid_levels']:
                    hits=[i for i,r in enumerate(rr) if val(r,f)==level]
                    levelw=math.fsum(ww[i] for i in hits)
                    baseline.append({'pair':pair,'race_group':race,'field':f,'type':'categorical',
                        'level':level,'label':labels[f][level],'n':len(rr),'sum_W':sw,
                        'weighted_mean':'','weighted_sd':'','level_n':len(hits),'level_sum_W':levelw,
                        'weighted_proportion':levelw/sw,'weight_ESS':''})
# The White complete-record rows must agree in both pairs under identical rules.
white_ids=[[r['DUPERSID'] for r in cohort_records[p] if group(r)=='White'] for p in cohort_records]
assert white_ids[0]==white_ids[1]
write_csv('cohort_flow.csv',flow)
write_csv('row_disposition_ledger.csv.gz',base_ledgers,compressed=True)
write_csv('covariate_exclusion_ledger.csv.gz',reasons,compressed=True)
write_csv('per_field_exclusions.csv',field_summary)
write_csv('weighted_baseline_summary.csv',baseline)
(OUT/'common_design_schema.json').write_text(json.dumps(schemas,indent=2)+'\n')
summary={'source_rows':len(rows), 'cohorts':{}, 'same_White_IDs_both_pairs':True,
 'original_individual_weights_preserved':True,'zero_outcomes_retained':True,
 'no_models_matching_or_effects':True,'exact_author_reproduction':False,
 'target':spec['target']}
for pair,complete in cohort_records.items():
 summary['cohorts'][pair]={'n':len(complete),'race_counts':dict(Counter(group(r) for r in complete)),
  'zero_expenditure_n':sum(val(r,'TOTEXP09')==0 for r in complete),
  'sum_W':math.fsum(val(r,'SAQWT09F') for r in complete),
  'design_columns_excluding_intercept':len(schemas[pair]['columns'])}
receipt={'started_utc':started,'finished_utc':datetime.now(timezone.utc).isoformat(),
 'frozen_spec_sha256':SPEC_SHA,'script_sha256':sha(Path(__file__)),
 'source_sha256':sha(SOURCE),
 'source_decoded_sha256':source_decoded_sha256,
 'accepted_source_gzip_sha256':spec['source']['sha256'],'summary':summary,
 'outputs':[{'path':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':sha(p)}
            for p in sorted(OUT.iterdir()) if p.is_file()]}
(OUT/'PREPROCESSING_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(summary,indent=2))
