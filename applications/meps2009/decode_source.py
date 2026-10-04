"""Decode official HC-129 selected fields and prepare source domains, without fitting.

Use explicit input paths and a new private output directory. The declared SAQ-adult
domain is not a claim to reproduce unpublished author preprocessing. Official
missing covariate codes remain unchanged until complete-record preparation.
"""
from pathlib import Path
from decimal import Decimal
from collections import Counter
import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import zipfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--raw-zip', type=Path, required=True,
                    help='Official HC-129 h129dat.zip ASCII archive')
parser.add_argument('--sas-program', type=Path, required=True,
                    help='Official HC-129 h129su.txt layout and formats')
parser.add_argument('--output-dir', type=Path, required=True,
                    help='New private directory for decoded data and receipt')
args = parser.parse_args()
if not __debug__:
    parser.error('Run without -O: source identity and format assertions are required.')
RAW_ZIP = args.raw_zip.resolve()
SAS_PROGRAM = args.sas_program.resolve()
OUT = args.output_dir.resolve()
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(RAW_ZIP) == "a1615e9d0b271e0a605dcdfac84091fc61c6e9b6cc577df9d960bc86fe2a1788"
assert sha(SAS_PROGRAM) == "84a9d48d443619cc85c79ed1b8b02faddaa4eb633bf803df669b6814e146e20e"
sas = SAS_PROGRAM.read_text()
OUT.mkdir(parents=True, exist_ok=False)
fields = """DUPERSID PANEL AGE31X AGE42X AGE53X AGE09X SEX MARRY31X MARRY42X MARRY53X MARRY09X
RACEX RACETHNX RACEWX RACEAX RACEBX HISPANX HISPCAT TOTEXP09 PERWT09F SAQWT09F SAQELIG SFFLAG42
PCS42 MCS42 BMINDX53 CHECK53 RTHLTH31 RTHLTH42 RTHLTH53 MNHLTH31 MNHLTH42 MNHLTH53
HIBPDX CHDDX ANGIDX MIDX OHRTDX STRKDX EMPHDX CHBRON31 CHBRON53 CHOLDX CANCERDX DIABDX ARTHDX ASTHDX ADHDADDX
IADLHP31 IADLHP42 IADLHP53 ADLHLP31 ADLHLP42 ADLHLP53 WLKLIM31 WLKLIM53 ACTLIM31 ACTLIM53 SOCLIM31 SOCLIM53 COGLIM31 COGLIM53
POVCAT09 EDUCYR HIDEG INSCOV09 REGION31 REGION42 REGION53 REGION09 VARSTR VARPSU""".split()
layouts = {name: (int(start), char == "$", int(width), int(dec))
           for start, name, char, width, dec in re.findall(
               r"^\s*@([0-9]+)\s+(\w+)\s+(\$?)([0-9]+)\.([0-9]+)\s*$", sas, re.M)}
labels = dict(re.findall(r"^\s*(\w+)\s*=\s*'([^']*)'\s*$", sas, re.M))
format_names = dict(re.findall(r"^\s+(\w+)\s+(\$?\w+)\.\s*$", sas, re.M))
format_blocks = {name: body.strip() for name, body in re.findall(
    r"^VALUE\s+(\$?\w+)[^\n]*\n(.*?)(?=^\s*;)", sas, re.M | re.S)}
assert len(fields) == len(set(fields)) and set(fields) <= set(layouts)
dictionary = {}
for name in fields:
    start, char, width, decimals = layouts[name]
    dictionary[name] = dict(start=start, end=start+width-1,
        type="CHAR" if char else "NUM", width=width, implied_decimals=decimals,
        label=labels.get(name), format_name=format_names.get(name),
        official_format=format_blocks.get(format_names.get(name, "")),
        codebook_url=f"https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_codebook.jsp?PUFId=H129&varName={name}")
assert all(v['label'] and v['official_format'] for v in dictionary.values())
(OUT / "source_selected_dictionary.json").write_text(json.dumps(dictionary, indent=2)+"\n")

def decode(line, name):
    start, char, width, decimals = layouts[name]
    value = line[start-1:start-1+width].decode("ascii").strip()
    if char:
        return value
    if not value:
        raise ValueError(f"Unexpected blank in {name}")
    parsed = Decimal(value)
    if "." not in value and decimals:
        parsed /= 10**decimals
    return int(parsed) if parsed == parsed.to_integral_value() else float(parsed)

def race(row):
    if row['HISPANX'] == 1: return 'Hispanic'
    if row['HISPANX'] == 2:
        return {1:'White',2:'Black',4:'Asian'}.get(row['RACEX'], 'Other_non_Hispanic')
    return 'Unclassified'

def gz_writer(path, columns):
    # Reproducible gzip container (no filename or wall-clock header).
    stream = io.TextIOWrapper(gzip.GzipFile(filename="", mode="wb",
        fileobj=path.open("wb"), mtime=0), encoding="utf-8", newline="")
    writer = csv.DictWriter(stream, fieldnames=columns)
    writer.writeheader()
    return stream, writer

rows = []
with zipfile.ZipFile(RAW_ZIP) as archive:
    assert archive.namelist() == ['h129.dat']
    with archive.open('h129.dat') as source:
        for row_number, line in enumerate(source, 1):
            line = line.rstrip(b'\r\n')
            assert len(line) == 5628
            row = {name: decode(line, name) for name in fields}
            row['source_row'] = row_number
            row['race_group'] = race(row)
            rows.append(row)
assert len(rows) == 36855
assert len({r['DUPERSID'] for r in rows}) == len(rows)
assert sum(r['PERWT09F'] > 0 for r in rows) == 34920
assert sum(r['SAQWT09F'] > 0 for r in rows) == 23171

stream, writer = gz_writer(OUT/'source_selected.csv.gz', ['source_row']+fields+['race_group'])
with stream:
    writer.writerows(rows)

def base_eligible(row):
    return row['AGE42X'] >= 18 and row['SAQWT09F'] > 0 and row['TOTEXP09'] >= 0

base = [r for r in rows if base_eligible(r)]
counts = {
    'source_rows':len(rows),
    'positive_PERWT09F':sum(r['PERWT09F'] > 0 for r in rows),
    'positive_SAQWT09F':sum(r['SAQWT09F'] > 0 for r in rows),
    'SAQ_adult_source_domain':len(base),
    'SAQ_adult_source_domain_races':dict(Counter(r['race_group'] for r in base)),
    'pair_domains':{},
    'missing_or_inapplicable_covariate_counts':{},
    'source_claimed_race_counts_not_a_filter':{'White':9830,'Asian':1446,'Black':4020,'Hispanic':5150},
    'exact_author_cohort_reproduced':False,
    'status':'Source-domain data preparation only; no covariate completion, model fit, estimate or inference.'
}
for group in ['Asian','Hispanic']:
    cohort=[r for r in base if r['race_group'] in ['White',group]]
    columns=['source_row','DUPERSID','race_group','Z','Y','W']+[f for f in fields if f!='DUPERSID']
    stream,writer=gz_writer(OUT/f'white_{group.lower()}_saq_adult_source_domain.csv.gz', columns)
    with stream:
        for row in cohort:
            writer.writerow(dict(row, Z=int(row['race_group']=='White'),
                                 Y=row['TOTEXP09'], W=row['SAQWT09F']))
    counts['pair_domains'][f'White-{group}']={'n':len(cohort),
        'group_counts':dict(Counter(r['race_group'] for r in cohort)),
        'outcome_zero_count':sum(r['TOTEXP09']==0 for r in cohort)}
    excluded={'DUPERSID','PANEL','RACEX','RACETHNX','RACEWX','RACEAX','RACEBX','HISPANX','HISPCAT',
              'TOTEXP09','PERWT09F','SAQWT09F','SAQELIG','SFFLAG42','VARSTR','VARPSU'}
    counts['missing_or_inapplicable_covariate_counts'][f'White-{group}']={
        f:dict(Counter(str(r[f]) for r in cohort if r[f]<0)) for f in fields if f not in excluded}

(OUT/'SOURCE_DOMAIN_RECEIPT.json').write_text(json.dumps(counts,indent=2)+"\n")
print(json.dumps({k:v for k,v in counts.items() if k!='missing_or_inapplicable_covariate_counts'},indent=2))
