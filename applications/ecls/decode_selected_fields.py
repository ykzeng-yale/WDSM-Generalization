#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
"""Portable ECLS selected39 extraction/decoding; no analysis at import time.

Official TSV -> selected39 token CSV -> decoded/reason/companion views. The
official revised disability ZIP/DAT is loaded by exact character CHILDID. All
21409 original rows survive these steps; eligibility is a separate R operation.
"""
import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from decimal import Decimal
from pathlib import Path
import zipfile

OFFICIAL_TSV_MD5 = 'ebb102e230b0ebd64b506d152ad8972a'
OFFICIAL_RP_ZIP_SHA256 = '0004e7d4330e45eb2d2ef0b047cafe018992ac0d2d19d62c820949030bf9e636'
OFFICIAL_RP_MEMBER_SHA256 = 'e192801ca944703b08214a035aa491892cded6dd1b78615cd5e3eadf7317308d'
DEFAULT_POLICY = Path(__file__).resolve().with_name('selected39_variable_handling.csv')

STRINGS = {'CHILDID', 'S1_ID', 'S2_ID'}
MISSING = {-1: 'skip_not_applicable', -7: 'refused',
           -8: 'dont_know', -9: 'not_ascertained'}
ENUM = {
    'P2SPECND': (1, 2), 'GENDER': (1, 2), 'WKWHITE': (1, 2),
    'S2KPUPRI': (1, 2), 'P1EXPECT': tuple(range(1, 7)),
    'P1FIRKDG': (1, 2), 'P1HSEVER': (1, 2), 'FKCHGSCH': (0, 1),
    'S2KMINOR': tuple(range(1, 6)), 'P1FSTAMP': (1, 2),
    'P1HFAMIL': tuple(range(1, 6)), 'P1HPARNT': tuple(range(1, 10)),
    'WKCAREPK': (1, 2), 'P1HSCALE': tuple(range(1, 6)),
    'P1ATTENI': tuple(range(1, 5)), 'P1SOLVE': tuple(range(1, 5)),
    'P1PRONOU': tuple(range(1, 5)), 'P1DISABL': (1, 2),
    'P1WEIGHP': tuple(range(1, 14)), 'P1WEIGHO': tuple(range(0, 17)),
    'P1PREMAT': (1, 2), 'P1EARLY': tuple(range(1, 32)),
    'P1EARDAY': (1, 2), 'C1FMOTOR': tuple(range(0, 10)),
    'C1GMOTOR': tuple(range(0, 9)), 'P1AGEENT': tuple(range(54, 80))}
# These are substantive scale definitions, not rounded observed codebook ranges.
SCALE_BOUNDS = {'C1R4RSCL': (0, 212), 'C1R4MSCL': (0, 174),
                'C6R4MSCL': (0, 174), 'T1LEARN': (1, 4),
                'P1SADLON': (1, 4), 'P1IMPULS': (1, 4)}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for b in iter(lambda: stream.read(1024**2), b''):
            h.update(b)
    return h.hexdigest()


def decode(variable, token, missing_codes):
    token = token.strip()
    if token in ('', '.', 'NA', 'NaN'):
        return None, 'system_missing'
    if variable in STRINGS:
        return token, 'observed'
    value = Decimal(token)
    if not value.is_finite():
        return None, 'invalid_nonfinite'
    if value in missing_codes:
        return None, MISSING[int(value)]
    if value < 0 and variable != 'WKSESL':
        return None, 'undocumented_negative_code'
    if variable in ENUM and value not in ENUM[variable]:
        return None, 'invalid_documented_category'
    if variable in SCALE_BOUNDS:
        lo, hi = SCALE_BOUNDS[variable]
        if not lo <= value <= hi:
            return None, 'outside_substantive_scale'
    if variable == 'P1NUMSIB' and (value < 0 or value != value.to_integral_value()):
        return None, 'invalid_nonnegative_count'
    return token, 'observed'


def number(decoded, name):
    return None if decoded[name] is None else Decimal(decoded[name])


def companion(raw, decoded, disability):
    lead = number(decoded, 'P1PREMAT')
    amount = number(decoded, 'P1EARLY')
    unit = number(decoded, 'P1EARDAY')
    early_days = None
    payload_days = amount * (7 if unit == 1 else 1) if amount is not None and unit in (1, 2) else None
    flags = []
    if lead == 1:
        if payload_days is not None:
            early_days = payload_days
            status = 'reported_quantity_with_affirmative_lead'
            if early_days <= 14:
                flags.append('lead_yes_but_reported_days_le14')
        else:
            status = 'affirmative_lead_missing_number_or_unit'
            if raw['P1EARLY'] == '-1' or raw['P1EARDAY'] == '-1':
                flags.append('affirmative_lead_with_skipped_followup')
    elif lead == 2:
        status = 'not_more_than_two_weeks_early_exact_days_unknown'
        if not (raw['P1EARLY'] == '-1' and raw['P1EARDAY'] == '-1'):
            flags.append('negative_lead_with_nonstandard_followup')
    else:
        status = 'missing_lead_exact_days_unresolved'
        if payload_days is not None:
            flags.append('observed_followup_with_missing_lead')
    pounds = number(decoded, 'P1WEIGHP')
    ounces = number(decoded, 'P1WEIGHO')
    total_ounces = 16 * pounds + ounces if pounds is not None and ounces is not None else None
    identifier = raw['CHILDID']
    revised = disability.get(identifier)
    if revised is None:
        revised_value, revised_reason = None, 'no_ID_in_official_revised_sidecar'
    else:
        # Official RP1 member is keyed CHILDID/rp1disab, with no original column.
        # There is no original-P1DISABL fallback or positional alignment.
        revised_value, revised_reason = decode('P1DISABL', revised['rp1disab_official_revised'], {-1, -9})
    treatment = number(decoded, 'P2SPECND')
    weight = number(decoded, 'C1_6FC0')
    return {
        'CHILDID': identifier,
        'services_parent_report_kindergarten_schoolyear': None if treatment is None else int(treatment == 1),
        'wm_positive_supplied_weight_eligible': int(weight is not None and weight > 0),
        'birth_weight_total_ounces': total_ounces,
        'birth_weight_status': 'valid_two_components' if total_ounces is not None else 'component_missing_or_invalid',
        'preterm_more_than_two_weeks': None if lead is None else int(lead == 1),
        'prematurity_reported_days_affirmative_lead': early_days,
        'prematurity_number_unit_payload_days': payload_days,
        'prematurity_status': status,
        'prematurity_consistency_flags': '|'.join(flags) if flags else 'none',
        'parental_expectation_missing': int(decoded['P1EXPECT'] is None),
        'RP1DISAB_official_revised': revised_value,
        'RP1DISAB_missing_reason': revised_reason,
        'disability_original_and_revised_raw_codes_differ': None if revised is None else int(
            raw['P1DISABL'] != revised['rp1disab_official_revised'])}



def read_policy(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as stream:
        policies = list(csv.DictReader(stream))
    fields = [p['variable'] for p in policies]
    if len(fields) != 39 or len(set(fields)) != 39 or not STRINGS <= set(fields):
        raise ValueError('Exact selected39 policy required')
    missing = {p['variable']: {int(v) for v in p['explicit_negative_codes_listed_in_codebook'].split(';') if v}
               for p in policies}
    return fields, missing


def extract_selected_fields(source_tsv, output_csv, policy=DEFAULT_POLICY):
    """Single streaming TSV pass; no dense18949-column dataset and no row filter."""
    fields, _ = read_policy(policy)
    target = Path(output_csv)
    if target.exists() or not target.parent.is_dir():
        raise ValueError('Fresh selected39 CSV in an existing local directory required')
    md5 = hashlib.md5()
    count = 0
    seen = set()
    with Path(source_tsv).open('rb') as stream, target.open('x', newline='', encoding='utf-8') as sink:
        def lines():
            first = True
            for raw in stream:
                md5.update(raw)
                yield raw.decode('utf-8-sig' if first else 'utf-8')
                first = False
        reader = csv.DictReader(lines(), delimiter='\t')
        if len(reader.fieldnames or []) != 18949 or len(set(reader.fieldnames)) != 18949 or not set(fields) <= set(reader.fieldnames):
            raise ValueError('Official DS0001 TSV header differs')
        writer = csv.DictWriter(sink, fields)
        writer.writeheader()
        for row in reader:
            if None in row or any(row[v] is None for v in fields):
                raise ValueError('TSV row has malformed column count')
            identifier = row['CHILDID']
            if not identifier or identifier != identifier.strip() or identifier in seen:
                raise ValueError('Nonempty unique exact character CHILDID required')
            seen.add(identifier)
            writer.writerow({v: row[v] for v in fields})
            count += 1
    if count != 21409 or md5.hexdigest() != OFFICIAL_TSV_MD5:
        raise ValueError('Official21409-row TSV identity differs; selected file retained, no success receipt')
    return dict(status='SELECTED39_TOKEN_EXTRACTION_COMPLETE_DATA_ONLY', source_md5=md5.hexdigest(),
                rows=count, fields=39, policy_sha256=sha(Path(policy)),
                selected_csv_sha256=sha(target), no_row_filter=True,
                no_numeric_ID_cast=True, no_models_matching_RNG=True)


def read_revised_disability(path):
    """Read official NCES2007031 ZIP or extracted RP1DISAB.dat; no ID casting."""
    source = Path(path)
    if source.suffix.lower() == '.zip':
        if sha(source) != OFFICIAL_RP_ZIP_SHA256:
            raise ValueError('Official revised-disability ZIP SHA differs')
        with zipfile.ZipFile(source) as archive:
            info = archive.getinfo('RP1DISAB.dat')
            if info.file_size != 255167 or info.CRC != 0x847b993c:
                raise ValueError('Official RP1 member size/CRC differs')
            payload = archive.read(info)  # full read verifies ZIP CRC; only this small member
    else:
        if source.name != 'RP1DISAB.dat':
            raise ValueError('Supply official ZIP or exact RP1DISAB.dat member')
        payload = source.read_bytes()
    member_sha = hashlib.sha256(payload).hexdigest()
    if member_sha != OFFICIAL_RP_MEMBER_SHA256:
        raise ValueError('Official revised RP1 member SHA differs')
    reader = csv.DictReader(io.StringIO(payload.decode('ascii')), delimiter='\t')
    if reader.fieldnames != ['CHILDID', 'rp1disab']:
        raise ValueError('Official revised RP1 header differs')
    revised = {}
    for row in reader:
        identifier = row['CHILDID']
        if not identifier or identifier != identifier.strip() or identifier in revised or None in row:
            raise ValueError('Official revised table must have unique exact CHILDIDs')
        revised[identifier] = {'rp1disab_official_revised': row['rp1disab']}
    if len(revised) != 21260:
        raise ValueError('Official revised-disability row count differs')
    return revised, dict(source_sha256=sha(source), member='RP1DISAB.dat', member_sha256=member_sha,
                         rows=21260, actual_field='rp1disab', document_label='RP1DISABL')


def decode_selected_fields(tokens_csv, disability_source, output_directory,
                           tokens_sha256, policy=DEFAULT_POLICY):
    fields, missing = read_policy(policy)
    if sha(Path(tokens_csv)) != tokens_sha256:
        raise ValueError('Selected39 token input SHA differs from supplied extraction provenance')
    revised, rp_provenance = read_revised_disability(disability_source)
    target = Path(output_directory)
    if target.exists() or not target.parent.is_dir():
        raise ValueError('Fresh local decoder destination required')
    target.mkdir()
    names = {'decoded': 'codebook_decoded_selected39.csv', 'reasons': 'selected39_missing_reasons.csv',
             'companions': 'scientific_companion_derivations.csv'}
    field_status = {v: Counter() for v in fields}
    patterns = Counter()
    seen = set()
    receipt = dict(status='DECODING_IN_PROGRESS', tokens_sha256=tokens_sha256,
                   policy_sha256=sha(Path(policy)), revised_disability=rp_provenance,
                   no_models_matching_RNG=True, no_weight_estimation=True, scientific_acceptance=False)
    try:
        with Path(tokens_csv).open(newline='', encoding='utf-8-sig') as stream, \
                (target/names['decoded']).open('x', newline='') as a, \
                (target/names['reasons']).open('x', newline='') as b, \
                (target/names['companions']).open('x', newline='') as c:
            reader = csv.DictReader(stream)
            if reader.fieldnames != fields:
                raise ValueError('Selected39 token names/order differ from policy')
            decoded_writer = csv.DictWriter(a, fields); decoded_writer.writeheader()
            reason_writer = csv.DictWriter(b, fields); reason_writer.writeheader()
            companion_writer = None
            for raw in reader:
                identifier = raw['CHILDID']
                if not identifier or identifier != identifier.strip() or identifier in seen or None in raw:
                    raise ValueError('Unique exact all-row CHILDID binding required')
                seen.add(identifier)
                decoded, reasons = {}, {}
                for variable in fields:
                    decoded[variable], reasons[variable] = decode(variable, raw[variable], missing[variable])
                    field_status[variable][reasons[variable]] += 1
                if decoded['CHILDID'] != identifier:
                    raise ValueError('Character CHILDID changed')
                decoded_writer.writerow(decoded)
                reason_writer.writerow(dict(reasons, CHILDID=identifier))
                derived = companion(raw, decoded, revised)
                if companion_writer is None:
                    companion_writer = csv.DictWriter(c, list(derived)); companion_writer.writeheader()
                companion_writer.writerow(derived)
                for k in ('birth_weight_status', 'prematurity_status', 'prematurity_consistency_flags', 'RP1DISAB_missing_reason'):
                    patterns[f'{k}:{derived[k]}'] += 1
        if len(seen) != 21409 or not set(revised) <= seen:
            raise ValueError('All-row21409 / official21260 ID join differs')
        receipt.update(status='DECODED_ALL_ROWS_COMPANION_VIEWS_REVIEW_REQUIRED', records_retained=len(seen),
                       fields_decoded=39, revised_IDs_matched=len(set(revised)&seen),
                       revised_IDs_unmatched=len(seen-set(revised)), field_reason_counts=field_status,
                       companion_pattern_counts=patterns, negative_SES_preserved=True,
                       unknown_prematurity_not_zero=True, supplied_C1_6FC0_unchanged=True,
                       no_original_disability_fallback=True, no_school_means_or_cohort_selection=True,
                       outputs={k: dict(file=name, sha256=sha(target/name), bytes=(target/name).stat().st_size)
                                for k,name in names.items()})
    except BaseException as error:
        receipt.update(status='DECODING_FAILED_PARTIAL_OUTPUTS_RETAINED', error=repr(error))
        raise
    finally:
        with (target/'semantic_decoding_receipt.json').open('x') as stream:
            json.dump(receipt, stream, indent=2); stream.write('\n')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    extract = sub.add_parser('extract')
    extract.add_argument('--tsv', required=True); extract.add_argument('--selected-csv', required=True)
    extract.add_argument('--receipt', required=True)
    decoding = sub.add_parser('decode')
    decoding.add_argument('--selected-csv', required=True); decoding.add_argument('--tokens-sha256', required=True)
    decoding.add_argument('--revised-disability', required=True); decoding.add_argument('--output', required=True)
    args = parser.parse_args()
    if args.action == 'extract':
        receipt_path = Path(args.receipt)
        if receipt_path.exists() or not receipt_path.parent.is_dir():
            raise ValueError('Fresh extraction receipt required')
        result = extract_selected_fields(args.tsv, args.selected_csv)
        with receipt_path.open('x') as stream: json.dump(result, stream, indent=2); stream.write('\n')
    else:
        result = decode_selected_fields(args.selected_csv, args.revised_disability,
                                        args.output, args.tokens_sha256)
    print(json.dumps(dict(status=result['status'], no_models_matching_RNG=True)))


if __name__ == '__main__':
    main()
