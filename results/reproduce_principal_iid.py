#!/usr/bin/env python3
"""Reaggregate the published principal iid CSV records, without fits or RNG.

Uses the existing formal_metrics implementation and the archived fixed target
metadata. This verifies reporting arithmetic; it does not recreate raw samples,
fits, exact roots, or the truth quadrature. No private execution files are needed.
"""
import csv
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys


def need(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csv_rows(path):
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def main(argv):
    need(len(argv) == 1, "Usage: reproduce_principal_iid.py NEW_OUTPUT_DIRECTORY")
    output = Path(argv[0]).expanduser().resolve()
    need(not output.exists() and output.parent.is_dir(), "Use a fresh output directory under an existing parent")
    root = Path(__file__).resolve().parents[1]
    bundle = root / "results/principal_iid"
    manifest_path = root / "results/publication_results_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    pins = [x for x in manifest["files"] if x["file"].startswith("results/principal_iid/")
            and x["file"].endswith((".csv", ".csv.gz"))]
    need(len({x["file"] for x in pins}) == len(pins), "Duplicate archive input")
    for pin in pins:
        need(sha(root / pin["file"]) == pin["sha256"], "Published artifact changed: " + pin["file"])
    metric_path = root / "simulations/principal_iid/formal_metrics.py"
    inventory_path = metric_path.with_name("source_inventory.json")
    inventory = json.loads(inventory_path.read_text())
    metric_pin = [x for x in inventory["files"] if x["path"] == str(metric_path.relative_to(root))]
    need(len(metric_pin) == 1 and sha(metric_path) == metric_pin[0]["sha256"], "Reporting-source identity changed")
    loader = importlib.util.spec_from_file_location("principal_archive_metrics", metric_path)
    metrics = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(metrics)
    _, records = csv_rows(bundle / "bound_records.csv.gz")
    integer_fields = ("group_index", "n", "replicate", "data_seed", "count_seed", "M", "B", "p")
    number_fields = ("target", "estimate", "raw_estimate", "analytic_variance", "analytic_lower",
                     "analytic_upper", "refit_variance", "refit_lower", "refit_upper")
    boolean_fields = ("point_available", "analytic_available", "contribution_available", "refit_available")
    cases = {}
    for row in records:
        for field in integer_fields:
            row[field] = int(row[field])
        for field in number_fields:
            row[field] = float(row[field]) if row[field] else None
            need(row[field] is None or math.isfinite(row[field]), "Nonfinite recorded number")
        for field in boolean_fields:
            need(row[field] in ("True", "False"), "Invalid availability flag")
            row[field] = row[field] == "True"
        for field in ("count_fit_failures", "count_draw_failures"):
            row[field] = json.loads(row[field])
            ix = row[field]
            need(isinstance(ix, list) and all(type(i) is int and 1 <= i <= 200 for i in ix)
                 and ix == sorted(set(ix)), "Invalid count failure positions")
        g, n, rep, overlap = (row[k] for k in ("group_index", "n", "replicate", "overlap"))
        need((g, overlap, n) in metrics.GROUPS and 1 <= rep <= metrics.R, "Unexpected sample key")
        need(row["family"] in metrics.FAMILIES and row["estimand"] in metrics.ESTIMANDS, "Unexpected method")
        need(row["M"] == 3 and row["B"] == 200, "Changed M or B")
        p = {"PS": (26, 18), "DSM": (42, 26), "X6": (68, 40)}[row["family"]][row["estimand"] == "PATT"]
        need(row["p"] == p, "Changed complete nuisance dimension")
        need(row["case_id"] == f"g{g}_{overlap}_n{n}_rep{rep:04d}", "Changed case identifier")
        need(row["data_seed"] == 104000000 + 100000*g + rep and
             row["count_seed"] == 204000000 + 100000*g + rep, "Changed prescribed seed")
        need(row["refit_variance_divisor"] == "B" and row["empirical_variance_divisor"] == "R-1", "Changed variance definition")
        need(not row["refit_available"] or not(row["count_fit_failures"] or row["count_draw_failures"]),
             "An available refit conceals a failed column")
        need(not row["contribution_available"] or row["analytic_available"], "Inconsistent contribution availability")
        for kind in ("analytic", "refit"):
            if row[kind + "_available"]:
                vv, lo, hi = (row[kind + "_" + x] for x in ("variance", "lower", "upper"))
                need(row["point_available"] and all(x is not None for x in (vv, lo, hi))
                     and vv >= 0 and lo <= hi, "Invalid available interval")
        case = cases.setdefault((g, rep), {})
        key = (row["family"], row["estimand"])
        need(key not in case, "Duplicate method row")
        case[key] = row
    expected_methods = {(f, e) for f in metrics.FAMILIES for e in metrics.ESTIMANDS}
    need(set(cases) == {(g, r) for g, _, _ in metrics.GROUPS for r in range(1, metrics.R+1)}, "Incomplete sample universe")
    for case in cases.values():
        need(set(case) == expected_methods, "Missing method row")
        need(len({(r["sample_sha256"], r["counts_sha256"]) for r in case.values()}) == 1, "Methods are not paired")

    # Only fixed target metadata are read here; reported performance is recomputed.
    _, sensitivity = csv_rows(bundle / "target_sensitivity.csv")
    truth_values = {}
    for row in sensitivity:
        key = (row["overlap"], row["estimand"])
        slot = truth_values.setdefault(key, {})
        for field, value in (("target", float(row["target"])), ("truth_method", row["truth_method"]),
                             ({-1e-6: "sensitivity_lower", 0.: "target", 1e-6: "sensitivity_upper"}
                              [float(row["truth_offset"])], float(row["perturbed_target"]))):
            need(field not in slot or slot[field] == value, "Inconsistent fixed target metadata")
            slot[field] = value
    need(set(truth_values) == {(o, e) for _, o, _ in metrics.GROUPS for e in metrics.ESTIMANDS}, "Missing target")
    reports = metrics.summarize(records, truth_values)
    comparisons = {}
    saved_headers = {}
    for name, fresh in reports.items():
        fields, saved = csv_rows(bundle / (name + ".csv"))
        saved_headers[name] = fields
        need(len(fresh) == len(saved), "Changed row count: " + name)
        largest = 0.
        for actual, expected in zip(fresh, saved):
            need(set(actual) == set(expected), "Changed fields: " + name)
            for field, value in actual.items():
                old = expected[field]
                if value is None:
                    need(old == "", "Changed missing value: " + name + "/" + field)
                elif isinstance(value, bool):
                    need(old == str(value), "Changed Boolean: " + name + "/" + field)
                elif isinstance(value, (int, float)):
                    gap = abs(value - float(old)); largest = max(largest, gap)
                    need(gap <= 1e-12 + 1e-10*abs(float(old)), "Changed numeric metric: " + name + "/" + field)
                elif isinstance(value, (dict, list)):
                    need(json.loads(old) == value, "Changed structured field")
                else:
                    need(old == value, "Changed label: " + name + "/" + field)
        comparisons[name] = dict(rows=len(fresh), maximum_absolute_numeric_difference=largest)
    output.mkdir()
    for name, fresh in reports.items():
        with (output / (name + ".csv")).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=saved_headers[name]); writer.writeheader()
            for row in fresh:
                writer.writerow({k: json.dumps(v, sort_keys=True, allow_nan=False) if isinstance(v, (dict, list)) else v
                                 for k, v in row.items()})
    need(all(sha(root / x["file"]) == x["sha256"] for x in pins), "Archive changed during reporting")
    need(sha(metric_path) == metric_pin[0]["sha256"], "Reporting source changed during execution")
    receipt = dict(status="COMPLETE_ARCHIVE_REAGGREGATION", cases=len(cases), method_records=len(records),
                   comparisons=comparisons, inputs=pins, reporting_source=metric_pin[0],
                   script_sha256=sha(Path(__file__)),
                   scope="Reporting arithmetic from published synthetic records and fixed target metadata; no regeneration, fitting, inference recomputation or proof of population assumptions.")
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: receipt[k] for k in ("status", "cases", "method_records", "comparisons")}, indent=2))


if __name__ == "__main__":
    main(sys.argv[1:])
