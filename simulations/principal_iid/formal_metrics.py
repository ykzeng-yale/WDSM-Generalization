"""Pure reporting arithmetic for the predeclared complete principal iid study.

No files, fitting, simulation or RNG. Definitions transcribe existing project
paired_comparison_metrics.R and summarize_saved_augmentation.R. None values are
undefined/withheld; a failed interval has operational hit zero, not variance zero.
"""
import math

R = 1000
FAMILIES = ("PS", "DSM", "X6")
ESTIMANDS = ("PATE", "PATT")
GROUPS = ((1, "GoodOverlap", 1000), (2, "GoodOverlap", 5000),
          (3, "PoorOverlap", 1000), (4, "PoorOverlap", 5000))


def average(x):
    return math.fsum(x) / len(x) if x else None


def variance(x):
    if len(x) < 2:
        return None
    m = average(x)
    return math.fsum((v - m) ** 2 for v in x) / (len(x) - 1)


def mcse(x):
    v = variance(x)
    return math.sqrt(v / len(x)) if v is not None else None


def bernoulli_se(rate, denominator):
    return math.sqrt(rate * (1 - rate) / denominator)


def point_values(x, target):
    if not x:
        return dict(mean_estimate=None, signed_bias=None, signed_RB_percent=None,
                    bias_mcse=None, RB_mcse=None, empirical_variance=None,
                    empirical_sd=None, EV_mcse=None, mse=None, rmse=None,
                    rmse_mcse=None)
    mean = average(x)
    bias = mean - target
    ev = variance(x)
    squared = [(v - target) ** 2 for v in x]
    mse = average(squared)
    rmse = math.sqrt(mse)
    q = [len(x) / (len(x) - 1) * (v - mean) ** 2 for v in x] if len(x) > 1 else []
    return dict(mean_estimate=mean, signed_bias=bias,
                signed_RB_percent=100 * bias / abs(target) if target != 0 else None,
                bias_mcse=mcse(x),
                RB_mcse=100 * mcse(x) / abs(target) if target != 0 and len(x) > 1 else None,
                empirical_variance=ev, empirical_sd=math.sqrt(ev) if ev is not None else None,
                EV_mcse=mcse(q), mse=mse, rmse=rmse,
                rmse_mcse=mcse([v / (2 * rmse) for v in squared]) if rmse > 0 else None)


def interval_arrays(rows, kind, target):
    available = [bool(r[kind + "_available"]) for r in rows]
    hit = [int(ok and r[kind + "_lower"] <= target <= r[kind + "_upper"])
           for r, ok in zip(rows, available)]
    reported = [r[kind + "_variance"] for r, ok in zip(rows, available) if ok]
    width = [r[kind + "_upper"] - r[kind + "_lower"]
             for r, ok in zip(rows, available) if ok]
    return available, hit, reported, width


def finite_output(x):
    if isinstance(x, float) and not math.isfinite(x):
        raise ArithmeticError("Nonrepresentable reporting arithmetic; no clipping or replacement")
    if isinstance(x, dict):
        for y in x.values():
            finite_output(y)
    elif isinstance(x, (tuple, list)):
        for y in x:
            finite_output(y)


def summarize(records, truths):
    """Input is 24000 fully authenticated point rows, never a partial sample."""
    bycell = {}
    for row in records:
        key = (row["group_index"], row["estimand"], row["family"])
        bycell.setdefault(key, []).append(row)
    expected = {(g, e, f) for g, _, _ in GROUPS for e in ESTIMANDS for f in FAMILIES}
    if set(bycell) != expected or len(records) != 24000:
        raise ValueError("Exact24 point-cell universe required before metrics")
    points, diagnostics, inference, paired, sensitivity, failures = [], [], [], [], [], []
    stored = {}
    for g, overlap, n in GROUPS:
        for estimand in ESTIMANDS:
            target_record = truths[(overlap, estimand)]
            target = target_record["target"]
            for family in FAMILIES:
                key = (g, estimand, family)
                x = sorted(bycell[key], key=lambda r: r["replicate"])
                if [r["replicate"] for r in x] != list(range(1, R + 1)):
                    raise ValueError("Missing/duplicate prescribed replicate")
                if any(r["target"] != target for r in x):
                    raise ValueError("Cell target is not the authenticated fixed truth")
                common = dict(group_index=g, overlap=overlap, n=n, estimand=estimand,
                              family=family, matching_dimension={"PS": 1, "DSM": 2, "X6": 6}[family],
                              M=3, B=200, p=x[0]["p"], requested=R, received=R,
                              target=target, truth_method=target_record["truth_method"])
                available_x = [r["estimate"] for r in x if r["point_available"]]
                complete_points = len(available_x) == R
                pv = point_values(available_x, target)
                full = dict(pv) if complete_points else {k: None for k in pv}
                point = dict(common, point_available=len(available_x), point_failures=R-len(available_x),
                             point_metric_scope="all_requested" if complete_points else "withheld_point_failures",
                             relative_bias_status="withheld_point_failures" if not complete_points else
                             "undefined_zero_target" if target == 0 else "defined", **full)
                diagnostic = dict(common, diagnostic_scope="available_points_only_not_full_estimator_performance",
                                  denominator=len(available_x), **pv)
                points.append(point)
                diagnostics.append(diagnostic)
                q = [R/(R-1)*(r["estimate"]-pv["mean_estimate"])**2 for r in x] if complete_points else None
                ev = full["empirical_variance"]
                bykind = {}
                for kind in ("analytic", "refit"):
                    av, hit, vv, width = interval_arrays(x, kind, target)
                    all_intervals = all(av)
                    coverage = sum(hit)/R
                    availability = sum(av)/R
                    full_mean_var = average(vv) if all_intervals else None
                    calibration = full_mean_var/ev if all_intervals and complete_points and ev > 0 else None
                    calibration_mcse = mcse([(v-calibration*qr)/ev for v, qr in zip(vv,q)]) if calibration is not None else None
                    summary = dict(point, inference=kind,
                        reported_variance_kind="complete_fitted_row" if kind == "analytic" else "fixed_reuse_prediction_refit",
                        empirical_variance_divisor="R-1", refit_variance_divisor="B" if kind == "refit" else "not_applicable",
                        interval_available=sum(av), interval_failures=R-sum(av),
                        availability=availability, availability_mcse=bernoulli_se(availability,R),
                        coverage_rule="inclusive", operational_coverage=coverage,
                        coverage_mcse=bernoulli_se(coverage,R), operational_denominator=R,
                        successful_intervals_coverage=sum(hit)/sum(av) if any(av) else None,
                        successful_intervals_denominator=sum(av),
                        mean_reported_variance=full_mean_var,
                        mean_reported_variance_mcse=mcse(vv) if all_intervals else None,
                        reported_variance_scope="all_requested" if all_intervals else "withheld_interval_failures",
                        mean_reported_variance_available_only=average(vv),
                        mean_interval_width=average(width) if all_intervals else None,
                        mean_interval_width_available_only=average(width),
                        mean_variance_over_EV=calibration, variance_calibration_mcse=calibration_mcse,
                        calibration_status="defined" if calibration is not None else
                        "withheld_point_or_interval_failures" if not(complete_points and all_intervals) else "undefined_zero_EV",
                        contribution_available=sum(r["contribution_available"] for r in x),
                        contribution_MC_is_separate_sampling_variance_method=False)
                    inference.append(summary)
                    bykind[kind] = (av, hit, vv)
                    # The frozen truth artifact supplies this +/-1e-6 arithmetic display
                    # for BOTH targets. PATE remains analytically known, without GH error.
                    for offset, shifted in ((-1e-6,target_record["sensitivity_lower"]),
                                            (0.0,target), (1e-6,target_record["sensitivity_upper"])):
                        _, hs, _, _ = interval_arrays(x,kind,shifted)
                        shifted_cov = sum(hs)/R
                        rb = 100*(pv["mean_estimate"]-shifted)/abs(shifted) if complete_points and shifted != 0 else None
                        sensitivity.append(dict(common,inference=kind,truth_offset=offset,
                            perturbed_target=shifted,operational_coverage=shifted_cov,
                            coverage_mcse=bernoulli_se(shifted_cov,R),operational_denominator=R,
                            signed_RB_percent=rb,point_metric_scope=point["point_metric_scope"],
                            coverage_changes_from_primary=sum(a!=b for a,b in zip(hs,hit)),
                            sensitivity_scope="arithmetic_display_analytic_PATE_no_quadrature_uncertainty" if estimand=="PATE" else
                            "operational_numerical_target_sensitivity_not_rigorous_enclosure",
                            primary_target_unchanged=True))
                aa, ah, avar = bykind["analytic"]
                ra, rh, rvar = bykind["refit"]
                differences = [a-b for a,b in zip(ah,rh)]
                ratio = average(avar)/average(rvar) if all(aa) and all(ra) and average(rvar)>0 else None
                pair = dict(common, paired_cases=R,
                    same_points=True, pairing="same sample_sha256 and counts_sha256 within case",
                    both_available=sum(a and b for a,b in zip(aa,ra)),
                    analytic_only_available=sum(a and not b for a,b in zip(aa,ra)),
                    refit_only_available=sum(b and not a for a,b in zip(aa,ra)),
                    neither_available=sum(not a and not b for a,b in zip(aa,ra)),
                    both_cover=sum(a and b for a,b in zip(ah,rh)),
                    analytic_only_cover=sum(a and not b for a,b in zip(ah,rh)),
                    refit_only_cover=sum(b and not a for a,b in zip(ah,rh)),
                    neither_cover=sum(not a and not b for a,b in zip(ah,rh)),
                    coverage_difference_analytic_minus_refit=average(differences),
                    paired_coverage_difference_mcse=mcse(differences),
                    ratio_of_mean_variances_analytic_over_refit=ratio,
                    paired_ratio_mcse=mcse([(a-ratio*b)/average(rvar) for a,b in zip(avar,rvar)]) if ratio is not None else None,
                    mean_variance_difference=average([a-b for a,b in zip(avar,rvar)]) if all(aa) and all(ra) else None,
                    paired_variance_difference_mcse=mcse([a-b for a,b in zip(avar,rvar)]) if all(aa) and all(ra) else None,
                    paired_variance_status="defined" if ratio is not None else
                    "withheld_interval_failures" if not(all(aa) and all(ra)) else "undefined_zero_mean_refit_variance",
                    finite_B_variation_retained=True)
                paired.append(pair)
                stored[key] = (x,point,q)
                for row in x:
                    if not all(row[z] for z in ("point_available","analytic_available","contribution_available","refit_available")):
                        failures.append(dict(row))
    # Paired empirical efficiency uses WM-PS within the SAME overlap/n/estimand/M.
    for key,(x,point,q) in stored.items():
        ref_x,ref_point,ref_q = stored[(key[0],key[1],"PS")]
        if [(r["sample_sha256"],r["counts_sha256"]) for r in x] != [(r["sample_sha256"],r["counts_sha256"]) for r in ref_x]:
            raise ValueError("Empirical-RE reference rows are not paired by actual inputs")
        ev = point["empirical_variance"]
        ev_ref = ref_point["empirical_variance"]
        eligible = point["point_metric_scope"]==ref_point["point_metric_scope"]=="all_requested"
        re = ev_ref/ev if eligible and ev>0 and ev_ref>0 else None
        extra = dict(reference_family="PS",reference_empirical_variance=ev_ref,
                     relative_efficiency=re,
                     relative_efficiency_mcse=mcse([(a-re*b)/ev for a,b in zip(ref_q,q)]) if re is not None else None,
                     relative_efficiency_status="defined" if re is not None else
                     "withheld_paired_point_failures" if not eligible else "undefined_zero_EV")
        point.update(extra)
        for item in inference:
            if (item["group_index"],item["estimand"],item["family"])==key:
                item.update(extra)
    # Wide24-row main table and long48-row figure data keep denominators/failures.
    main=[]
    for point in points:
        out=dict(point)
        for kind in ("analytic","refit"):
            inf=next(x for x in inference if (x["group_index"],x["estimand"],x["family"],x["inference"])==
                     (point["group_index"],point["estimand"],point["family"],kind))
            for field in ("interval_available","interval_failures","availability","operational_coverage",
                          "coverage_mcse","mean_reported_variance","mean_variance_over_EV","variance_calibration_mcse",
                          "reported_variance_scope","calibration_status"):
                out[kind+"_"+field]=inf[field]
        main.append(out)
    result=dict(points=points,point_success_diagnostics=diagnostics,inference=inference,paired=paired,
                target_sensitivity=sensitivity,failures=failures,main_table=main,figure_data=inference)
    if not (len(points)==24 and len(inference)==48 and len(paired)==24 and len(sensitivity)==144):
        raise ValueError("Unexpected complete-study output cell universe")
    finite_output(result)
    return result
