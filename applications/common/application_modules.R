# SPDX-License-Identifier: GPL-3.0-only
# Companion to the extracted application functions; no work runs at source time.
# Loading this provider explicitly checks the installed WeightedMatching API. It does
# not install packages, fit models, create counts, or choose a private library.
wm_application_modules <- function() {
  if (!requireNamespace("WeightedMatching", quietly = TRUE))
    stop("Install the validated WeightedMatching source before using this application module.")
  if (!requireNamespace("digest", quietly = TRUE))
    stop("The extracted count-column checker requires the digest package.")
  ns <- asNamespace("WeightedMatching")
  required <- c(".wm_ns_ols", ".wm_ns_basis", "wm_match",
                "wm_bootstrap_refit", "wm_wdsm_fit")
  available <- vapply(required, function(name)
    exists(name, envir = ns, inherits = FALSE) &&
      is.function(get(name, envir = ns, inherits = FALSE)), logical(1))
  if (!all(available))
    stop("Installed WeightedMatching lacks required application functions: ",
         paste(required[!available], collapse = ", "))
  formals_required <- list(
    wm_match = c("Y", "Z", "weights", "scores0", "scores1", "M",
                 "estimand", "mean0", "mean1", "variance", "tie_rule", "tie_seed"),
    wm_bootstrap_refit = c("object", "counts", "refit", "conf.level"),
    wm_wdsm_fit = c("Y", "Z", "weights", "ps_design", "pg0_design", "pg1_design",
                   "M", "estimand", "inference", "ps_weighting", "tie_rule", "tie_seed"),
    .wm_ns_ols = c("x", "y", "w", "label", "tolerance"),
    .wm_ns_basis = "s")
  for (name in names(formals_required)) {
    missing <- setdiff(formals_required[[name]],
                      names(formals(get(name, envir = ns, inherits = FALSE))))
    if (length(missing))
      stop("Installed WeightedMatching API differs for ", name, ": ",
           paste(missing, collapse = ", "))
  }
  if (!exists("wdsm_case_fit_propensity", mode = "function", inherits = TRUE))
    stop("Source case_propensity_solver.R in the same environment first.")
  list(source = list(wdsm_case_fit_propensity = wdsm_case_fit_propensity),
       common = ns)
}
