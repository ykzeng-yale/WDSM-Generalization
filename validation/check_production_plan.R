# Deterministic, source-only production-plan checks. No simulations or jobs.
# Rscript validation/check_production_plan.R
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this validation with Rscript")
script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
simulation_dir <- file.path(dirname(script), "..", "simulations")
source(file.path(simulation_dir, "production.R"))
config <- utils::read.csv(file.path(simulation_dir, "production_config.csv"), stringsAsFactors = FALSE)
plan <- wm_production_plan(config)
stopifnot(nrow(plan) == 244L, identical(plan$task_id, 1:244),
  !anyDuplicated(plan$result_file), sum(plan$replications) == 48000L,
  sum(plan$requested_records) == 576000L,
  all(plan$requested_records == 12 * plan$replications),
  identical(which(plan$stage == "timing"), 1:4),
  all(plan$first[1:4] == 1), all(plan$last[1:4] == 10),
  all(plan$n[1:4] == 2000), all(plan$M[1:4] == 3),
  setequal(plan$d[1:4], c(1, 5)), setequal(plan$design[1:4], c("strong", "weak")))
for (scenario in config$scenario) {
  p <- plan[plan$scenario == scenario, , drop = FALSE]
  ids <- unlist(lapply(seq_len(nrow(p)), function(i) seq.int(p$first[i], p$last[i])), use.names = FALSE)
  stopifnot(identical(sort(as.integer(ids)), 1:1000), !anyDuplicated(ids),
    sum(p$requested_records) == 12000L,
    all(p$seed == config$seed[config$scenario == scenario]))
}
stopifnot(sum(plan$n == 200) == 16L, sum(plan$n == 800) == 64L,
  sum(plan$n == 2000) == 164L,
  all(plan$replications[plan$n == 200] == 1000),
  all(plan$replications[plan$n == 800] == 250),
  sum(plan$replications == 90) == 4L)
fails <- function(expression, pattern) {
  message <- tryCatch({force(expression); NULL}, error = conditionMessage)
  if (is.null(message) || !grepl(pattern, message)) stop("Expected error: ", pattern)
}
fails(wm_production_plan(config[-1, ]), "48-scenario")
bad <- config; bad$seed[2] <- bad$seed[1]
fails(wm_production_plan(bad), "unique")
bad <- config; bad$n[1] <- 201
fails(wm_production_plan(bad), "full primary grid")
bad <- config; bad$scenario[1] <- "misleading_name"
fails(wm_production_plan(bad), "exact configuration")
bad <- config; bad$seed[1] <- 1.5
fails(wm_production_plan(bad), "integer")
cat("Production plan checks passed: 244 disjoint tasks, 48,000 datasets, 576,000 records; no work was launched.\n")
