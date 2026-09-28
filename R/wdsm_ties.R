# Recipient-specific uniform selection at the Mth-distance boundary.
# The tie stream is independent of sampling/replication and is restored on exit.
wdsm_make_matches <- function(A, D0, D1 = NULL, M = 5L,
    estimand = c("PATE", "PATT"), tie_seed = 20260917L,
    tie_rule = c("random", "row_order"), tie_tolerance = 64 * .Machine$double.eps) {
  estimand <- match.arg(estimand); tie_rule <- match.arg(tie_rule)
  n <- length(A)
  if (!all(A %in% 0:1) || length(M) != 1L || !is.finite(M) || M < 1 || M != as.integer(M) ||
      length(tie_seed) != 1L || !is.finite(tie_seed) || tie_seed < 0 ||
      tie_seed > .Machine$integer.max || tie_seed != as.integer(tie_seed) ||
      length(tie_tolerance) != 1L || !is.finite(tie_tolerance) || tie_tolerance < 0)
    stop("Invalid matching or tie configuration")
  rng_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(rng_kind))
    if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(tie_seed))
  result <- list(matches_0 = vector("list", n), matches_1 = vector("list", n))
  diagnostics <- list()
  for (arm in if (estimand == "PATE") 0:1 else 0L) {
    D <- if (arm == 0L) D0 else D1
    if (!is.matrix(D) || !identical(dim(D), c(as.integer(n), 2L)) || any(!is.finite(D)))
      stop("Matching requires two finite score coordinates per row")
    donors <- which(A == arm); recipients <- which(A != arm)
    if (length(donors) < M || !length(recipients)) stop("Insufficient opposite-arm donors")
    matches <- vector("list", n); usage <- integer(n)
    boundary_ties <- exact_boundary_ties <- randomized <- 0L
    maximum_boundary_group <- 1L
    for (i in recipients) {
      distance <- (D[donors, 1L] - D[i, 1L])^2 + (D[donors, 2L] - D[i, 2L])^2
      ordering <- order(distance, method = "radix")
      cutoff <- distance[ordering[M]]
      tolerance <- tie_tolerance * max(1, abs(cutoff))
      closer <- which(distance < cutoff - tolerance)
      tied <- which(abs(distance - cutoff) <= tolerance)
      needed <- M - length(closer)
      stopifnot(needed >= 1L, length(tied) >= needed)
      if (sum(distance == cutoff) > M - sum(distance < cutoff))
        exact_boundary_ties <- exact_boundary_ties + 1L
      if (length(tied) > needed) boundary_ties <- boundary_ties + 1L
      maximum_boundary_group <- max(maximum_boundary_group, length(tied))
      chosen <- if (tie_rule == "row_order") ordering[seq_len(M)] else {
        selected <- if (length(tied) > needed) {
          randomized <- randomized + 1L
          tied[sample.int(length(tied), needed, replace = FALSE)]
        } else tied
        positions <- c(closer, selected)
        positions[order(distance[positions], positions, method = "radix")]
      }
      matches[[i]] <- donors[chosen]
      usage[matches[[i]]] <- usage[matches[[i]]] + 1L
    }
    result[[paste0("matches_", arm)]] <- matches
    diagnostics[[as.character(arm)]] <- list(donor_arm = arm, recipients = length(recipients),
      boundary_tie_recipients = boundary_ties, exact_boundary_tie_recipients = exact_boundary_ties,
      randomized_recipients = randomized, maximum_boundary_group = maximum_boundary_group,
      donors_used = sum(usage[donors] > 0L), maximum_donor_reuse_count = max(usage[donors]),
      total_selected_pairs = sum(usage))
  }
  result$tie_diagnostics <- list(rule = tie_rule, seed = as.integer(tie_seed),
    tolerance = tie_tolerance, tolerance_scale = "max(1, Mth squared distance)",
    independent_recipient_selection = tie_rule == "random",
    original_match_sets_fixed_in_replication = TRUE, arms = diagnostics)
  result
}
