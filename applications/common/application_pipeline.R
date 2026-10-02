# SPDX-License-Identifier: GPL-3.0-only
# Thin in-memory application assembly. No function runs when this file is sourced.
# The statistical equations are supplied by application_statistics.R / wdsmatch.
wm_application_function_identity <- function(fn) {
  if (!is.function(fn)) stop("Runtime dependency must be a function")
  text <- function(x) paste(deparse(x, width.cutoff = 500L,
    control = c("keepNA", "keepInteger")), collapse = "\n")
  syntax <- if (is.primitive(fn)) list(primitive = text(fn)) else
    list(formals = text(formals(fn)), body = text(body(fn)))
  digest::digest(syntax, algo = "sha256")
}

wm_application_source_definitions <- function(path, expected_sha256, wanted) {
  path <- normalizePath(path, mustWork = TRUE)
  if (dir.exists(path) || !identical(digest::digest(file = path, algo = "sha256",
      serialize = FALSE), expected_sha256)) stop("Accepted application source bytes changed")
  definitions <- list()
  for (node in parse(file = path, keep.source = FALSE)) {
    if (is.call(node) && identical(node[[1L]], as.name("<-")) &&
        is.symbol(node[[2L]]) && as.character(node[[2L]]) %in% wanted &&
        is.call(node[[3L]]) && identical(node[[3L]][[1L]], as.name("function"))) {
      name <- as.character(node[[2L]])
      if (name %in% names(definitions)) stop("Duplicate accepted source definition")
      # Evaluate only a function expression to obtain syntax; no body is called.
      definitions[[name]] <- eval(node[[3L]], envir = baseenv())
    }
  }
  if (!setequal(names(definitions), wanted)) stop("Accepted source definitions missing")
  definitions[wanted]
}

wm_application_default_modules <- function(component_source, modules = NULL) {
  directory <- dirname(normalizePath(component_source, mustWork = TRUE))
  expected_solver <- wm_application_source_definitions(
    file.path(directory, "case_propensity_solver.R"),
    "35a7cad805577c79ac1e66078fc838959310d8e3b43f65eb2e95d4020278890e",
    "wdsm_case_fit_propensity")[[1L]]
  expected_provider <- wm_application_source_definitions(
    file.path(directory, "application_modules.R"),
    "23233ec83ff188bd800e43e96a6db01ee2900fe5934a33745bafd51a1dcb896a",
    "wm_application_modules")[[1L]]
  provider <- get("wm_application_modules", envir = environment(wm_application_run), inherits = TRUE)
  if (!identical(wm_application_function_identity(provider),
      wm_application_function_identity(expected_provider))) stop("Loaded default provider differs from accepted source")
  solver <- get("wdsm_case_fit_propensity", envir = environment(provider), inherits = TRUE)
  if (!identical(wm_application_function_identity(solver),
      wm_application_function_identity(expected_solver))) stop("Loaded default solver differs from accepted source")
  if (is.null(modules)) modules <- provider()
  if (!is.list(modules) || !identical(names(modules), c("source", "common")) ||
      !is.list(modules$source) || !identical(names(modules$source), "wdsm_case_fit_propensity") ||
      !identical(wm_application_function_identity(modules$source$wdsm_case_fit_propensity),
        wm_application_function_identity(expected_solver)) ||
      !identical(modules$common, asNamespace("wdsmatch")))
    stop("Default modules must use the accepted solver and actual wdsmatch namespace; custom modules require a custom callback/recipe/dependencies")
  modules
}

wm_application_package_identity <- function(ns) {
  if (!isNamespace(ns)) stop("Package implementation must be an actual namespace")
  package <- getNamespaceName(ns)
  if (!identical(ns, getNamespace(package))) stop("Supplied namespace identity changed")
  directory <- normalizePath(getNamespaceInfo(ns, "path"), mustWork = TRUE)
  files <- c("DESCRIPTION", "NAMESPACE")
  for (part in c("R", "libs")) {
    subdirectory <- file.path(directory, part)
    if (dir.exists(subdirectory)) {
      paths <- list.files(subdirectory, recursive = TRUE, all.files = TRUE, full.names = FALSE)
      paths <- paths[!file.info(file.path(subdirectory, paths))$isdir]
      files <- c(files, paste0(part, "/", paths))
    }
  }
  files <- sort(unique(files), method = "radix")
  if (anyNA(file.info(file.path(directory, files))$size) ||
      any(file.info(file.path(directory, files))$isdir)) stop("Installed package payload missing")
  payload <- setNames(vapply(files, function(f) digest::digest(file = file.path(directory, f),
    algo = "sha256", serialize = FALSE), character(1)), files)
  keys <- sort(ls(ns, all.names = TRUE), method = "radix")
  if (any(vapply(keys, bindingIsActive, logical(1), env = ns)))
    stop("Active namespace bindings require a separately declared custom runtime")
  functions <- keys[vapply(keys, function(k) is.function(get(k, envir = ns, inherits = FALSE)), logical(1))]
  signatures <- setNames(vapply(functions, function(k)
    wm_application_function_identity(get(k, envir = ns, inherits = FALSE)), character(1)), functions)
  list(package = package, loaded_version = as.character(getNamespaceVersion(ns)),
    installed_version = as.character(utils::packageVersion(package, lib.loc = dirname(directory))),
    payload_sha256 = as.list(payload), loaded_function_sha256 = as.list(signatures))
}

wm_application_runtime_identity <- function(modules, statistics_source, custom_dependencies = NULL,
                                            default = TRUE) {
  names_statistical <- c("app_capture", "app_count_columns", "app_ps", "app_components",
    "app_predictions", "app_row", "app_match", "app_replicate")
  expected <- wm_application_source_definitions(statistics_source,
    "a42b2fffc76bc7c6875c24031b032f88a33f0967f151494cb8925d274eab11cc", names_statistical)
  signatures <- vapply(expected, wm_application_function_identity, character(1))
  loaded <- lapply(names_statistical, get, envir = environment(wm_application_run), inherits = TRUE)
  if (!identical(vapply(loaded, wm_application_function_identity, character(1)), unname(signatures)))
    stop("Loaded application statistical definitions differ from accepted source")
  # Bind actual lexical helper resolution, without hashing closure environments.
  for (fn in loaded) for (name in names_statistical) {
    resolved <- get(name, envir = environment(fn), inherits = TRUE)
    if (!identical(wm_application_function_identity(resolved), signatures[[name]]))
      stop("Application helper lexical resolution differs from accepted source")
  }
  if (!is.list(modules) || !is.list(modules$source) ||
      !is.function(modules$source$wdsm_case_fit_propensity)) stop("Explicit source solver module required")
  required <- c(".wm_ns_ols", ".wm_ns_basis", "wm_match", "wm_bootstrap_refit", "wm_wdsm_fit")
  if (isNamespace(modules$common)) {
    if (!all(vapply(required, function(name) exists(name, envir = modules$common, inherits = FALSE) &&
        is.function(get(name, envir = modules$common, inherits = FALSE)), logical(1))))
      stop("Actual package namespace lacks required application functions")
    implementation <- wm_application_package_identity(modules$common)
  } else {
    if (default || !is.list(modules$common) || !setequal(names(modules$common), required) ||
        anyDuplicated(names(modules$common)) || !all(vapply(modules$common, is.function, logical(1))))
      stop("Custom common module must be an actual package namespace or an explicit required-function list")
    implementation <- list(custom_loaded_function_sha256 = as.list(
      vapply(modules$common[required], wm_application_function_identity, character(1))))
  }
  dependencies <- list()
  if (!default) {
    if (!is.list(custom_dependencies) || !length(custom_dependencies) ||
        is.null(names(custom_dependencies)) || anyNA(names(custom_dependencies)) ||
        anyDuplicated(names(custom_dependencies)) ||
        any(!grepl("^[A-Za-z0-9_.:-]+$", names(custom_dependencies))))
      stop("Custom callbacks/modules require explicit named source-pinned dependencies")
    for (key in sort(names(custom_dependencies), method = "radix")) {
      pin <- custom_dependencies[[key]]
      if (!is.list(pin) || !setequal(names(pin), c("source_file", "source_sha256")) ||
          anyDuplicated(names(pin)) || !is.character(pin$source_file) || length(pin$source_file) != 1L ||
          is.na(pin$source_file) || !is.character(pin$source_sha256) || length(pin$source_sha256) != 1L ||
          is.na(pin$source_sha256) || !grepl("^[0-9a-f]{64}$", pin$source_sha256)) stop("Malformed custom dependency source pin")
      path <- normalizePath(pin$source_file, mustWork = TRUE)
      if (dir.exists(path) || !identical(digest::digest(file = path, algo = "sha256", serialize = FALSE),
          pin$source_sha256)) stop("Custom runtime dependency source bytes changed")
      dependencies[[key]] <- pin$source_sha256
    }
  } else if (!is.null(custom_dependencies)) stop("Default dependency specification cannot be overridden")
  list(R_version = as.character(getRversion()), R_platform = R.version$platform,
    application_helper_sha256 = as.list(signatures),
    actual_solver_sha256 = wm_application_function_identity(modules$source$wdsm_case_fit_propensity),
    implementation = implementation, custom_dependency_source_sha256 = dependencies)
}

wm_application_recipe_binding <- function(component_callback, component_recipe,
                                          component_source, modules = NULL,
                                          application_statistics_source = NULL,
                                          component_dependencies = NULL) {
  # Source paths are used only to verify bytes and are never included in bindings.
  if (!is.character(component_source) || length(component_source) != 1L ||
      is.na(component_source) || !nzchar(component_source))
    stop("Supply the nuisance implementation source file explicitly")
  source_file <- normalizePath(component_source, mustWork = TRUE)
  if (dir.exists(source_file) || file.access(source_file, 4L) != 0L)
    stop("Nuisance implementation source must be a readable file")
  source_sha256 <- digest::digest(file = source_file, algo = "sha256", serialize = FALSE)
  if (is.null(component_callback)) {
    modules <- wm_application_default_modules(component_source, modules)
    if (!is.null(application_statistics_source)) stop("Default statistics source cannot be overridden")
    application_statistics_source <- component_source
    if (!is.null(component_recipe)) stop("Default nuisance recipe cannot be overridden")
    recipe <- list(recipe_id = "accepted_app_components_v1",
      source_sha256 = "a42b2fffc76bc7c6875c24031b032f88a33f0967f151494cb8925d274eab11cc",
      metadata = list(propensity_solver_sha256 =
        "35a7cad805577c79ac1e66078fc838959310d8e3b43f65eb2e95d4020278890e"))
  } else {
    if (is.null(modules) || is.null(application_statistics_source))
      stop("Custom callbacks require explicitly supplied modules and canonical application statistics source")
    if (!is.function(component_callback) || !is.list(component_recipe) ||
        !setequal(names(component_recipe), c("recipe_id", "source_sha256", "metadata")) ||
        anyDuplicated(names(component_recipe)))
      stop("Custom callbacks require explicit recipe_id/source_sha256/metadata")
    normalize_metadata <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.list(x) && !is.object(x)) {
        if (length(setdiff(names(attributes(x)), "names")))
          stop("Recipe metadata lists cannot carry non-name attributes")
        if (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
            any(!nzchar(names(x))) || anyDuplicated(names(x))))
          stop("Recipe metadata lists must have unique explicit names")
        if (!length(x)) return(list())
        return(lapply(x[order(names(x), method = "radix")], normalize_metadata))
      }
      if (!is.null(attributes(x)) || anyNA(x))
        stop("Recipe metadata cannot contain objects, closures, attributes or missing values")
      if (is.character(x)) {
        if (any(grepl("^(/|[A-Za-z]:[/\\\\]|file://)", x)))
          stop("Recipe metadata must not include local filesystem paths")
        return(enc2utf8(x))
      }
      if (is.logical(x)) return(x)
      if (is.numeric(x) && !is.complex(x) && all(is.finite(x))) return(as.double(x))
      stop("Recipe metadata must be named lists of plain finite numeric/logical/character vectors")
    }
    recipe <- list(recipe_id = component_recipe$recipe_id,
      source_sha256 = component_recipe$source_sha256,
      metadata = normalize_metadata(component_recipe$metadata))
    if (!is.character(recipe$recipe_id) || length(recipe$recipe_id) != 1L ||
        is.na(recipe$recipe_id) || !grepl("^[A-Za-z0-9_.:-]+$", recipe$recipe_id) ||
        !is.character(recipe$source_sha256) || length(recipe$source_sha256) != 1L ||
        is.na(recipe$source_sha256) || !grepl("^[0-9a-f]{64}$", recipe$source_sha256))
      stop("Invalid declared nuisance recipe/source identity")
    recipe$recipe_id <- enc2utf8(recipe$recipe_id)
  }
  if (!identical(source_sha256, recipe$source_sha256))
    stop("Declared nuisance implementation source bytes changed")
  recipe$runtime_dependencies <- wm_application_runtime_identity(modules,
    application_statistics_source, component_dependencies, default = is.null(component_callback))
  if (is.null(component_callback)) {
    recipe$runtime_dependencies$default_source_sha256 <- list(
      statistics = "a42b2fffc76bc7c6875c24031b032f88a33f0967f151494cb8925d274eab11cc",
      solver = "35a7cad805577c79ac1e66078fc838959310d8e3b43f65eb2e95d4020278890e",
      provider = "23233ec83ff188bd800e43e96a6db01ee2900fe5934a33745bafd51a1dcb896a")
    recipe$runtime_dependencies$actual_provider_sha256 <- wm_application_function_identity(
      get("wm_application_modules", envir = environment(wm_application_run), inherits = TRUE))
  }
  list(specification = recipe, sha256 = digest::digest(recipe, algo = "sha256"))
}

wm_application_count_input <- function(a, counts, count_ids, provenance) {
  if (!identical(count_ids, a$id)) stop("Count rows must bind exact prepared IDs and order")
  if (!is.matrix(counts) || !is.numeric(counts) || is.complex(counts) ||
      !identical(dim(counts), c(as.integer(a$n), 200L))) stop("Explicit numeric n-by-200 counts required")
  if (!is.list(provenance) || !length(provenance)) stop("Explicit supplied-count provenance required")
  validation <- app_count_columns(a, counts)
  list(matrix = counts, id = count_ids, B = 200L, validation = validation,
    sha256 = digest::digest(counts, algo = "sha256"), provenance = provenance,
    all_columns_valid = all(vapply(validation, `[[`, logical(1), "ok")),
    scope = "Supplied counts only; no count generation, redraw or reduced-B interval")
}

wm_application_shared_refits <- function(a, modules, count_input,
                                         component_callback = NULL) {
  if (is.null(component_callback)) component_callback <- function(a, modules, m)
    app_components(a, modules, m)
  if (!is.function(component_callback)) stop("component_callback must be a function(a, modules, m)")
  counts <- count_input$matrix
  stopifnot(identical(count_input$id, a$id),
    identical(count_input$sha256, digest::digest(counts, algo = "sha256")),
    identical(count_input$validation, app_count_columns(a, counts)))
  # These allocation, arm/component stripping and column assignments follow app_refits.
  means <- list(DSM = list(mean0 = matrix(NA_real_, a$n, 200L), mean1 = matrix(NA_real_, a$n, 200L)),
    linear = list(mean0 = matrix(NA_real_, a$n, 200L), mean1 = matrix(NA_real_, a$n, 200L)))
  diagnostics <- vector("list", 200L)
  for (b in seq_len(200L)) {
    captured <- app_capture(function() {
      if (!count_input$validation[[b]]$ok) stop(count_input$validation[[b]]$error)
      component_callback(a, modules, counts[, b])
    })
    if (!captured$ok) {diagnostics[[b]] <- captured; next}
    c <- captured$value
    column <- app_capture(function() {
      if (!is.list(c)) stop("Nuisance callback must return the captured component structure")
      if (length(captured$warnings)) c$callback_warnings <- captured$warnings
      predictions <- list(DSM = list(mean0 = rep(NA_real_, a$n), mean1 = rep(NA_real_, a$n)),
                          linear = list(mean0 = rep(NA_real_, a$n), mean1 = rep(NA_real_, a$n)))
      for (kind in c("DSM", "linear")) for (z in 0:1) {
        v <- c[[kind]][[as.character(z)]]
        if (!is.list(v) || !is.logical(v$ok) || length(v$ok) != 1L || is.na(v$ok))
          stop("Malformed captured nuisance component: ", kind, " arm", z)
        if (v$ok) {
          if (!is.numeric(v$value$mean) || is.complex(v$value$mean) ||
              !is.null(dim(v$value$mean)) || length(v$value$mean) != a$n ||
              any(!is.finite(v$value$mean))) stop("Invalid nuisance mean prediction")
          predictions[[kind]][[paste0("mean", z)]] <- v$value$mean
        }
        v$value <- if (v$ok) list(coefficients = v$value$coefficients) else NULL
        c[[kind]][[as.character(z)]] <- v
      }
      if (isTRUE(c$PS$ok)) c$PS$value$probability <- NULL
      if (isTRUE(c$PS_score$ok)) c$PS_score$value <- NULL
      for (z in 0:1) if (isTRUE(c$PG[[as.character(z)]]$ok)) c$PG[[as.character(z)]]$value$mean <- NULL
      list(diagnostics = c, predictions = predictions)
    })
    if (!column$ok) {diagnostics[[b]] <- column; next}
    for (kind in c("DSM", "linear")) for (z in 0:1)
      means[[kind]][[paste0("mean", z)]][, b] <- column$value$predictions[[kind]][[paste0("mean", z)]]
    diagnostics[[b]] <- column$value$diagnostics
  }
  list(means = means, diagnostics = diagnostics, B = 200L,
    counts_sha256 = count_input$sha256,
    scope = "Shared supplied-count predictions; all 200 success/failure columns retained")
}

wm_application_check_refits <- function(a, count_input, refits) {
  if (!is.list(refits) || !is.numeric(refits$B) || length(refits$B) != 1L ||
      is.na(refits$B) || refits$B != 200L ||
      !identical(refits$counts_sha256, count_input$sha256) ||
      !is.list(refits$diagnostics) || length(refits$diagnostics) != 200L)
    stop("Shared refits must retain all 200 diagnostics and bind the supplied count matrix")
  for (kind in c("DSM", "linear")) for (arm in c("mean0", "mean1")) {
    value <- refits$means[[kind]][[arm]]
    if (!is.matrix(value) || !is.numeric(value) || is.complex(value) ||
        !identical(dim(value), c(as.integer(a$n), 200L))) stop("Malformed shared prediction matrix")
    # Nonfinite columns remain failure evidence; app_replicate checks required arms.
  }
  invisible(TRUE)
}

wm_application_check_point <- function(a, components, point, family, M, estimand) {
  if (!is.list(point) || !is.logical(point$ok) || length(point$ok) != 1L || is.na(point$ok))
    stop("Supplied point must retain its captured success/failure record")
  if (!point$ok) {
    if (!is.character(point$error) || length(point$error) != 1L ||
        is.na(point$error) || !nzchar(point$error)) stop("Missing point failure reason")
    return(invisible(TRUE))
  }
  p <- point$value
  if (!inherits(p, "wm_match") || !identical(p$method, "self_normalized") ||
      p$n != a$n || p$M != M || !identical(p$estimand, estimand) ||
      !identical(as.numeric(p$data$Y), as.numeric(a$Y)) ||
      !identical(as.numeric(p$data$Z), as.numeric(a$Z)) ||
      !identical(as.numeric(p$weights), as.numeric(a$W)) ||
      !identical(as.numeric(p$analysis_weights), as.numeric(a$W/p$weight_scale)) ||
      !isTRUE(p$info$corrected) || isTRUE(p$info$nuisance_correction) ||
      !identical(p$graph$fold_id, rep.int("all", a$n)) ||
      !identical(p$graph$strata, rep.int("all", a$n)) ||
      !identical(p$graph$tie_diagnostics$rule, "source_random") ||
      !identical(p$graph$tie_diagnostics$seed, as.integer(20260917L + M)) ||
      !identical(p$graph$tie_diagnostics$tolerance, 64 * .Machine$double.eps))
    stop("Supplied point data, method, M, target, W or tie convention differs")
  predictions <- app_predictions(components, family, estimand)
  scores0 <- switch(family, PS = components$PS_score$value,
    DSM = components$DSM[["0"]]$value$scores, full_X = components$full_X_scores$value)
  scores1 <- if (estimand == "PATE") switch(family, PS = scores0,
    DSM = components$DSM[["1"]]$value$scores, full_X = scores0) else NULL
  canonical_score <- function(x, n, name) {
    if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
        nrow(x) != n || ncol(x) < 1L || anyNA(x) || any(!is.finite(x))) {
      stop(name, " must be a finite numeric matrix with ", n,
           " rows and at least one column.", call. = FALSE)
    }
    storage.mode(x) <- "double"
    unname(x)
  }
  scores0 <- canonical_score(scores0, a$n, "scores0")
  if (!is.null(scores1)) scores1 <- canonical_score(scores1, a$n, "scores1")
  if (!identical(p$graph$scores0, scores0) || !identical(p$graph$scores1, scores1) ||
      !identical(as.numeric(p$predictions$mean0), as.numeric(predictions$mean0)) ||
      (estimand == "PATE" && !identical(as.numeric(p$predictions$mean1), as.numeric(predictions$mean1))))
    stop("Supplied point maps or predictions differ from the bound components")
  invisible(TRUE)
}

wm_application_run <- function(a, counts, count_ids, count_provenance,
                               modules = NULL, components = NULL, refits = NULL,
                               precomputed_points = NULL, reuse_binding = NULL,
                               component_callback = NULL, component_recipe = NULL,
                               component_source, application_statistics_source = NULL,
                               component_dependencies = NULL, target_scope,
                               inference_scope) {
  # The caller supplies all prepared input objects; this driver never selects raw rows.
  verified <- wm_application_input(a$data, a$id, a$Y, a$Z, a$W, a$X, a$design,
    a$variables, a$study, a$outcome, a$baseline, a$units, a$ps_weighting, a$provenance)
  if (!identical(verified$sample_sha256, a$sample_sha256) ||
      !identical(verified$dataset_sha256, a$dataset_sha256) ||
      !identical(verified$nuisance_sha256, a$nuisance_sha256))
    stop("Prepared input hash/order or nuisance specification changed")
  if (!is.character(target_scope) || !identical(names(target_scope), c("PATE", "PATT")) ||
      anyNA(target_scope) || any(!nzchar(target_scope)) || length(inference_scope) != 1L ||
      !is.character(inference_scope) || is.na(inference_scope) || !nzchar(inference_scope))
    stop("Declare target-specific scope and an inference scope; no causal validity is inferred")
  if (is.null(component_callback)) modules <- wm_application_default_modules(component_source, modules)
  else if (is.null(modules)) stop("Custom callbacks require explicitly supplied modules")
  nuisance_recipe <- wm_application_recipe_binding(component_callback, component_recipe, component_source,
    modules, application_statistics_source, component_dependencies)
  count_input <- wm_application_count_input(a, counts, count_ids, count_provenance)
  callback_declared <- component_callback
  uses_reuse <- !is.null(components) || !is.null(refits) || !is.null(precomputed_points)
  if (uses_reuse) {
    if (!is.list(reuse_binding) || !identical(reuse_binding$sample_sha256, a$sample_sha256) ||
        !identical(reuse_binding$dataset_sha256, a$dataset_sha256) ||
        !identical(reuse_binding$counts_sha256, count_input$sha256) ||
        !identical(reuse_binding$nuisance_sha256, a$nuisance_sha256) ||
        !identical(reuse_binding$component_recipe_sha256, nuisance_recipe$sha256))
      stop("Explicit accepted input/count/nuisance/recipe reuse binding required")
    check_hash <- function(value, hash, label) {
      if (!is.character(hash) || length(hash) != 1L ||
          !identical(digest::digest(value, algo = "sha256"), hash)) stop("Reuse object hash differs: ", label)
    }
    if (!is.null(components)) check_hash(components, reuse_binding$components_sha256, "point components")
    if (!is.null(refits)) check_hash(refits, reuse_binding$refits_sha256, "shared refits")
    if (!is.null(precomputed_points)) check_hash(precomputed_points, reuse_binding$points_sha256, "point records")
  }
  if (is.null(component_callback)) component_callback <- function(a, modules, m)
    app_components(a, modules, m)
  if (!is.function(component_callback)) stop("Explicit nuisance component callback required")
  if (is.null(components)) components <- component_callback(a, modules, rep(1, a$n))
  if (is.null(refits)) refits <- wm_application_shared_refits(a, modules, count_input, component_callback)
  wm_application_check_refits(a, count_input, refits)
  requested <- expand.grid(M = c(1L, 3L, 5L), family = c("PS", "DSM", "full_X"),
    estimand = c("PATE", "PATT"), stringsAsFactors = FALSE)
  keys <- paste0(requested$family, "_M", requested$M, "_", requested$estimand)
  if (!is.null(precomputed_points) && (!is.list(precomputed_points) ||
      !setequal(names(precomputed_points), keys) || anyDuplicated(names(precomputed_points))))
    stop("Reuse must supply exactly all 18 captured point records; no implicit rematching")
  rows <- points <- replications <- list()
  for (i in seq_len(nrow(requested))) {
    family <- requested$family[i]; M <- requested$M[i]; estimand <- requested$estimand[i]; key <- keys[i]
    point <- if (!is.null(precomputed_points)) precomputed_points[[key]] else
      app_capture(function() app_match(a, modules, components, family, M, estimand))
    wm_application_check_point(a, components, point, family, M, estimand)
    replication <- app_capture(function() {if (!point$ok) stop(point$error)
      app_replicate(point$value, modules, counts, refits, family, estimand)})
    row <- app_row(a, family, M, estimand, point, replication)
    row$target_scope <- target_scope[[estimand]]; row$inference_scope <- inference_scope
    rows[[key]] <- row; points[[key]] <- point; replications[[key]] <- replication
  }
  if (!identical(nuisance_recipe,
      wm_application_recipe_binding(callback_declared, component_recipe, component_source,
        modules, application_statistics_source, component_dependencies)))
    stop("Declared nuisance implementation source/recipe changed during execution")
  list(rows = do.call(rbind, rows), requested = requested, points = points,
    replications = replications, components = components, refits = refits,
    count_input = count_input, sample = a, reuse_binding = reuse_binding,
    nuisance_recipe = nuisance_recipe,
    status = "APPLICATION_ROWS_AND_FAILURES_RETAINED_REVIEW_REQUIRED",
    scientific_acceptance = FALSE)
}
