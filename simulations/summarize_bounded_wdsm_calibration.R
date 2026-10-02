# Read-only aggregation of the frozen bounded-WDSM calibration records.
# No fits, RNG calls, experiments, or package loading occur on source().

.wm_ca_key <- function(x, fields = c("dataset_id", "branch", "model", "estimand", "M"))
  do.call(paste, c(x[fields], sep = "|"))
.wm_ca_mean <- function(x) if (length(x)) mean(x) else NA_real_
.wm_ca_mcse <- function(x) if (length(x) > 1L) stats::sd(x)/sqrt(length(x)) else NA_real_
.wm_ca_var <- function(x) if (length(x) > 1L) stats::var(x) else NA_real_
.wm_ca_bind <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
.wm_ca_check <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

# This computes conditional count-law diagnostics, not outer sampling accuracy.
wm_cal_aggregate_counts <- function(draws, B, estimate, analytic_V, n, truth) {
  .wm_ca_check(B >= 2L && B == as.integer(B), "Invalid prescribed B.")
  .wm_ca_check(all(c("draw","status","estimate","one_step_root") %in% names(draws)),
               "Missing draw columns.")
  .wm_ca_check(!anyDuplicated(draws$draw) && all(draws$draw %in% seq_len(B)),
               "Duplicate or unrequested count draw.")
  full <- data.frame(draw = seq_len(B), status = "missing", estimate = NA_real_,
                     one_step_root = NA_real_, error = "Missing prescribed draw")
  j <- match(draws$draw, full$draw)
  for (field in intersect(names(full), names(draws))) full[j,field] <- draws[[field]]
  .wm_ca_check(all(full$status %in% c("completed","failed","pending","missing","terminal_unproduced")),
               "Unknown count status.")
  .wm_ca_check(all(is.finite(full$estimate[full$status == "completed"])),
               "Completed count draw has a nonfinite estimate.")
  good <- full$status == "completed"
  refit_ok <- all(good) && is.finite(estimate)
  one_ok <- all(is.finite(full$one_step_root))
  if (is.finite(estimate) && "root_n" %in% names(draws)) {
    complete_rows <- draws$status=="completed"
    expected <- sqrt(n)*(draws$estimate[complete_rows]-estimate)
    .wm_ca_check(all(is.finite(draws$root_n[complete_rows])) &&
      all(abs(draws$root_n[complete_rows]-expected)<=1e-9*(1+abs(expected))),
      "Saved count roots disagree with point-to-root-n conversion.")
  }
  out <- list(B = B, successful_draws = sum(good), missing_draws = sum(full$status == "missing"),
    terminal_unproduced_draws = sum(full$status == "terminal_unproduced"),
    all_refits_succeeded = refit_ok, all_one_step_available = one_ok,
    source_point_variance_B = NA_real_, root_n_variance = NA_real_,
    interval_status = "unavailable", se = NA_real_, lower = NA_real_, upper = NA_real_,
    covered = NA, one_step_exact_root_variance = analytic_V,
    one_step_root_variance_Bminus1 = NA_real_, refit_root_variance_mcse = NA_real_,
    one_step_root_variance_mcse = NA_real_, paired_root_variance_difference = NA_real_,
    paired_variance_difference_mcse = NA_real_, root_difference_mean = NA_real_,
    root_difference_mean_mcse = NA_real_, root_difference_rms = NA_real_,
    root_difference_rms_mcse = NA_real_)
  if (one_ok) {
    qo <- (full$one_step_root-mean(full$one_step_root))^2
    out$one_step_root_variance_Bminus1 <- stats::var(full$one_step_root)
    out$one_step_root_variance_mcse <- B/(B-1)*.wm_ca_mcse(qo)
  }
  if (refit_ok) {
    roots <- sqrt(n)*(full$estimate-estimate)
    qr <- (roots-mean(roots))^2
    out$source_point_variance_B <- mean((full$estimate-mean(full$estimate))^2)
    out$root_n_variance <- n*B/(B-1)*out$source_point_variance_B
    out$refit_root_variance_mcse <- B/(B-1)*.wm_ca_mcse(qr)
    if (out$source_point_variance_B > 0) {
      out$interval_status <- "completed"; out$se <- sqrt(out$source_point_variance_B)
      out$lower <- estimate-stats::qnorm(.975)*out$se
      out$upper <- estimate+stats::qnorm(.975)*out$se
      out$covered <- truth >= out$lower && truth <= out$upper
    }
    if (one_ok) {
      delta <- roots-full$one_step_root
      out$paired_root_variance_difference <- stats::var(roots)-stats::var(full$one_step_root)
      out$paired_variance_difference_mcse <- B/(B-1)*.wm_ca_mcse(qr-qo)
      out$root_difference_mean <- mean(delta)
      out$root_difference_mean_mcse <- .wm_ca_mcse(delta)
      out$root_difference_rms <- sqrt(mean(delta^2))
      out$root_difference_rms_mcse <- if (out$root_difference_rms > 0)
        .wm_ca_mcse(delta^2)/(2*out$root_difference_rms) else 0
    }
  }
  list(summary = as.data.frame(out, stringsAsFactors = FALSE), draws = full)
}

# Full requested grid is authoritative; a missing output is retained as missing.
# For shards, pass every available records root. Missing shard outputs remain
# missing on the common global grid, and never shrink denominators.
wm_cal_collect_records <- function(records_roots) {
  roots <- normalizePath(records_roots, mustWork = TRUE)
  .wm_ca_check(!anyDuplicated(roots) && length(roots) > 0L, "Duplicate or empty records roots.")
  configs <- lapply(roots, function(r) jsonlite::read_json(file.path(r,"config.json"), simplifyVector=TRUE))
  provenance <- lapply(roots,function(r) jsonlite::read_json(file.path(r,"provenance.json"),simplifyVector=TRUE))
  for(i in seq_along(roots)) {
    .wm_ca_check(identical(digest::digest(file=file.path(roots[i],"config.json"),algo="sha256",serialize=FALSE),
                          provenance[[i]]$config_sha256),"Saved configuration hash disagrees with provenance.")
    for(field in c("source_sha256","package_sha256","package_source_sha256","package_version","mode"))
      .wm_ca_check(!is.null(configs[[i]][[field]]) && isTRUE(all.equal(configs[[i]][[field]],provenance[[i]][[field]])),
                    paste("Configuration/provenance disagreement:",field))
  }
  shard_counts <- vapply(provenance,function(p)if(is.null(p$shard_count))1 else p$shard_count,numeric(1L))
  .wm_ca_check(length(unique(shard_counts))==1L,"Cannot combine different shard allocation schemes.")
  mode <- unique(vapply(configs, function(c) c$mode, character(1L)))
  .wm_ca_check(length(mode) == 1L && mode %in% c("runtime_pilot","calibration"), "Mixed or unknown modes.")
  scientific <- c("index","radius","sample_sizes","models","estimands","M","datasets_per_n",
    "count_draws","master_seed","refit_modulus","refit_remainder","reviewed_plan_sha256",
    "source_sha256","package_source_sha256")
  for (cfg in configs) .wm_ca_check(isTRUE(all.equal(cfg[scientific],configs[[1L]][scientific])),
                                   "Scientific configuration or source pins differ.")
  required <- c("dataset_id","n","replication_id","requested_count_draws","branch","model","estimand","M")
  globals <- lapply(roots, function(r) {
    f <- file.path(r,"global_requested_grid.csv")
    if (!file.exists(f)) f <- file.path(r,"requested_grid.csv")
    x <- utils::read.csv(f, stringsAsFactors=FALSE)
    .wm_ca_check(all(required %in% names(x)), "Requested grid lacks required columns.")
    x[required]
  })
  grid <- globals[[1L]]; keys <- .wm_ca_key(grid)
  .wm_ca_check(!anyDuplicated(keys), "Duplicate global requested record.")
  cfg <- configs[[1L]]
  pilot <- mode=="runtime_pilot"
  .wm_ca_check(identical(cfg$index,"source_good") && cfg$radius==.25 &&
    identical(as.numeric(cfg$sample_sizes),c(500,2000)) &&
    identical(as.character(cfg$models),c("CorCor","MisCor")) &&
    identical(as.character(cfg$estimands),c("PATE","PATT")) && identical(as.numeric(cfg$M),c(1,3)) &&
    cfg$datasets_per_n==(if(pilot)4 else 1000) && cfg$count_draws==(if(pilot)39 else 399) &&
    cfg$master_seed==(if(pilot)29092026 else 29092027) && cfg$refit_modulus==(if(pilot)1 else 4) &&
    cfg$refit_remainder==1 && cfg$max_candidate_multiplier==16,
    "Configuration does not match the frozen calibration design.")
  expected_grid <- .wm_ca_bind(lapply(cfg$sample_sizes,function(n)
    .wm_ca_bind(lapply(seq_len(cfg$datasets_per_n),function(id) {
      cells <- rbind(expand.grid(branch=c("actual_fitted","oracle_means_fitted_map"),
        model=c("CorCor","MisCor"),estimand=c("PATE","PATT"),M=c(1,3),stringsAsFactors=FALSE),
        expand.grid(branch="known_map_oracle",model="shared",estimand=c("PATE","PATT"),M=c(1,3),stringsAsFactors=FALSE))
      cbind(dataset_id=sprintf("%s-n%d-id%04d",mode,n,id),n=n,replication_id=id,
        requested_count_draws=if(mode=="runtime_pilot")39 else if((id-1L)%%4L==0L)399 else 0,cells)
    }))))
  j <- match(.wm_ca_key(expected_grid),keys)
  .wm_ca_check(nrow(grid)==nrow(expected_grid) && !anyNA(j) &&
    isTRUE(all.equal(expected_grid[required],grid[j,required],check.attributes=FALSE)),
    "Global requested grid differs from the full predeclared design.")
  for (x in globals) {
    j <- match(keys,.wm_ca_key(x))
    .wm_ca_check(nrow(x)==nrow(grid) && !anyNA(j) &&
      isTRUE(all.equal(grid,x[j,,drop=FALSE],check.attributes=FALSE)), "Global grids differ.")
  }
  owners <- .wm_ca_bind(lapply(seq_along(roots),function(i) {
    x <- utils::read.csv(file.path(roots[i],"requested_grid.csv"),stringsAsFactors=FALSE)
    .wm_ca_check(all(.wm_ca_key(x) %in% keys),"Shard requests outside global grid.")
    wanted <- grid[grid$dataset_id %in% x$dataset_id,,drop=FALSE]
    .wm_ca_check(nrow(x)==nrow(wanted) && !anyDuplicated(.wm_ca_key(x)) &&
      setequal(.wm_ca_key(x),.wm_ca_key(wanted)),"Shard omits requested records for an owned dataset.")
    shard_count <- if(is.null(provenance[[i]]$shard_count))1L else provenance[[i]]$shard_count
    shard_index <- if(is.null(provenance[[i]]$shard_index))1L else provenance[[i]]$shard_index
    .wm_ca_check(shard_count>=1 && shard_count==as.integer(shard_count) &&
      shard_index>=1 && shard_index<=shard_count && shard_index==as.integer(shard_index),"Invalid shard ownership metadata.")
    assigned <- unique(grid$dataset_id[1+floor((grid$replication_id-1)*shard_count/cfg$datasets_per_n)==shard_index])
    .wm_ca_check(setequal(unique(x$dataset_id),assigned),"Requested shard ownership disagrees with fixed allocation.")
    data.frame(dataset_id=unique(x$dataset_id),root=rep(roots[i],length(unique(x$dataset_id))),stringsAsFactors=FALSE)
  }))
  .wm_ca_check(!anyDuplicated(owners$dataset_id), "Multiple roots claim one dataset.")
  streams <- lapply(roots,function(r)readRDS(file.path(r,"streams.rds")))
  for(s in streams) .wm_ca_check(isTRUE(all.equal(s,streams[[1L]])),"Global stream allocations differ.")
  allocation <- streams[[1L]]
  allocation_ids <- vapply(allocation,function(a)a$dataset_id,character(1L))
  .wm_ca_check(!anyDuplicated(allocation_ids) && setequal(allocation_ids,unique(grid$dataset_id)),
               "Global stream allocation IDs differ from requested datasets.")
  for(a in allocation) {
    g <- grid[grid$dataset_id==a$dataset_id,,drop=FALSE]
    .wm_ca_check(a$n==g$n[1L] && a$replication_id==g$replication_id[1L] &&
      identical(a$refit,g$requested_count_draws[1L]>0) &&
      length(a$data_seed)==7L && length(a$count_seed)==7L &&
      all(is.finite(a$data_seed)) && all(is.finite(a$count_seed)),"Invalid stream dataset/n/refit allocation.")
  }
  targets <- lapply(roots,function(r) readRDS(file.path(r,"target.rds")))
  .wm_ca_check(all(vapply(targets,function(t) isTRUE(all.equal(t,targets[[1L]])),logical(1L))),
               "Target definitions differ.")
  truth <- targets[[1L]]
  execution <- data.frame(dataset_id=unique(grid$dataset_id),status="missing_or_active",terminal_error="",stringsAsFactors=FALSE)
  states <- setNames(vector("list",nrow(execution)),execution$dataset_id)
  terminal_points <- terminal_draws <- states
  dataset_files <- states
  for (j in seq_len(nrow(execution))) {
    owner <- owners$root[match(execution$dataset_id[j],owners$dataset_id)]
    if(!is.na(owner)) dataset_files[[execution$dataset_id[j]]] <- list.files(file.path(owner,execution$dataset_id[j]))
    f <- if(is.na(owner)) NA_character_ else file.path(owner,execution$dataset_id[j],"completion.json")
    if(!is.na(f) && file.exists(f)) {
      state <- jsonlite::read_json(f,simplifyVector=TRUE)
      .wm_ca_check(state$status %in% c("completed","failed"),"Unknown dataset completion status.")
      .wm_ca_check(identical(state$dataset_id,execution$dataset_id[j]),"Dataset completion ID disagrees with directory/allocation.")
      .wm_ca_check(is.character(state$error) && length(state$error)==1L &&
        ((state$status=="completed" && state$error=="") || (state$status=="failed" && nzchar(state$error))),
        "Invalid dataset terminal error record.")
      ap <- file.path(owner,execution$dataset_id[j],"allocation.rds")
      if(file.exists(ap)) .wm_ca_check(isTRUE(all.equal(readRDS(ap),allocation[[match(state$dataset_id,allocation_ids)]])),
                                      "Dataset allocation disagrees with global streams.")
      execution$status[j] <- state$status
      execution$terminal_error[j] <- state$error; states[[state$dataset_id]] <- state
      if(state$status=="failed") {
        directory <- dirname(f); g <- grid[grid$dataset_id==state$dataset_id,,drop=FALSE]
        file <- file.path(directory,"terminal_record_accounting.csv")
        .wm_ca_check(file.exists(file),"Terminal failure lacks requested-record accounting.")
        acct <- utils::read.csv(file,stringsAsFactors=FALSE)
        .wm_ca_check(all(c(required,"record_file","produced","terminal_error") %in% names(acct)) &&
          nrow(acct)==nrow(g) && !anyDuplicated(.wm_ca_key(acct)),"Malformed terminal point accounting.")
        index <- match(.wm_ca_key(g),.wm_ca_key(acct))
        .wm_ca_check(!anyNA(index) && isTRUE(all.equal(g[required],acct[index,required],check.attributes=FALSE)),
                      "Terminal point accounting differs from requested records.")
        prefix <- ifelse(acct$branch=="known_map_oracle","known",acct$model)
        suffix <- ifelse(acct$branch=="known_map_oracle","-point.csv",
          ifelse(acct$branch=="actual_fitted","-actual-point.csv","-oracle-fitted-point.csv"))
        wanted_files <- paste0(prefix,"-",acct$estimand,"-",acct$M,suffix)
        .wm_ca_check(identical(acct$record_file,wanted_files) && !anyNA(acct$produced) &&
          all(acct$produced==file.exists(file.path(directory,wanted_files))) &&
          all(acct$terminal_error==state$error),"Terminal point production/error accounting is inconsistent.")
        terminal_points[[state$dataset_id]] <- acct
        B <- g$requested_count_draws[1L]
        if(B>0L) {
          file <- file.path(directory,"terminal_draw_accounting.csv")
          .wm_ca_check(file.exists(file),"Terminal failure lacks prescribed-draw accounting.")
          acct <- utils::read.csv(file,stringsAsFactors=FALSE)
          expected <- expand.grid(model=c("CorCor","MisCor"),estimand=c("PATE","PATT"),M=c(1,3),draw=seq_len(B),stringsAsFactors=FALSE)
          df <- c("model","estimand","M","draw")
          .wm_ca_check(all(c(df,"record_file","produced","terminal_error") %in% names(acct)) &&
            nrow(acct)==nrow(expected) && !anyDuplicated(.wm_ca_key(acct,df)) &&
            setequal(.wm_ca_key(acct,df),.wm_ca_key(expected,df)) && !anyNA(acct$produced) &&
            all(acct$terminal_error==state$error) &&
            identical(acct$record_file,sprintf("%s-%s-%d-count%03d.csv",acct$model,acct$estimand,acct$M,acct$draw)),
            "Malformed terminal draw production/error accounting.")
          terminal_draws[[state$dataset_id]] <- acct
        }
      }
    }
  }
  shard_complete <- logical(length(roots))
  for(i in seq_along(roots)) {
    file <- file.path(roots[i],"completion.json")
    if(!file.exists(file)) next
    top <- jsonlite::read_json(file,simplifyVector=FALSE)
    own <- owners$dataset_id[owners$root==roots[i]]; inventory <- top$dataset_status
    ids <- vapply(inventory,function(s)s$dataset_id,character(1L))
    .wm_ca_check(top$status %in% c("requested_grid_processed","empty_shard_processed") &&
      identical(top$mode,mode) && top$datasets_requested==length(own) &&
      !anyDuplicated(ids) && setequal(ids,own),"Shard completion inventory disagrees with owned allocation.")
    if(!is.null(top$global_datasets_requested)) .wm_ca_check(top$global_datasets_requested==nrow(execution),"Shard global dataset count mismatch.")
    for(field in c("shard_index","shard_count")) if(!is.null(provenance[[i]][[field]]))
      .wm_ca_check(identical(top[[field]],provenance[[i]][[field]]),"Shard completion index/count mismatch.")
    for(s in inventory) .wm_ca_check(!is.null(states[[s$dataset_id]]) &&
      identical(s$status,states[[s$dataset_id]]$status) && identical(s$error,states[[s$dataset_id]]$error),
      "Shard inventory disagrees with dataset completion.")
    .wm_ca_check(top$datasets_terminal_errors==sum(vapply(inventory,function(s)s$status=="failed",logical(1L))),
                  "Shard terminal-error count mismatch.")
    shard_complete[i] <- TRUE
  }
  rows <- list(); counts <- list(); draw_status <- list()
  for (i in seq_len(nrow(grid))) {
    g <- grid[i,,drop=FALSE]; owner <- owners$root[match(g$dataset_id,owners$dataset_id)]
    directory <- if (is.na(owner)) NA_character_ else file.path(owner,g$dataset_id)
    tag <- paste(if (g$branch=="known_map_oracle") "known" else g$model,g$estimand,g$M,sep="-")
    suffix <- if (g$branch=="oracle_means_fitted_map") "-oracle-fitted-point.csv" else
      if (g$branch=="actual_fitted") "-actual-point.csv" else "-point.csv"
    path <- if (is.na(directory)) NA_character_ else file.path(directory,paste0(tag,suffix))
    state <- states[[g$dataset_id]]
    if(!is.null(state) && state$status=="completed")
      .wm_ca_check(!is.na(path) && file.exists(path),"Completed dataset is missing a requested point CSV.")
    wanted_truth <- truth[[g$estimand]]
    .wm_ca_check(length(wanted_truth)==1L && is.finite(wanted_truth),"Invalid stored target.")
    row <- data.frame(branch=g$branch,model=g$model,estimand=g$estimand,M=g$M,n=g$n,
      truth=wanted_truth,point_status="missing",interval_status="unavailable",
      interval_unavailable_reason="Missing requested point record",error="Missing requested point record",
      estimate=NA_real_,root_n_variance=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_,
      covered=NA,V0=NA_real_,cross_term=NA_real_,nuisance_variance=NA_real_,
      naive_se=NA_real_,naive_lower=NA_real_,naive_upper=NA_real_,stringsAsFactors=FALSE)
    if(!is.null(state) && state$status=="failed" && (is.na(path)||!file.exists(path))) {
      row$point_status <- "terminal_unproduced"; row$error <- state$error
      row$interval_unavailable_reason <- paste("Terminal dataset error before point record:",state$error)
    }
    if (!is.na(path) && file.exists(path)) {
      x <- utils::read.csv(path,stringsAsFactors=FALSE)
      .wm_ca_check(nrow(x)==1L && all(names(row) %in% names(x)),"Malformed point record.")
      row <- x[names(row)]
      .wm_ca_check(all(row[c("branch","model","estimand","M","n")]==g[c("branch","model","estimand","M","n")]) &&
        abs(row$truth-wanted_truth)<1e-12,"Point record disagrees with requested key/target.")
      .wm_ca_check(row$point_status %in% c("completed","failed") &&
        row$interval_status %in% c("completed","unavailable"),"Unknown point/interval status.")
      if (row$point_status=="completed") .wm_ca_check(is.finite(row$estimate),"Nonfinite completed point.")
      if (row$interval_status=="completed") {
        .wm_ca_check(row$point_status=="completed" && is.finite(row$root_n_variance) && row$root_n_variance>0 &&
          is.finite(row$se) && row$se>0 && abs(row$se-sqrt(row$root_n_variance/g$n))<1e-10 &&
          abs(row$lower-(row$estimate-stats::qnorm(.975)*row$se))<1e-10 &&
          abs(row$upper-(row$estimate+stats::qnorm(.975)*row$se))<1e-10 &&
          identical(as.logical(row$covered),row$truth>=row$lower && row$truth<=row$upper),
          "Completed analytic interval fails scale/endpoint/coverage audit.")
      }
    }
    row$dataset_id <- g$dataset_id; row$replication_id <- g$replication_id
    row$requested_count_draws <- g$requested_count_draws; row$inner_variance_mcse <- NA_real_
    rows[[length(rows)+1L]] <- row
    if (g$branch=="actual_fitted") {
      naive <- row; naive$branch <- "fitted_naive"; naive$root_n_variance <- row$V0
      okay <- row$point_status=="completed" && is.finite(row$V0) && row$V0>0
      naive$interval_status <- if (okay) "completed" else "unavailable"
      naive$interval_unavailable_reason <- if (okay) "" else "Point or positive V0 unavailable"
      naive$se <- if (okay) sqrt(row$V0/g$n) else NA_real_
      naive$lower <- naive$estimate-stats::qnorm(.975)*naive$se
      naive$upper <- naive$estimate+stats::qnorm(.975)*naive$se
      naive$covered <- if (okay) wanted_truth>=naive$lower && wanted_truth<=naive$upper else NA
      rows[[length(rows)+1L]] <- naive
      if (g$requested_count_draws>0L) {
        B <- g$requested_count_draws
        empty <- data.frame(draw=integer(),status=character(),estimate=numeric(),one_step_root=numeric(),error=character())
        combined <- if (is.na(directory)) NA_character_ else file.path(directory,paste0(tag,"-draws.csv"))
        if(!is.null(state) && state$status=="completed")
          .wm_ca_check(!is.na(combined) && file.exists(combined),"Completed dataset lacks a final prescribed-draw table.")
        if (!is.na(combined) && file.exists(combined)) d <- utils::read.csv(combined,stringsAsFactors=FALSE) else {
          pieces <- lapply(seq_len(B),function(b) {
            f <- if (is.na(directory)) NA_character_ else file.path(directory,sprintf("%s-count%03d.csv",tag,b))
            if (!is.na(f) && file.exists(f)) utils::read.csv(f,stringsAsFactors=FALSE) else NULL
          }); pieces <- Filter(Negate(is.null),pieces)
          d <- if (length(pieces)) do.call(rbind,pieces) else empty
        }
        if(!is.na(combined) && file.exists(combined)) {
          files <- dataset_files[[g$dataset_id]]
          files <- files[startsWith(files,paste0(tag,"-count")) & grepl("-count[0-9]+[.]csv$",files)]
          for(file in files) {
            checkpoint <- utils::read.csv(file.path(directory,file),stringsAsFactors=FALSE)
            k <- match(checkpoint$draw,d$draw)
            .wm_ca_check(nrow(checkpoint)==1L && !anyNA(k) &&
              all(names(checkpoint) %in% names(d)) &&
              isTRUE(all.equal(checkpoint,d[k,names(checkpoint),drop=FALSE],check.attributes=FALSE)),
              "Final draw table disagrees with a retained per-draw checkpoint.")
          }
        }
        if(!is.null(state) && state$status=="completed") .wm_ca_check(nrow(d)==B &&
          !anyDuplicated(d$draw) && setequal(d$draw,seq_len(B)) && all(d$status %in% c("completed","failed")),
          "Completed dataset has missing/duplicate/nonterminal prescribed draw IDs.")
        if(!is.null(state) && state$status=="failed") {
          acct <- terminal_draws[[g$dataset_id]]
          acct <- acct[acct$model==g$model & acct$estimand==g$estimand & acct$M==g$M,,drop=FALSE]
          .wm_ca_check(all(acct$produced==(acct$draw %in% d$draw)),"Terminal draw accounting disagrees with available draw rows.")
          absent <- acct[!acct$produced,,drop=FALSE]
          if(nrow(absent)) {
            pad <- if(nrow(d))d[rep(NA_integer_,nrow(absent)),,drop=FALSE] else
              data.frame(draw=absent$draw,status="terminal_unproduced",estimate=NA_real_,one_step_root=NA_real_,error=state$error)
            pad$draw <- absent$draw;pad$status <- "terminal_unproduced";pad$estimate <- NA_real_
            pad$one_step_root <- NA_real_;pad$error <- state$error
            d <- rbind(d,pad)
          }
        }
        out <- wm_cal_aggregate_counts(d,B,row$estimate,row$root_n_variance,g$n,wanted_truth)
        id <- g[c("dataset_id","n","replication_id","model","estimand","M")]
        counts[[length(counts)+1L]] <- cbind(id,out$summary)
        draw_status[[length(draw_status)+1L]] <- cbind(id[rep(1L,B),,drop=FALSE],out$draws)
        refit <- row; refit$branch <- "original_refit"
        for (field in c("interval_status","se","lower","upper","covered","root_n_variance"))
          refit[[field]] <- out$summary[[field]]
        refit$inner_variance_mcse <- out$summary$refit_root_variance_mcse
        refit$interval_unavailable_reason <- if (refit$interval_status=="completed") "" else
          if (!out$summary$all_refits_succeeded) "Missing/failed prescribed draw or original point" else
            "Source count variance nonpositive"
        if(any(out$draws$status=="terminal_unproduced")) {
          refit$error <- state$error
          refit$interval_unavailable_reason <- paste("Terminal dataset error before prescribed count record:",state$error)
        }
        rows[[length(rows)+1L]] <- refit
      }
    }
  }
  list(mode=mode,records=.wm_ca_bind(rows),count_summaries=.wm_ca_bind(counts),
       draw_status=.wm_ca_bind(draw_status),requested_grid=grid,records_roots=roots,
       dataset_execution=execution,shard_completion_valid=shard_complete,
       execution_complete=all(execution$status %in% c("completed","failed")) && all(shard_complete),
       assessment_ready=mode=="calibration" && all(execution$status %in% c("completed","failed")) && all(shard_complete),
       note="Original-refit variance is qualified by all-prescribed-draw success; no failed-draw completion.")
}

.wm_ca_binomial <- function(success,n) {
  if (!n) return(c(rate=NA_real_,mcse=NA_real_,lower=NA_real_,upper=NA_real_))
  p <- success/n; ci <- stats::binom.test(success,n)$conf.int
  c(rate=p,mcse=sqrt(p*(1-p)/n),lower=ci[1L],upper=ci[2L])
}

# Same-set root-variance comparison. v's influence includes divisor R-1.
wm_cal_variance_comparison <- function(x,V,n,inner_mcse=rep(NA_real_,length(x))) {
  R <- length(x)
  out <- list(common_n=R,sampling_root_variance=NA_real_,sampling_root_variance_mcse=NA_real_,
    mean_reported_root_variance=.wm_ca_mean(V),mean_reported_variance_outer_mcse=.wm_ca_mcse(V),
    mean_reported_variance_inner_mcse=if(R && all(is.finite(inner_mcse))) sqrt(sum(inner_mcse^2))/R else NA_real_,
    variance_ratio=NA_real_,variance_ratio_outer_mcse=NA_real_,variance_ratio_inner_mcse=NA_real_,
    variance_difference=NA_real_,variance_difference_outer_mcse=NA_real_)
  if (R<2L) return(out)
  q <- (x-mean(x))^2; v <- n*stats::var(x); iv <- n*R/(R-1)*(q-mean(q)); iV <- V-mean(V)
  out$sampling_root_variance <- v; out$sampling_root_variance_mcse <- .wm_ca_mcse(iv)
  out$variance_difference <- mean(V)-v
  out$variance_difference_outer_mcse <- .wm_ca_mcse(iV-iv)
  if (v>0) {
    out$variance_ratio <- mean(V)/v
    out$variance_ratio_outer_mcse <- .wm_ca_mcse(iV/v-mean(V)*iv/v^2)
    out$variance_ratio_inner_mcse <- out$mean_reported_variance_inner_mcse/v
  }
  out
}

wm_cal_summarize_cell <- function(d) {
  .wm_ca_check(length(unique(d$n))==1L && length(unique(d$truth))==1L,"Mixed n or truth in cell.")
  p <- d$point_status=="completed" & is.finite(d$estimate)
  valid <- p & d$interval_status=="completed" & is.finite(d$root_n_variance) & d$root_n_variance>0
  x <- d$estimate[p]; e <- x-d$truth[1L]; R <- length(x); requested <- nrow(d)
  bias <- .wm_ca_mean(e); rmse <- if(R) sqrt(mean(e^2)) else NA_real_
  point_var_mcse <- if(R>1L)d$n[1L]*R/(R-1)*.wm_ca_mcse((x-mean(x))^2) else NA_real_
  coverage <- .wm_ca_binomial(sum(d$covered[valid]),sum(valid))
  reliable <- .wm_ca_binomial(sum(d$covered[valid]),requested)
  out <- list(requested=requested,point_successes=R,interval_successes=sum(valid),
    mean_estimate=.wm_ca_mean(x),mean_estimate_mcse=.wm_ca_mcse(x),signed_bias=bias,bias_mcse=.wm_ca_mcse(x),
    signed_relative_bias_percent=if(d$truth[1L]!=0)100*bias/d$truth[1L] else NA_real_,
    relative_bias_mcse_percent=if(d$truth[1L]!=0)100*.wm_ca_mcse(x)/abs(d$truth[1L]) else NA_real_,
    rmse=rmse,rmse_mcse=if(is.finite(rmse)&&rmse>0).wm_ca_mcse(e^2)/(2*rmse) else NA_real_,
    empirical_sd=sqrt(.wm_ca_var(x)),all_point_sampling_root_variance=d$n[1L]*.wm_ca_var(x),
    all_point_sampling_root_variance_mcse=point_var_mcse,
    empirical_sd_mcse=if(R>1L && .wm_ca_var(x)>0) point_var_mcse/(2*d$n[1L]*sqrt(.wm_ca_var(x))) else NA_real_,
    mean_reported_se=.wm_ca_mean(d$se[valid]),mean_reported_se_mcse=.wm_ca_mcse(d$se[valid]),
    mean_interval_length=.wm_ca_mean(d$upper[valid]-d$lower[valid]),
    mean_interval_length_mcse=.wm_ca_mcse(d$upper[valid]-d$lower[valid]),
    coverage=coverage[["rate"]],coverage_mcse=coverage[["mcse"]],coverage_lower=coverage[["lower"]],coverage_upper=coverage[["upper"]],
    success_and_cover=reliable[["rate"]],success_and_cover_mcse=reliable[["mcse"]],
    success_and_cover_lower=reliable[["lower"]],success_and_cover_upper=reliable[["upper"]],
    point_summary_scope=if(all(p))"all_requested" else "conditional_on_point_success",
    interval_summary_scope=if(all(valid))"all_requested" else "conditional_on_interval_success")
  c(out,wm_cal_variance_comparison(d$estimate[valid],d$root_n_variance[valid],d$n[1L],d$inner_variance_mcse[valid]))
}

wm_cal_summarize <- function(collected) {
  .wm_ca_check(identical(collected$mode,"calibration"),
               "Performance aggregation is prohibited for runtime-pilot records.")
  d <- collected$records
  fields <- c("branch","model","estimand","M","n")
  cells <- split(seq_len(nrow(d)),.wm_ca_key(d,fields))
  summary <- .wm_ca_bind(lapply(cells,function(j) cbind(d[j[1L],fields,drop=FALSE],
    as.data.frame(wm_cal_summarize_cell(d[j,,drop=FALSE]),stringsAsFactors=FALSE))))
  membership <- d[c("dataset_id",fields)]
  membership$point_success <- d$point_status=="completed" & is.finite(d$estimate)
  membership$interval_and_variance_common <- membership$point_success & d$interval_status=="completed" &
    is.finite(d$root_n_variance) & d$root_n_variance>0
  failures <- list()
  for (j in seq_len(nrow(d))) {
    if (!membership$point_success[j]) failures[[length(failures)+1L]] <- cbind(d[j,c("dataset_id",fields)],stage="point",reason=d$error[j])
    if (!membership$interval_and_variance_common[j]) failures[[length(failures)+1L]] <- cbind(d[j,c("dataset_id",fields)],stage="interval",reason=d$interval_unavailable_reason[j])
  }
  pairs <- list(); pair_members <- list()
  base <- d[d$branch=="actual_fitted",,drop=FALSE]
  for (other in c("fitted_naive","original_refit")) {
    alt <- d[d$branch==other,,drop=FALSE]
    jf <- c("dataset_id","model","estimand","M","n")
    joined <- match(.wm_ca_key(alt,jf),.wm_ca_key(base,jf)); .wm_ca_check(!anyNA(joined),"Unpaired comparison record.")
    a <- base[joined,,drop=FALSE]
    valid <- a$interval_status=="completed" & alt$interval_status=="completed" &
      a$point_status=="completed" & alt$point_status=="completed"
    group <- split(seq_len(nrow(alt)),.wm_ca_key(alt,c("model","estimand","M","n")))
    for (ids in group) {
      j <- ids[valid[ids]]; R <- length(j)
      row <- alt[ids[1L],c("model","estimand","M","n"),drop=FALSE]
      left <- wm_cal_variance_comparison(a$estimate[j],a$root_n_variance[j],row$n)
      right <- wm_cal_variance_comparison(alt$estimate[j],alt$root_n_variance[j],row$n,alt$inner_variance_mcse[j])
      names(left) <- paste0("analytic_",names(left)); names(right) <- paste0("comparison_",names(right))
      vdiff <- alt$root_n_variance[j]-a$root_n_variance[j]
      cdiff <- as.numeric(alt$covered[j])-as.numeric(a$covered[j])
      pairs[[length(pairs)+1L]] <- cbind(row,comparison=other,requested_pairs=length(ids),joint_valid=R,
        coverage_difference=.wm_ca_mean(cdiff),coverage_difference_mcse=.wm_ca_mcse(cdiff),
        mean_variance_difference=.wm_ca_mean(vdiff),mean_variance_difference_outer_mcse=.wm_ca_mcse(vdiff),
        as.data.frame(c(left,right),stringsAsFactors=FALSE))
      if (R) pair_members[[length(pair_members)+1L]] <- cbind(alt[j,c("dataset_id","model","estimand","M","n")],comparison=other)
    }
  }
  inner <- collected$count_summaries; inner_summary <- list()
  if (nrow(inner)) for (j in split(seq_len(nrow(inner)),.wm_ca_key(inner,c("model","estimand","M","n")))) {
    x <- inner[j,,drop=FALSE]; good <- x$all_refits_succeeded & x$all_one_step_available
    v <- x$paired_root_variance_difference[good]; s <- x$paired_variance_difference_mcse[good]; R <- length(v)
    inner_summary[[length(inner_summary)+1L]] <- cbind(x[1L,c("model","estimand","M","n")],
      requested_datasets=nrow(x),joint_complete_count_sets=R,
      mean_refit_minus_onestep_variance=.wm_ca_mean(v),outer_mcse_including_inner=.wm_ca_mcse(v),
      inner_mcse_component=if(R)sqrt(sum(s^2))/R else NA_real_,
      mean_root_difference=.wm_ca_mean(x$root_difference_mean[good]),
      mean_dataset_root_difference_rms=.wm_ca_mean(x$root_difference_rms[good]))
  }
  list(summary=summary,membership=membership,failures=.wm_ca_bind(failures),
    paired_comparisons=.wm_ca_bind(pairs),paired_membership=.wm_ca_bind(pair_members),
    count_pair_summary=.wm_ca_bind(inner_summary),count_summaries=inner,
    draw_status=collected$draw_status,
    dataset_execution=collected$dataset_execution,
    execution_status=if(isTRUE(collected$execution_complete))"all_requested_datasets_terminal" else
      "provisional_missing_or_active_datasets",
    assessment_ready=isTRUE(collected$assessment_ready),
    notes=c("No PSM reference was allocated; relative efficiency is not defined.",
      "All successful-set memberships are explicit. Missing requested records remain in requested denominators.",
      "Variance ratios compare each point estimator with its own sampling variance on the identical dataset set.",
      "Outer MCSE includes inner-count noise; the separately reported inner component is not added to it.",
      "One-step exact variance is an algebraic identity, not an independently validated sampling method.",
      "assessment_ready denotes complete internally consistent record accounting only, not scientific acceptance.",
      "No unconditional variance convergence or failure completion is inferred from these summaries."))
}
