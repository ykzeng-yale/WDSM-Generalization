script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)),mustWork=TRUE)
root <- normalizePath(file.path(dirname(script),".."))
for (name in c("first_stage.R","supplement.R","extension_production.R")) source(file.path(root,"simulations",name))
config <- utils::read.csv(file.path(root,"simulations","first_stage_production.csv"),stringsAsFactors=FALSE)
first <- wm_extension_plan("first_stage",config)
supp <- wm_extension_plan("supplement",wm_supplement_config())
for (plan in list(first,supp)) {
  stopifnot(!anyDuplicated(plan$result_file),identical(plan$task_id,seq_len(nrow(plan))))
  for (scenario in unique(plan$scenario)) {
    x <- plan[plan$scenario == scenario,]
    ids <- unlist(Map(seq.int,x$first,x$last),use.names=FALSE)
    stopifnot(!anyDuplicated(ids),identical(sort(ids),seq_len(x$replications[1])))
  }
}
stopifnot(sum(first$requested_records)==234000L,sum(supp$requested_records)==10800L,
  sum(first$stage=="timing")==4L,sum(supp$stage=="imported_pilot")==8L,
  sum(supp$requested_records[supp$stage=="imported_pilot"])==270L)
bad <- config;bad$seed[1] <- bad$seed[1]+1000L
stopifnot(inherits(try(wm_extension_plan("first_stage",bad),silent=TRUE),"try-error"))
bad <- config;bad$d[1] <- 1L
stopifnot(inherits(try(wm_extension_plan("first_stage",bad),silent=TRUE),"try-error"))
bad <- wm_supplement_config();bad$replications[1] <- 201L
stopifnot(inherits(try(wm_extension_plan("supplement",bad),silent=TRUE),"try-error"))
cat("PASS: first-stage",nrow(first),"tasks and supplemental",nrow(supp),"tasks cover every planned replication exactly once, preserve pilot IDs, and reject altered seeds/scope.\n")
