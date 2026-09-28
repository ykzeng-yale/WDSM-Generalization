# Deterministic integrity fixtures; no estimator, simulator, cluster or job runs.
script<-normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)),mustWork=TRUE)
root<-normalizePath(file.path(dirname(script),".."))
for (name in c("first_stage.R","supplement.R","extension_production.R","extension_collect.R"))
  source(file.path(root,"simulations",name))
fails<-function(expression,pattern) {
  error<-tryCatch({force(expression);NULL},error=conditionMessage)
  if (is.null(error) || !grepl(pattern,error)) stop("Expected error: ",pattern,"; observed: ",error)
}
check<-function() {
  work<-tempfile("extension-collector-");dir.create(work)
  on.exit(unlink(work,recursive=TRUE))
  # Same path-to-hash mapping under different host collation is identical;
  # changed/missing sources remain different and duplicate paths are invalid.
  hashes<-setNames(c(paste(rep("a",32),collapse=""),paste(rep("b",32),collapse="")),
    c("R/wm_geometry.R","R/wm_geometry_alpha.R"))
  stopifnot(identical(.wm_extension_collect_hashes(hashes),.wm_extension_collect_hashes(rev(hashes))))
  changed_hashes<-hashes;changed_hashes[1]<-paste(rep("c",32),collapse="")
  stopifnot(!identical(.wm_extension_collect_hashes(hashes),.wm_extension_collect_hashes(changed_hashes)),
    !identical(.wm_extension_collect_hashes(hashes),.wm_extension_collect_hashes(hashes[-1])))
  fails(.wm_extension_collect_hashes(c(hashes,hashes[1])),"Invalid frozen source hash")
  first<-utils::read.csv(file.path(root,"simulations","first_stage_production.csv"),stringsAsFactors=FALSE)
  for (family in c("first_stage","supplement")) {
    config<-if (family=="first_stage") first else wm_supplement_config()
    plan<-wm_extension_plan(family,config)
    csv<-file.path(work,paste0(family,".csv"))
    utils::write.csv(plan,csv,row.names=FALSE)
    frozen<-list(family=family,config=config,plan=plan,plan_md5=unname(tools::md5sum(csv)))
    stopifnot(identical(.wm_extension_collect_plan(family,frozen,csv),plan))
    # Regression: omitting all tasks of one first-stage cell previously passed
    # the CSV-hash check and the per-present-scenario completeness checks.
    omitted<-plan$scenario==tail(config$scenario,1L)
    shortened<-frozen;shortened$plan<-plan[!omitted,,drop=FALSE]
    fails(.wm_extension_collect_plan(family,shortened,csv),"RDS plan.*reconstructed")
    shortened<-frozen;shortened$plan<-plan[-1L,,drop=FALSE]
    fails(.wm_extension_collect_plan(family,shortened,csv),"RDS plan.*reconstructed")
    shuffled<-frozen;shuffled$plan<-plan[rev(seq_len(nrow(plan))),,drop=FALSE]
    fails(.wm_extension_collect_plan(family,shuffled,csv),"RDS plan.*reconstructed")
    # An updated checksum cannot make a shortened CSV a complete declared plan.
    utils::write.csv(plan[!omitted,,drop=FALSE],csv,row.names=FALSE)
    changed<-frozen;changed$plan_md5<-unname(tools::md5sum(csv))
    fails(.wm_extension_collect_plan(family,changed,csv),"CSV differs.*reconstructed")
    fails(.wm_extension_collect_plan(family,frozen,csv),"CSV hash differs")
    utils::write.csv(plan,csv,row.names=FALSE)
    altered<-frozen;altered$config$seed[1]<-altered$config$seed[1]+1000L
    fails(.wm_extension_collect_plan(family,altered,csv),if (family=="first_stage") "seeds changed" else "pre-specified")
  }
  # Imported pilot files must still be byte-identical after transfer/retrieval;
  # a valid completion status alone must not waive these recorded copy hashes.
  path<-file.path(work,"pilot.csv")
  writeLines("replication,estimate\n1,0.5",path)
  saveRDS(list(completed=TRUE),paste0(path,".metadata.rds"))
  saveRDS(list(retained=1L),paste0(path,".diagnostics.rds"))
  files<-paste0(path,c("",".metadata.rds",".diagnostics.rds"))
  metadata<-list(imported=TRUE,source_md5=tools::md5sum(files))
  task<-data.frame(stage="imported_pilot",stringsAsFactors=FALSE)
  stopifnot(identical(unname(.wm_extension_collect_import(path,task,metadata)),unname(metadata$source_md5)))
  writeLines("replication,estimate\n1,0.6",path)
  fails(.wm_extension_collect_import(path,task,metadata),"differs from its original hashes")
  writeLines("replication,estimate\n1,0.5",path)
  saveRDS(list(retained=2L),files[3])
  fails(.wm_extension_collect_import(path,task,metadata),"differs from its original hashes")
  unlink(files[3]);fails(.wm_extension_collect_import(path,task,metadata),"missing")
  missing<-metadata;missing$source_md5<-NULL
  fails(.wm_extension_collect_import(path,task,missing),"missing")
  missing<-metadata;missing$imported<-FALSE
  fails(.wm_extension_collect_import(path,task,missing),"missing")
  ordinary<-data.frame(stage="main",stringsAsFactors=FALSE)
  stopifnot(is.null(.wm_extension_collect_import(path,ordinary,list())))
  fails(.wm_extension_collect_import(path,ordinary,metadata),"Unexpected")
}
check()
cat("PASS: collector reconstructs complete grids, binds both plan representations, rejects omitted cells/changed seeds and verifies every retained pilot artifact hash.\n")
