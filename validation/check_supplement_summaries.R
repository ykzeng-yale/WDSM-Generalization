script<-normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)))
root<-normalizePath(file.path(dirname(script),".."))
for (name in c("first_stage.R","first_stage_summarize.R","supplement.R","supplement_summarize.R")) source(file.path(root,"simulations",name))
equal<-function(x,y) {
  result<-all.equal(x,y,tolerance=1e-10,check.attributes=FALSE)
  if (!isTRUE(result)) stop(paste(result,collapse="; "))
}
fails<-function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
cfg<-wm_supplement_config()[1,,drop=FALSE]
x<-expand.grid(replication=1:4,method=c("self_normalized","stabilized"),estimand=c("PATE","PATT"),correction=c("oracle","feasible"),stringsAsFactors=FALSE)
x$scenario<-cfg$scenario;x$branch<-cfg$branch;x$n<-x$analysis_n<-cfg$n;x$d<-cfg$d;x$M<-cfg$M
x$fraction1<-NA_real_;x$base_seed<-cfg$seed;x$target<-1.5;x$status<-"ok"
x$estimate<-1.5+c(-.3,-.1,.1,.3)[x$replication]+.1*(x$correction=="feasible")
x$sampling_variance<-.08+.01*x$replication;x$root_n_variance<-x$analysis_n*x$sampling_variance
x$lower<-x$estimate-qnorm(.975)*sqrt(x$sampling_variance);x$upper<-x$estimate+qnorm(.975)*sqrt(x$sampling_variance)
x$clipped<-0L;x$rng_state<-paste0("manual/",x$replication);x$warning<-x$error<-""
s<-wm_supplement_summarize(x)
stopifnot(nrow(s)==8L,all(s$requested==4L),all(s$paired_completed==4L))
equal(s$bias,ifelse(s$correction=="feasible",.1,0));equal(s$empirical_root_n_variance,rep(800/15,8))
equal(s$mean_root_n_variance,rep(84,8));equal(s$feasible_minus_oracle_point,rep(.1,8))
equal(s$feasible_minus_oracle_point_mcse,rep(0,8));equal(s$coverage_exact95_lower,rep(.025^.25,8))
stopifnot(all(s$coverage_exact95_upper==1))
fails(wm_supplement_summarize(x[-1,]));fails(wm_supplement_summarize(rbind(x,x[1,])))
y<-x;y$analysis_n[1]<-799;fails(wm_supplement_summarize(y))
y<-x;y$rng_state[1]<-"different";fails(wm_supplement_summarize(y))
y<-x;bad<-y$replication==1 & y$method=="stabilized" & y$estimand=="PATE" & y$correction=="feasible"
y$status[bad]<-"error";y$error[bad]<-"intentional fixture failure"
for (f in c("estimate","root_n_variance","sampling_variance","lower","upper")) y[[f]][bad]<-NA_real_
z<-wm_supplement_summarize(y)
stopifnot(z$completed[z$method=="stabilized" & z$estimand=="PATE" & z$correction=="feasible"]==3,
  all(z$paired_completed[z$method=="stabilized" & z$estimand=="PATE"]==3))
directory<-tempfile("supplement-fixture-");dir.create(directory);path<-file.path(directory,"records.csv")
write.csv(x,path,row.names=FALSE)
hash<-setNames(paste(rep("a",32),collapse=""),"/test/a.R")
meta<-list(schema_version="supplement_v1",config=cfg,first=1L,last=4L,completed=Sys.time(),code_unchanged=TRUE,
  requested_records=32L,records_written=32L,columns=names(x),last_completed_replication=4L,code_md5=hash,design_md5=hash)
saveRDS(meta,paste0(path,".metadata.rds"))
read<-wm_supplement_read_batches(path,FALSE);equal(wm_supplement_summarize(read),s)
fails(wm_supplement_read_batches(path,TRUE));fails(wm_supplement_read_batches(c(path,path),FALSE))
meta$code_unchanged<-FALSE;saveRDS(meta,paste0(path,".metadata.rds"));fails(wm_supplement_read_batches(path,FALSE))
unlink(directory,recursive=TRUE)
cat("PASS: supplemental summary manual moments, exact coverage interval, paired failures, complete record grids, sample-size/RNG guards and CSV/provenance round trip.\n")
