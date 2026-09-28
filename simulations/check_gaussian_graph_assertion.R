#!/usr/bin/env Rscript
# Focused deterministic harness check; no production dataset or random draws.
check_gaussian_graph_assertion <- function(helper_dir) {
source(file.path(helper_dir, "first_stage.R"), local = TRUE)
checks <- 0L
check <- function(value) { stopifnot(isTRUE(value)); checks <<- checks + 1L }
graph <- function(scores, arm = 0L) {
  delta <- scores[2L,,drop=FALSE] - scores[1L,,drop=FALSE]
  largest <- max(abs(delta)); distance <- largest * sqrt(sum((delta/if(largest>0) largest else 1)^2))
  list(edges=data.frame(query=1L,donor=2L,arm=arm,distance=distance,share=1,outcome_share=1),
       neighbors=matrix(2L,1L,1L),unordered=FALSE,cell=rep(1L,2),fold_id=rep(1L,2),strata=rep(1L,2),
       scores0=scores,scores1=if(arm==1L) scores else NULL,
       transform="Supplied coordinates; no automatic scaling",tie_rule="Exact distance, then original row index")
}
a <- graph(matrix(c(.7768372525744135,.777271207325254),ncol=1L))
s <- a$scores0; s[1L,1L] <- s[1L,1L] + .Machine$double.eps / 2
b <- graph(s)
check(!isTRUE(all.equal(a,b,tolerance=1e-13)))
check(.wm_fs_gaussian_same_graph(a,b))
check(.wm_fs_gaussian_same_graph(b,a))
check(.wm_fs_gaussian_same_graph(a,a))
for(name in c("query","donor","arm","share","outcome_share")) {
  bad <- b; bad$edges[[name]][1L] <- bad$edges[[name]][1L]+1L
  check(!.wm_fs_gaussian_same_graph(a,bad))
}
for(name in c("neighbors","unordered","cell","fold_id","strata","transform","tie_rule")) {
  bad <- b; bad[[name]] <- NULL
  check(!.wm_fs_gaussian_same_graph(a,bad))
}
bad <- b; bad$edges$distance <- bad$edges$distance + .Machine$double.eps
check(!.wm_fs_gaussian_same_graph(a,bad))
bad <- b; bad$edges$distance <- -bad$edges$distance
check(!.wm_fs_gaussian_same_graph(a,bad))
bad <- b; bad$edges$share <- as.integer(bad$edges$share)
check(!.wm_fs_gaussian_same_graph(a,bad))
s <- a$scores0 + .001; bad <- graph(s)
check(!.wm_fs_gaussian_same_graph(a,bad))
# Equal scores/ties and arm1/full-dimensional paths remain supported.
tie <- graph(matrix(c(0,0),ncol=1L));check(.wm_fs_gaussian_same_graph(tie,tie))
a2 <- graph(cbind(c(.7768372525744135,.777271207325254),c(.2,.20003)),arm=1L)
s <- a2$scores0;s[1L,1L] <- s[1L,1L] + .Machine$double.eps/2
b2 <- graph(s,arm=1L)
check(.wm_fs_gaussian_same_graph(a2,b2))
bad <- b2;bad$scores1 <- bad$scores1 + .01
check(!.wm_fs_gaussian_same_graph(a2,bad))
bad <- b; colnames(bad$scores0) <- "renamed"
check(!.wm_fs_gaussian_same_graph(a,bad))
bad <- b; rownames(bad$edges) <- "renamed"
check(!.wm_fs_gaussian_same_graph(a,bad))
bad <- b; names(bad$edges$distance) <- "renamed"
check(!.wm_fs_gaussian_same_graph(a,bad))
for(nonfinite in c(NA_real_, Inf, -Inf, NaN)) {
  both <- a;both$scores0 <- rbind(both$scores0,nonfinite)
  check(!.wm_fs_gaussian_same_graph(both,both))
}
cat("PASS", checks, "Gaussian graph assertion checks; no random draws or production rerun.\n")

invisible(checks)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  stopifnot(length(args) == 1L)
  check_gaussian_graph_assertion(args[1L])
}
