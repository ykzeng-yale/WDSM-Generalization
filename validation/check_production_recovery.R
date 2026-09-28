# Only deterministic complete-prefix safety fixtures; no simulation/job execution.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(script), "..", "simulations", "production_recover.R"))
task <- data.frame(scenario = "fixture", n = 20, d = 1, M = 1, design = "strong", seed = 123, first = 11L, last = 15L)
x <- expand.grid(replication = 11:15, method = c("self_normalized", "stabilized"),
  estimand = c("PATE", "PATT"), correction = c("oracle", "polynomial", "raw"), stringsAsFactors = FALSE)
for (f in c("scenario", "n", "d", "M", "design", "seed")) x[[f]] <- task[[f]]
x$status <- "ok"; x$status[3] <- "error"
stopifnot(.wm_recovery_grid(x, task) == 15L,
          .wm_recovery_grid(x[x$replication < 14, ], task, TRUE) == 13L)
fails <- function(z, pattern, allow = TRUE) {
  message <- tryCatch({.wm_recovery_grid(z, task, allow); NULL}, error = conditionMessage)
  stopifnot(is.character(message), grepl(pattern, message))
}
fails(x[-1, ], "incomplete replication")
fails(x[x$replication != 12, ], "contiguous prefix")
fails(x[x$replication < 14, ], "contiguous prefix", FALSE)
fails(rbind(x, x[1, ]), "duplicated")
y <- x; y$seed[1] <- 456; fails(y, "configuration mismatch")
y <- x; y$replication[y$replication == 15] <- 16; fails(y, "contiguous prefix")
cat("PASS: exact contiguous prefixes, retained estimator failures, missing interior/tail records, duplicate IDs and configuration guards.\n")
