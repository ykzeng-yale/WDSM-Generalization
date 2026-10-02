# Validate concrete historical count inputs for a compatibility replay only.
# This function does not generate counts or change the RNG state.
lenis_canonical_counts <- function(bundle, key, expected, n, B, seed, rng_kind, rng_before) {
  stopifnot(identical(bundle$purpose, "historical_count_input_replay"),
    is.character(key), length(key) == 1L, key %in% names(bundle$cases),
    !anyDuplicated(names(bundle$cases)))
  item <- bundle$cases[[key]]
  x <- item$matrix
  stopifnot(identical(item$metadata, expected),
    identical(as.integer(seed), as.integer(expected$seed)),
    identical(rng_kind, expected$rng_kind), identical(rng_before, expected$rng_before),
    identical(as.integer(n), as.integer(expected$n)),
    identical(as.integer(B), as.integer(expected$B)),
    is.matrix(x), is.integer(x), identical(dim(x), as.integer(c(n, B))),
    is.null(dimnames(x)), !anyNA(x), all(x >= 0L), all(colSums(x) == n),
    identical(digest::digest(x, algo = "sha256"), expected$sha256),
    identical(vapply(seq_len(B), function(b) digest::digest(x[, b], algo = "sha256"),
                     character(1)), expected$column_sha256))
  list(counts = x, origin_metadata = item$metadata,
    provenance = list(mode = "supplied_historical_count_matrix",
      case = key, origin = bundle$origin, origin_metadata_sha256 = item$metadata_file_sha256,
      matrix_sha256 = expected$sha256, generated_by_this_replay = FALSE,
      rng_after_action = "Restore the canonical historical post-count state for input coupling; not a claim of Linux count generation"))
}
