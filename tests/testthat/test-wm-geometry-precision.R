test_that("angular geometry marks zero empirical variance unresolved without changing estimates", {
  fixture <- NULL
  local_mocked_bindings(.wm_geometry_mc = function(...) fixture,
                        .package = "WeightedMatching")

  cases <- list(
    list(beta = 2, mcse = 0, hits = 2,
         status = "unresolved_zero_empirical_variance"),
    list(beta = 0, mcse = 0, hits = 0,
         status = "unresolved_zero_hits"),
    list(beta = 2, mcse = 0.5, hits = 2,
         status = "monte_carlo_estimate")
  )
  for (case in cases) {
    fixture <- list(beta = case$beta, beta_mcse = case$mcse,
                    covariance = matrix(case$mcse^2, 1, 1),
                    alpha_mcse = case$mcse, nonzero = case$hits,
                    max_value = case$beta)
    result <- wm_geometry(1, 2, draws = 2)
    expect_identical(unname(result$component_precision), case$status)
    expected_status <- if (grepl("^unresolved", case$status))
      "unresolved_components" else "monte_carlo_estimate"
    expect_identical(result$precision_status, expected_status)
    expect_identical(unname(result$beta), case$beta)
    expect_identical(result$alpha, case$beta)
    expect_identical(unname(result$beta_mcse), case$mcse)
    expect_identical(result$alpha_mcse, case$mcse)
    expect_identical(unname(result$covariance), fixture$covariance)
    expect_identical(unname(result$diagnostics$nonzero_draws), case$hits)
  }
})
