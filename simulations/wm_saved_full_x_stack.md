# Complete influence inputs from an existing matching fit

`wm_saved_full_x_stack.R` constructs inputs for the common
`wm_fitted_inference()` interface from recorded fits. It does not fit a
model, build a graph, generate data or draw bootstrap counts. PS matching
and matching on any fixed number of supplied covariates use the same helper.
The original corrected point, weights, donor fractions and graph are retained.

Source the helper after installing `WeightedMatching` and pass its namespace
explicitly. This source-pinned helper retains its historical default namespace;
the explicit argument selects the separate current WM package:

```r
source("simulations/wm_saved_full_x_stack.R")

stack <- wm_saved_full_x_stack(
  data = list(Y = Y, Z = Z, weights = W),
  fit = recorded_ps_fit,
  mean_models = list(
    mean0 = list(design = D0, coefficients = beta0),
    mean1 = list(design = D1, coefficients = beta1)),
  family = "PS",
  wm = asNamespace("WeightedMatching"),
  ps = list(design = Dps, probability = recorded_probability,
            weighting = "probability"))

inference <- do.call(WeightedMatching::wm_fitted_inference,
                     stack$inference_arguments)
```

`recorded_ps_fit` is an existing corrected `wm_match()` object. Its scalar
map must equal the pooled standardized recorded probabilities. Declare the
original PS fitting weights: `"probability"` uses the supplied weights and
`"unit"` uses ones. This adapter covers the source zero-offset logistic model.
Optional recorded PS coefficients are checked against the probabilities;
missing coefficients remain unavailable. The complete influence still includes
every original logistic coefficient, pooled center and variance, and all
prediction coefficients. Each recorded mean design and coefficient vector
must reproduce the fit's predictions in the same row order.

For covariate matching, use `family = "full_X"`, omit `ps`, and supply
`raw_X`. Both recorded arm maps must equal its pooled standardized coordinates.
The mean models must use the complete polynomial basis of those coordinates.
The default basis order is intercept, raw-coordinate order, then products
`X_j * X_k` for `j <= k`. To use another fixed degree, supply every monomial
through that total degree as the ordered `polynomial_exponents` matrix, with
matching coordinate and basis names.

The helper converts standardized polynomial coefficients to the raw basis
algebraically. Complete polynomial spaces are invariant under invertible
affine changes of coordinates: fitted predictions agree for the original
sample and every full-rank count refit. Coordinate scaling can still change
multidimensional matching distances. The complete raw parameterization retains
all coordinate centers and variances; it does not discard their joint
influence merely because their prediction derivatives are zero.

For PATT, supply only `mean0`; the returned common-inference arguments use
`covariance_scope = "patt"` with explicit zero transport from current centering.
PATE supplies both means and uses `"full_x"`. The helper returns the complete
equations, Jacobian, influence matrix, covariance, score and mean derivatives,
and their common parameter order. It uses the existing inference implementation.

Full-design correction correctness, target identification, regular statistical
roots, numerical localization, current score geometry and the applicable
moment/actual-row conditions must be justified separately. If weights depend
on observed marks omitted from `raw_X`, they can still vary within matching
coordinates; retain their marked variance law. The helper does not apply an
unmarked full-covariate variance formula or validate a full/partial refit
callback. Original Gaussian and dependent survey designs remain empirical
comparison settings under the currently stated theorem scope.
