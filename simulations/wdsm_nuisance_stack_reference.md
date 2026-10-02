# Original WDSM nuisance preparation and its joint influence

This standalone reference independently implements the original retrospective weighted logistic propensity fit, ordinary arm-specific prognostic regressions, pooled unweighted standardization, and weighted full quadratic bias correction. It does not import predecessor predictions or export a package API. Its full estimating-equation Jacobian includes generated predictors and every nuisance cross block.

The finite-dimensional system and parameter order are implemented in [wdsm_nuisance_stack_reference.R](wdsm_nuisance_stack_reference.R) and described below. The public [fitted-pipeline contract](wdsm_fitted_pipeline_reference.md) specifies the separately declared conditions required for inference; this nuisance preparation returns `assumptions_verified = FALSE`. In particular, correct PG gives the full-X centering route in the bounded iid source-form CorCor/MisCor analogue; accepting a design matrix is not verification of that route in other data.

## Inputs and outputs

Source `simulations/wdsm_nuisance_stack_reference.R` and call

```r
stack <- wm_wdsm_nuisance_stack_reference(
  Y, Z, weights, ps_design = X_ps,
  pg0_design = X_pg, pg1_design = X_pg, estimand = "PATE",
  ps_offset = offset, multiplicity = rep(1, length(Y)))
```

Supply raw model matrices with the intended intercept columns already present. The propensity model uses all rows; each prognostic model uses only its corresponding arm to estimate coefficients and predicts on every row. The two PG matrices may differ. For PATT omit `pg1_design` and `pg1_offset`; only the control PG and correction are fitted. Offsets are fixed known numeric vectors; `ps_offset`, `pg0_offset`, and `pg1_offset` default to zero. PG offsets are optional fixed regression offsets, not estimated score coordinates.

`multiplicity` accepts finite nonnegative real masses, including zeros, with positive mass in each arm. This permits supplied multinomial counts or deterministic empirical perturbations. The same mass multiplies every estimating equation. All pooled moments divide by total mass; with unit masses this is the original divisor n. The code neither generates counts nor certifies their probability law.

The result includes:

- `parameter`, `parameter_names`, and `blocks`: one fixed column order for all derivatives and influences;
- `raw_scores`, `scores0`, `scores1`, `mean0`, and `mean1`: raw PS/PG predictions, the two standardized arm maps, and quadratic correction predictions;
- `estimating_equations` (n by p), `mean_equation`, `jacobian` (p by p), `nuisance_influence` (n by p), and `nuisance_covariance` (p by p);
- `mean0_derivative`, `mean1_derivative`, and `score_derivatives`: derivatives with respect to the full parameter vector, evaluated with observed rows fixed;
- `weight_derivative`: an n-by-p zero matrix because individual weights are known here;
- fitting/root diagnostics and `assumptions_verified = FALSE`.

PATT's unused arm-1 outputs are NULL. The PS center/variance is shared across arms, so its moment equations appear once. Each PG has its own pooled center/variance. This removes duplicate, perfectly correlated moment coordinates without changing either standardized map or its total influence. Parameter dimension p counts fitted coefficients and moments; each matching map remains two-dimensional.

The arbitrary-parameter evaluator

```r
at_theta <- wm_wdsm_nuisance_evaluate_reference(stack, parameter)
```

returns all equations, their analytic averaged Jacobian, scores, predictions, and derivatives without fitting. Named parameter vectors must use exactly the stored order. This interface supports independent central differences of the *entire* stack, including generated-score terms.

## Equations and full chain rule

Let a_i be the supplied masses and `pi_i=a_i/sum(a)`. Put `P_a g=sum_i pi_i g_i`. A fixed common scale `c_W=max_i W_i` is used for arithmetic and `w_i=W_i/c_W`. Multiplying every propensity/correction equation by this common nonzero scalar does not change their root or `-J^{-1}psi`; it is held fixed when differentiating. It is not a new estimated-weight coordinate.

Write x_i for the PS design, g_zi for a PG design, and the fixed offsets as o_P,i and o_z,i. The raw coordinates are

\[
 e_i=\operatorname{expit}(o_{P,i}+x_i^\top\alpha),\qquad
 t_{zi}=o_{z,i}+g_{zi}^\top\beta_z.
\]

Let `T_i=(e_i,t_0i,t_1i)` for PATE, omitting t_1 for PATT. For each raw coordinate j, parameters c_j and v_j represent its pooled mean and variance. The own arm's standardized map is

\[
 D_{zi}=((e_i-c_P)/\sqrt{v_P},(t_{zi}-c_z)/\sqrt{v_z})^\top.
\]

With `B(D)=(1,D1,D2,D1²,D1 D2,D2²)'` and `f_zi=B(D_zi)'delta_z`, the unmultiplied row equations are

\[
\begin{aligned}
 \psi_{\alpha,i}&=w_i x_i(Z_i-e_i),\\
 \psi_{\beta_z,i}&=1(Z_i=z)g_{zi}(Y_i-t_{zi}),\\
 \psi_{c_j,i}&=T_{ji}-c_j,\\
 \psi_{v_j,i}&=(T_{ji}-c_j)^2-v_j,\\
 \psi_{\delta_z,i}&=1(Z_i=z)w_i B(D_{zi})(Y_i-f_{zi}).
\end{aligned}                                                   \tag{1}
\]

The fitted stack solves `P_a psi=0`. PG equations have no individual weight W. The propensity and quadratic correction do. The quadratic correction includes an intercept and all five nonconstant degree-at-most-two terms; it is a fixed six-term regression, not a growing sieve.

For any full parameter coordinate k,

\[
\begin{aligned}
 \partial_k e_i&=e_i(1-e_i)x_i^\top\partial_k\alpha,\\
 \partial_k D_{zi,r}
 &=\{\partial_k T_{zi,r}-\partial_k c_r\}/\sqrt{v_r}
          -D_{zi,r}\partial_k v_r/(2v_r),\\
 \partial_k B(D)&=(0,\partial_k D_1,\partial_k D_2,
 2D_1\partial_k D_1,D_2\partial_k D_1+D_1\partial_k D_2,
 2D_2\partial_k D_2)^\top,\\
 \partial_k f_{zi}&=\delta_z^\top\partial_k B(D_{zi})
                      +B(D_{zi})^\top\partial_k\delta_z.
\end{aligned}                                                   \tag{2}
\]

The variance parameter is v, not standard deviation; this explains the `-D/(2v)` derivative. The off-diagonal center/variance rows retain derivatives of the raw fitted scores. At arbitrary parameters, the variance-to-center derivative is `-2 P_a(T_j-c_j)`; it is evaluated rather than dropped merely because it vanishes at an exact moment root.

The correction block includes **both** generated-basis terms:

\[
 \partial_k\psi_{\delta_z,i}
 =1(Z_i=z)w_i\{(Y_i-f_{zi})\partial_k B(D_{zi})
                         -B(D_{zi})\partial_k f_{zi}\}.        \tag{3}
\]

The raw model, pooled-moment and correction order makes the Jacobian block triangular. Its diagonal regression blocks are negative weighted Gram/information matrices; moment blocks have negative identities. When these are nonsingular, the full Jacobian is invertible even if the resulting influence covariance is singular.

The implementation returns

\[
 \widehat J=P_a\partial_\theta\psi_i,\qquad
 \widehat\ell_i=-\widehat J^{-1}\psi_i,\qquad
 \widehat\Sigma=P_a\{(\widehat\ell_i-P_a\widehat\ell)
                 (\widehat\ell_i-P_a\widehat\ell)^\top\}.       \tag{4}
\]

This is a complete sandwich influence. Treatment-dependent individual weights do not justify substituting an ordinary logistic Fisher-information inverse for the covariance. The off-diagonal entries in (4) retain PS/PG/moment/correction covariance and the generated-regressor dependence. For unit masses, this is the usual empirical iid influence covariance. For nonunit masses, it is the functional influence and covariance under that weighted empirical measure, not a claim that the original n rows are iid from a new law.

For a perturbation `a_i(t)=a_i(1+t r_i)`, with positive nearby masses and `P_a r=0`, implicit differentiation gives

\[
 \left.\frac{d\widehat\theta(t)}{dt}\right|_{t=0}
    =P_a(\widehat\ell_i r_i).                                \tag{5}
\]

At unit masses this is `mean(ell_i*r_i)`. Equation (5) differentiates nuisance fitting, not nearest-neighbor selection. Central differences of refits under two distinct signed zero-mean directions check the complete influence chain; derivatives of predictions then multiply (5) by the full row prediction Jacobian in (2).

## Connection to fitted matching components

For unit-mass preparation, construct the original corrected point estimate with the returned maps and means:

```r
fit <- wdsmatch::wm_match(
  Y, Z, weights, scores0 = stack$scores0, scores1 = stack$scores1,
  mean0 = stack$mean0, mean1 = stack$mean1,
  estimand = stack$estimand, variance = FALSE)
smooth <- wm_graph_gradient_reference(
  fit, weight_derivative = stack$weight_derivative,
  mean0_derivative = stack$mean0_derivative,
  mean1_derivative = stack$mean1_derivative)
```

The matching gradient holds the full parameter coordinates independent while differentiating predictions. The fitted nuisance dependence then enters through `stack$nuisance_influence`. Do not substitute a derivative that has already multiplied by this influence, and do not add center/variance effects twice.

Supply `stack$nuisance_influence` and `smooth$sensitivity` to `wm_fitted_variance_reference`, with a separately justified graph sensitivity and covariance scope. In the proved source-form bounded iid CorCor/MisCor full-X branch, the graph sensitivity is zero. A correct target PS alone does not establish that branch when PG is misspecified. No code here estimates a nonzero graph-transport derivative, certifies arbitrary matching dimensions, or validates the original refit bootstrap outside its variance-agreement theorem.

## Numerical and validation boundaries

The fitting equations and model bases match the pinned original source. This independent reference uses a common weight scale and a default GLM convergence tolerance of `1e-12` rather than the predecessor's `1e-10`; these select the same statistical roots, not necessarily bitwise identical solver paths. Logistic fits begin at zero coefficients and use the same quasibinomial logit family. Rank failures, saturated scores, degenerate standardization, and failed roots are explicit errors. A convergence diagnostic does not prove the asymptotic requirement that numerical parameter error is negligible compared with n^(-1/2).

Required focused validation consists of central differences of `mean_equation` and both prediction derivatives, two deterministic mass perturbation/refit comparisons using (5), and independent preparation against pinned predecessor functions under aligned model/offset/mass conventions. These are algebra and source-fidelity checks. They do not estimate bias, coverage, sampling variance, or Monte Carlo performance. The historical finite-population Gaussian simulation still needs its own applicability assessment.
