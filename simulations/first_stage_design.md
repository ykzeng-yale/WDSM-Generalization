# Pre-specified first-stage validation extension

**Historical design, outside-scope branch retained:** the current WM framework
uses supplied, known probability weights and studies estimated scores,
coordinate scales and outcome corrections. The estimated-weight section below
records an earlier optional development experiment; it is not part of current
WM theory or its required validation. Preserve its original configurations,
seeds, failures and results without starting new runs for this branch.

2026-09-28. This design is fixed before inspecting first-stage simulation
performance. It supplements the primary design without changing its
scripts or result schema. The experiments use synthetic iid data and
the qualified parametric nuisance theorems, not dependent survey designs.

## Score and prediction first stage

Generate \(S\sim {\rm Unif}([0,1]^d)\), \(U\in\{-1,1\}\) fair and
independent of \(S\), and
\[
 Q(Z=1\mid U)=(1+rU)/2,\qquad r=0.6.
\]
Weights are one. Let \(\epsilon=\pm0.25\) independently with equal
probabilities, \(\sigma^2=1/16\), and
\[
 Y=S_1+Z\tau+\epsilon,\quad \tau=1,\quad
 D=US_1(1-S_1).
\]
The original covariate vector includes \(U\). Define a scalar prediction
and matching-map parameter, true value zero:
\[
 S_\theta=(S_1+\theta D,S_2,\ldots,S_d),\qquad
 f_0=S_1+\theta D,\quad f_1=f_0+\tau.
\]
On the declared parameter interval \([-0.9,0.9]\), the first-coordinate
Jacobian is \(1+\theta U(1-2S_1)\ge0.1\). Each map is a smooth
support-preserving bijection of the cube, with bounded derivatives and
uniform current densities. The other coordinates are unchanged, and the
displacement vanishes on the first-coordinate boundary. The discrete
conditional \(U\) marks and transformed smooth densities give the
required uniform interior marked-law regularity. Squared distance
comparisons are quadratic polynomials in \(\theta\).

Fit the known-offset least-squares calibration
\[
 \widehat\theta=
 \frac{\sum D_i(Y_i-S_{i1}-Z_i\tau)}{\sum D_i^2},\qquad
 \widehat U_i^\theta=
 \frac{D_i\{\epsilon_i-\widehat\theta D_i\}}
      {N^{-1}\sum D_i^2}.
\]
Here \(N=n\) for same-sample fitting and \(N=m\) for independent training.
The known offset and known constant effect deliberately give an
analytically checkable validation fixture; this is not an applied
procedure for estimating an unknown treatment effect. A nonpositive
denominator, nonfinite arithmetic, or \(|\widehat\theta|\ge0.8\) is a
first-stage failure. It is recorded, never clipped, refitted or replaced.
The 0.8 rule leaves a declared margin inside the 0.9 parameter interval.

Both potential-mean derivatives are \(D\) on the evaluation rows.
Use same-sample dimensions \(d=2,3,5\); no same-sample scalar claim is
made. The planar branch uses known weights and its distributional
transfer theorem. The higher-dimensional branches satisfy polynomial
graph comparisons. Independent training allows \(d=1,2,3,5\), with
\(m=n\) initially and its own sample-size scaling and independent RNG
substream.

### Independent benchmark derivation

\[
 E(D^2)=\int_0^1s^2(1-s)^2\,ds=1/30,\qquad
 \Sigma_\theta=30\sigma^2.
\]
Conditional means of \(U\) are \(r\) in treated rows and \(-r\) in
control rows. Since \(\int_0^1s(1-s)\,ds=1/6\), the PATE direction
imbalances are \(b_1=-r/6\) and \(b_0=r/6\); PATT has
\(b_0=r/3\). Therefore
\[
 b_A=b_T=-r/3.
\]
The arm score densities and probabilities agree, so the unweighted
incoming load has \(EK=1\) and
\(EK^2=1/M+\alpha_d(M)/M^2\). For PATE the residual contribution is
\((2Z-1)(1+K)\epsilon\). For PATT it is
\(2\{Z-(1-Z)K\}\epsilon\). Multiplying each by
\(U^\theta=30D\epsilon\) gives
\[
 C_A=C_T=10r\sigma^2,\quad
 V_{0,A}=\sigma^2\{3+1/M+\alpha_d(M)/M^2\},\quad
 V_{0,T}=2\sigma^2\{1+1/M+\alpha_d(M)/M^2\}.
\]
The same-sample limit is
\[
 V_{\rm same}=V_0+2bC+b^2\Sigma_\theta
             =V_0-(10/3)r^2\sigma^2.
\]
For independent training it is
\[
 V_{\rm independent}=V_0+(n/m)b^2\Sigma_\theta
                   =V_0+(n/m)(10/3)r^2\sigma^2.
\]
These formulas concern sampling variance. The deliberately unadjusted
fitted contribution variance converges to \(V_0\), so the benchmark table
must distinguish its reported variance limit from its actual sampling
variance.

## Estimated individual weights on fixed maps

Generate \(X\sim{\rm Unif}([0,1]^d)\), set
\(t=X_1-1/2\), \(T=(1,t)^\top\), \(q=\operatorname{expit}(3t)\), and
draw \(Z\mid X\sim{\rm Bernoulli}(q)\). The known target arm proportion is
\(\pi=1/2\), with true weights
\[
 W_1=\pi/q,\qquad W_0=(1-\pi)/(1-q).
\]
The weighted target has uniform \(X\) in both arms. Let
\[
 \mu_0=1+X_1+\sum_jX_j^2/(2d),\quad
 \mu_1=\mu_0+1+X_1,\quad
 Y=\mu_Z+\epsilon,\quad \epsilon=\pm0.25.
\]
Both target effects equal 1.5. Use \(X\) as the fixed matching map and
the oracle mean predictions. Fit logistic weights with design \(T\)
and declared coordinatewise parameter bound 5; the true parameter is
\((0,3)\), strictly interior. Do not clip probabilities or use ridge,
Firth correction, term dropping, or retries. All fitting failures remain.

The analytic weight derivative is \(W(q-Z)T\), and the estimated
influence comes from the fitted empirical logistic information.
For the true information and sensitivity,
\[
 I=\int_{-1/2}^{1/2}q(1-q)TT^\top\,dt,\quad
 g_A=\int(q-\pi)tT\,dt,\quad g_T=\int(q-1)tT\,dt.
\]
Full-\(X\) outcome centering makes the residual contribution orthogonal
to the treatment-model score. The heterogeneity contribution gives
\(C=-I^{-1}g\), so fitting weights changes variance to
\[
 V_{\rm fitted}=V_0-g^\top I^{-1}g.
\]
Put \(J(t)=1/q+1/(1-q)\). Direct conditional incoming-load moments give
\[
 V_{0,A}=\tfrac14\int t^2J(t)\,dt+
 \sigma^2\{\tfrac34+\tfrac1{4M}+
                   \tfrac{\alpha_d(M)}{4M^2}\}\int J(t)\,dt,
\]
\[
 V_{0,T}=\int t^2/q\,dt+
 \sigma^2\int\{(1+1/M)/q+\alpha_d(M)/(M^2(1-q))\}\,dt.
\]
Compute these one-dimensional integrals using stats::integrate with
tight tolerance and retain the reported quadrature errors. Symmetry
provides deterministic checks: \(I_{12}=0\), \(g_{A,2}=0\),
\(g_{T,2}=-1/24\), and equal intercept sensitivities. Solves use the
positive definite information matrix; no regularized substitute is used.

## Paired comparisons and initial configuration

Each replication produces six records: PATE/PATT crossed with
known-first-stage, fitted-naive and fitted-adjusted. The latter two share
exactly the same point estimate and realized matching graph. Known and
fitted estimators use the same evaluation data. If nuisance estimation
fails, retain any valid known-first-stage fits and four explicit fitted
failure records.

The bounded pilot consists of three scenarios, all \(n=200,M=3\) and
40 replications: same-sample scores in \(d=2\), independent training in
\(d=1,m=200\), and estimated weights in \(d=3\). A sourceable pilot-table
helper supplies unique fixed seeds and these declarations. Production
sizes and replication counts will be fixed in a later revision informed
by measured pilot runtime and failure capture, before production results.

Use one indexed L'Ecuyer-CMRG stream per scenario/replication and a
separate nextRNGSubStream for independent training. Record complete
evaluation/training seed states, the scenario/base seed, first-stage
parameters and fit errors. No failed replication is replaced.
CLI runs refuse existing result files or companion manifests and record
code/configuration hashes, runtime/session, host/resource metadata, and
the exact record schema. Batching must reproduce the same records apart
from timings.

## Geometry uncertainty and interpretation

For \(d=1\), \(\alpha_1(M)=M^2+M/2\) is exact. For higher dimensions,
accept separately computed geometry objects; do not run integration
inside estimator replications. Every variance formula is affine in
\(\alpha_d(M)\). Propagate its Monte Carlo standard error by multiplying
by the displayed nonnegative alpha coefficient. Missing geometry yields
symbolic coefficients and missing numerical benchmarks. Zero-hit or
zero-empirical-variance integration results remain unresolved, never
silently exact. Do not project an estimate onto the Jensen lower bound.

Report requested/successful/failed counts, bias, empirical and mean fitted
root-n variance, paired adjusted-minus-naive changes, interval coverage
and Monte Carlo uncertainty. The pilot checks calibration, pairing,
failure retention and runtime; 40 replications do not certify coverage.
These experiments assess the stated full-covariate parametric branches,
not unrestricted generated scores or estimated balancing weights.


## Qualified Gaussian same-sample scalar extension

This separately seeded branch uses the score DGP above with d=1, but changes
errors to independent N(0,1/16) and known weights to W=1+Z. The true conditional
means are S1 and S1+1; D=U S1(1-S1), theta0=0. The fitted score is S1+theta D,
and the same unrounded, same-sample OLS coefficient fits both potential means
with known offsets. The event abs(theta)<.8 is checked, never forced by clipping.
A fit outside this event is retained as a failed comparison. The known-theta
comparison is still run. The public wm_gaussian_match constructor performs the
adjustment; its finite linear Gaussian and uniform geometry assumptions are
satisfied by this support-preserving chart. Generic non-Gaussian changing scalar
scores remain outside the declared same-sample theorem.

Conditional on arm, W is constant, so each selected donor receives fraction 1/M.
The unconditional arm probabilities are 1/2, rho=(1,2), gamma=3/2. Direct plus
cross residual second moments are sigma2*(5/2+4), divided by gamma^2, giving
26 sigma2/9. The incoming diagonal and pair terms give
(10 sigma2/9)*(1/M+alpha/M^2). Thus

    V0_A = sigma2 [26/9+(10/9)(1/M+alpha/M^2)],
    V0_T = 2 sigma2 [1+1/M+alpha/M^2].

For both estimands, b=-r/3, Sigma=30 sigma2, and C=10 r sigma2, so the actual
fitted-estimator variance differs from the known-parameter value by
2 b C + b^2 Sigma = -(10/3) r^2 sigma2 = -.075. The naive fitted variance
still targets V0; corrected inference targets V0-.075. This extends the
existing paired six-record schema without changing the earlier three pilots.
A separate 40-replication d1,n200,M3 pilot uses seed 935101. It is diagnostic
and is not pooled with the production cells' independent seeds.
