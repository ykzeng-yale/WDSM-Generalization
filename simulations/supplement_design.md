# Pre-specified supplemental validation

Status: specification and source implementation, before any supplemental pilot.
These experiments supplement the primary strong/weak marked designs. They do
not change the primary grid, replace failed observations, or establish new theory.
All weights and matching maps here are known. Outcomes and errors are bounded.

## Fixed grid and reproducibility

The eight cells are spline (d,n)=(1,800),(1,2000),(2,800),(2,2000), strata
n=800,2000, and independent-block potential means n=2000 with potential-one
block fractions 1/2,1/3. Every cell uses M=3. The bounded pilot has five
replications per cell; the eventual specification is 200, conditional on review
of correctness and runtime. Seeds 933101 through 933108 index these cells in
this order. Replication r uses the rth L'Ecuyer-CMRG stream of its cell seed,
with Inversion/Rejection RNG conventions. Pilot observations remain IDs 1--5;
production continuation starts at 6 if the specification is unchanged.

There are eight records per spline/stratum replication (two methods, PATE/PATT,
oracle/feasible) and three records per split replication (the PATE and its two
component potential means). Thus the pilot requests 270 records, and 200
replications per cell request 10,800. The component records reuse one split fit.
Companion RDS files retain seeds, nuisance fits, training row IDs, ranks,
clipping counts, component sizes, errors and provenance. CSV records include
every requested comparison, including failed fits. Failed donors, rank tests,
nonfinite quantities and numerical guards are not repaired by resampling,
ridge, dropping terms, changing the target or retuning the specification.

## Shared bounded mark structure

Let U be a fair Bernoulli mark, independent of the continuous uniform scores
and of treatment conditional on any stratum. W0=.5+U and W1=1+2U, with
observed W=WZ. Potential errors have independent fair signs conditional on U,
with conditional variances v0(U)=(.25,1) and v1(U)=(.5,1.5). Thus residuals
are centered conditional on all score, stratum and weight marks. The potential
outcome means specified below are the ordinary and weighted arm means alike.
The mean weights are rho=(1,2), second weight moments nu=(5/4,5), and
E(Wz² epsilon_z²)=(37/32,7). Consequently the stabilized residual mark
eta_z=Wz epsilon_z/rho_z has variance t=(37/32,7/4) and covariance
E[(Wz-rho_z)eta_z]=0. Treatment is independent of the potential errors and
their weight mark given the stratum. This supplies an actual causal target,
as well as the algebraic weighted-mean identification conditions.

## Nonlinear feasible spline experiment

S is uniform on [0,1]^d and Z is Bernoulli(.4), independently. Set

    mu0(S) = 1 + sin(2*pi*S1) + (.3/d) sum_j cos(2*pi*Sj),
    mu1(S) = mu0(S) + 1 + S1.

Both weighted effects are 1.5; Var(mu1-mu0)=1/12. The means are analytic and
bounded, and are not correctly specified by a fixed complete quadratic basis.
Compare oracle correction with wm_fit regression='spline', order r=3,
mesh exponent a=1/6, moment_order=Inf, and supplied support [0,1]^d.
The conditional weight mean is fitted by unweighted least squares with the
same basis and clipped at the predeclared interval (.25,4). The true rho values
are strictly interior. Outcome regressions use individual W as WLS weights.

Use the fixed alternating two-fold assignment 1,2,1,2,... . Both methods fit
on the complementary fold and match wholly inside the evaluation fold.
Oracle comparisons use exactly these same within-fold graphs. For d=1,2,
a=1/6 lies strictly inside the bounded-moment spline windows, including the
stabilized requirements a>1/(4r) and a<1/(2d). Training sample size, rather
than the total sample size, determines the implemented mesh. The records retain
the resulting basis size/rank and fold-specific donor counts.

The effect and all residual/weight laws equal the primary strong design, so
wm_benchmark_strong with q=.4 gives both method-specific variance limits.
Pooling two fixed positive-proportion evaluation folds does not multiply the
root-n limit by two: their contributions sum as sum_l lambda_l V=V after
negligible fitted-nuisance remainders. The oracle and feasible estimators share
this limit, but need not have the same finite-sample point estimate.

The self-normalized d=2 benchmark needs the complete overlap vector beta;
alpha alone is insufficient with these nonconstant donor weights. An alpha-only
geometry object may supply the stabilized benchmark, leaving the original-rule
benchmark explicitly unavailable. Numeric geometry errors retain their status
and covariance; no low-quality geometric estimate is treated as exact.

## Finite random strata

Let L have probabilities pi=(13/20,7/20), Z|L have q=(3/10,3/5), and
S~Unif[0,1], independently of L and U. Set mu0=1+S+S²/2 and
mu1=mu0+c_L+S, with c=(.5,1.5). The local effect means are tau_l=(1,2)
and their within-stratum variance is 1/12. Matching and nuisance fits remain
inside the exact stratum. Compare oracle and correctly specified complete
degree-two polynomial fits, with the same fixed rho bounds (.25,4).

The arm and population denominators are

    gamma_l = 1+q_l = (13/10,8/5),       gamma = 281/200,
    gammaT_l = 2*q_l = (3/5,6/5),       gammaT = 81/100.

Hence the global PATE is 393/281 and the global PATT is 41/27. The global
targets are not the unweighted stratum average and do not condition on the
realized stratum counts. Put nu_l=5/4+15*q_l/4=(19/8,7/2) and
nuT_l=5*q_l=(3/2,3). If V_l is the local primary strong-design root-sample-size
variance at q_l, the correct global limits are

    V_A = sum_l pi_l [gamma_l² V_l,A + nu_l(tau_l-tau_A)²] / gamma²,
    V_T = sum_l pi_l [gammaT_l² V_l,T + nuT_l(tau_l-tau_T)²] / gammaT².

To see the offset term, the global unnormalized contribution is the local
centered contribution plus W(tau_l-tau_A), or ZW(tau_l-tau_T). Its local
cross moment vanishes by independent centered S and conditional residual
centering. The extra second moment is nu_l delta_l², which decomposes into
gamma_l² delta_l² (random stratum composition) plus
(nu_l-gamma_l²)delta_l² (within-stratum weight variation). Using only the
first piece drops a real variance term. The same argument applies to treated
weights. Dimension-one beta is exact: beta_l=2(l+1) for l<M-1 and
beta_(M-1)=3M/2, with alpha=M²+M/2 and zero geometric MCSE.

## Independent blocks with unequal score dimensions

Use independent uniform S=(S1,S2), q=.4 and the shared marks. Set
g0=1+S1, g1=2+S1+S2. Potential means are theta0=1.5 and theta1=3,
with variances 1/12 and 1/6. Their target contrast is 1.5. The control map
is scalar S1 and the treated map is two-dimensional S. All means and rho
values are supplied as known functions. Assign the first floor(n*lambda)
observations to block 1, the rest to block 0, independent of all data. Block z
estimates only potential mean z using wm_potential through wm_split_pate.
The two estimates therefore use independent observations.

Here gamma=7/5 and nu=11/4. Let q0=3/5,q1=2/5 and r_z=q_(1-z)/q_z.
The one-direction root-block-size variance is

    V_z = { nu Var(g_z) + q_z t_z [rho_z² + 2 r_z rho_(1-z) rho_z
            + r_z nu_(1-z)/M + alpha_(d_z) r_z² rho_(1-z)²/M²] } / gamma².

This follows by conditioning on the root score/arm and taking the incoming
load's first and second moments. The residual/heterogeneity cross term is zero;
there is no cross-direction reciprocal covariance because the blocks are
independent. Substituting the declared marks gives

    V0 = (25/49)[11/48 + 407/160 + 37/(16M) + 37 alpha1/(30M²)],
    V1 = (25/49)[11/24 + 7 + 21/(16M) + 63 alpha2/(40M²)].

For N1=floor(n*lambda), N0=n-N1, the total-sample root variance benchmark is
(n/N1)V1+(n/N0)V0; its limiting form is V1/lambda+V0/(1-lambda).
Both nominal and realized fractions are recorded. Component CSV rows use their
own analysis_n=Nz and variance Vz/Nz; the contrast uses analysis_n=n.
The d=1 geometric term is exact; only alpha2 is needed for potential one.
Its MCSE enters the contrast with coefficient
(n/N1)(25/49)63/(40M²). The components are reported for diagnosis, not as
additional independent simulated datasets or a full-sample efficiency claim.

## Assessment and failure reporting

The pilot establishes executability, donor/rank adequacy and runtime only.
Eventual summaries should report empirical bias/RMSE, empirical root variance,
mean estimated root variance, Wald coverage, clipping and failure frequencies,
plus geometric uncertainty/status. Use all requested replications as the
failure-rate denominator and clearly identify the successful subset for
estimator summaries. Paired oracle/feasible differences must retain pairing.
No pilot-based selection of signs, outcomes or favorable benchmark estimates
is allowed. All theoretical limits above are asymptotic, not exact finite-n
coverage or exact finite-n variance assertions.
