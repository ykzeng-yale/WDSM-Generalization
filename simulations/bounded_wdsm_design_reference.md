# Bounded iid counterpart of the original WDSM structural design

This is a new iid selection experiment with the original treatment, outcome and selection equations. It is not the historical Gaussian clustered/fixed-population experiment. Sourcing the file runs no simulation. An executed single dataset is an engineering check; empirical sampling variance and coverage require a separately reviewed replication plan.

The default radius is r=1/4. Draw the first three coordinates of a uniform point on the unit sphere in R4, multiply by r and add (0,0,1), giving (A,L,T). Independently draw a fair sign B, R uniformly on [1,2], and X5,X6 uniformly on [-1/2,1/2]. Set X1=BR, X2=BT/R and solve the two-by-two linear system so a'X=A and b'X=L. The source a,b coefficients, a7=log(1.10), treatment intercepts log(35/80) and log(20/80), and treatment slopes .6 and 2 are preserved. The original setting names describe which source index is used; they do not assert that the new bounded covariate law has the same overlap distribution.

Treatment probability under the population P is expit(a0+delta*(A+a7*T)). Means are mu0=.3L+.75T and mu1=1+.5L+1.05T. Two independent standard normals truncated at +/-2 and then rescaled to variance one give errors epsilon0 and epsilon1; potential outcomes use epsilon0 and epsilon0+epsilon1 respectively. Thus conditional variances remain 1 and 2 and their covariance remains 1. Drawing those independent errors after selection is equivalent and saves work.

Observed iid rows have law Q proportional to pi_Z(X)P, where the original selection equation is expit(log(.005/.995)+log(.9)*Z+c'X). The code derives a uniform bound pi_upper from interval bounds on reconstructed X and the sharper exact decomposition c'X=lambda_A*A+lambda_L*L+(c-lambda_A*a-lambda_L*b)'X, where lambda removes coordinates 3 and 4. Accepting a candidate with probability pi_Z/pi_upper gives exactly this Q. It avoids the unnecessary ~.005 acceptance rate while retaining the original weights W=1/pi_Z; the envelope is not incorporated into the weights. The first n accepted rows are iid Q. The function respects the existing R RNG state and never silently resets the seed. Candidate and numerical budgets fail explicitly.

The target weights are known and may depend on treatment. There are no strata or clusters. Correct feature matrices contain intercept, six X coordinates and X1X2; omitted-interaction matrices contain only intercept and the six X coordinates. PATT's unused treated PG remains the fitting routine's responsibility, not a second target calculation.

## Exact target reduction

The population contrast is 1+.2L+.3T. Symmetry gives E L=0, E T=1 and therefore PATE=1.3 exactly. Let kappa=sqrt(1+a7^2), eta_c=a0+delta*a7 and let x have semicircle density (2/pi)*sqrt(1-x^2) on [-1,1]. Rotational symmetry of the uniform sphere gives

\[
\rho=E_P Z=\int_{-1}^1\operatorname{expit}(\eta_c+\delta r\kappa x)
                 \frac2\pi\sqrt{1-x^2}\,dx,
\qquad
J=\int_{-1}^1 x\operatorname{expit}(\eta_c+\delta r\kappa x)
                 \frac2\pi\sqrt{1-x^2}\,dx.
\]

The conditional mean of the centered three-vector given its projection in direction (1,0,a7)/kappa equals that direction times r*x. Hence E[LZ]=0 and E[TZ]=rho+(r*a7/kappa)*J. The PATT is

\[
E_P[Y(1)-Y(0)\mid Z=1]
=1+.2\frac{E_P[LZ]}\rho+.3\frac{E_P[TZ]}\rho
=1.3+\frac{.3r a_7}{\kappa}\frac J\rho.
\]

This is the population-treated target before selection. Neither selected-sample truth nor historical fixed-population truth is substituted.

## Deterministic integration bound

Set x=cos(theta) to integrate smooth functions on [0,pi]. The rho integrand is (2/pi)sin^2(theta)expit(eta_c+B cos(theta)); J has one additional cos(theta), where B=delta*r*kappa.

Writing p=expit(eta), the absolute coefficient sums of its derivative polynomials of orders 1 through 4 are 2,6,26,150. These bound the derivative magnitudes because 0<=p<=1. The chain rule consequently bounds derivatives of expit(eta_c+B cos(theta)) of orders 0 through 4 by

\[
1,\quad 2|B|,\quad6B^2+2|B|,\quad26|B|^3+18B^2+2|B|,
\quad150B^4+156|B|^3+42B^2+2|B|.
\]

The corresponding derivative bounds of sin^2(theta) are 1,1,2,4,8. Since cos(theta)sin^2(theta)=(cos(theta)-cos(3theta))/4, its order-k derivative is bounded by (1+3^k)/4. Leibniz's rule yields the two explicit fourth-derivative bounds in the code. Composite Simpson with even interval count N and step h=pi/N then has error at most pi*h^4*sup|f''''|/180. The interval count doubles until both integral error bounds satisfy the requested tolerance.

The treatment probability is at least rho_min=expit(eta_c-|B|). If numerical integrals are rho_hat,J_hat with bounds e_rho,e_J, ratio error is at most e_J/rho_min+|J_hat|*e_rho/(rho_min*rho_hat). Multiplying by .3*r*a7/kappa gives the returned PATT bound. All bounds exclude floating-point roundoff; these are deterministic integrations, not Monte Carlo estimates of the targets.

## Inference boundary

The generator does not certify any estimator assumptions and returns that fact explicitly. The original correct-PG CorCor/MisCor population means satisfy full-X centering, but a theorem-controlled study also needs the reviewed support/current-density conditions, regular fitting roots, exact matching and complete fitted variance. The source-form geometry extension and the known-weight full-X pipeline are separate proofs. CorMis/MisMis remain stress cases unless their additional centering/identification conditions are separately proved. No variance or coverage claim follows from the generator or accurate target integration alone.
