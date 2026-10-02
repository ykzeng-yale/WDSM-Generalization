# Observable graph-transport reference

`wm_graph_transport_reference.R` implements the independently reviewed construction in the private research note `feasible_graph_transport_drift.md`. It is standalone; no package exports or automatic model-validity claim are added.

`wm_graph_transport_reference()` takes the actual raw matching chart for one donor direction, positive weights, arm labels, fitted residuals, effective chart tangents of shape n by d by p, and the evaluated values of a justified interior cutoff. It returns the unnormalized incoming graph slope j_z in original raw-weight units. For PATE use (j_1-j_0)/mean(W); for PATT use -j_0/mean(Z*W). Add the complete prediction sensitivity separately, with supplied known weights unchanged. Weight-estimation uncertainty is outside the current WM framework. A coordinate derivative is not necessarily the effective chart tangent: the latter removes the within-chart reparameterization as specified in the theorem.

The implementation uses the C2 product kernel with univariate factor (35/32)(1-u^2)^3 on |u|<1, the bandwidth n^(-1/[2(d+4)]) unless supplied, and an explicit vanishing donor-density threshold 0.01*n^(-1/4). It stops when that threshold is violated at an active donor. It does not silently trim an additional region or substitute a zero transport term. Coordinates, metric transformations, cutoff support and bandwidth units must be specified coherently by the caller.

The conditional donor weight distribution is the empirical kernel-weighted law. A one-dimensional Laplace integral evaluates the exact tagged donor factor for any fixed M, with no kernel in the weight coordinate. This accepts discrete, continuous and mixed positive weight samples. M=1 and a locally constant weight law are evaluated exactly. Other cases use deterministic composite Simpson quadrature with the analytic error bound below; exhausting the node budget raises an error. Temporary matrices are chunked. This is a correctness reference, not an optimized large-n implementation.

## Quadrature derivation

Fix a root weight w and a local donor law with probabilities pi_l, weights w_l and spatial derivatives g_l. Put m=M-1, a=w+m*min(w_l), b=w+m*max(w_l). The tagged factor and its spatial derivative are

\[
A=E\frac{w}{w+\sum_{r=1}^mW_r}
=\int_0^\infty we^{-wt}\Big(\sum_l\pi_le^{-tw_l}\Big)^m\,dt
=\int_0^\infty\frac wa e^{-wu/a}\Big(\sum_l\pi_le^{-uw_l/a}\Big)^m\,du,
\]

\[
\partial_k A
=m\int_0^\infty\frac wa e^{-wu/a}
 \Big(\sum_l\pi_le^{-uw_l/a}\Big)^{m-1}
 \Big(\sum_lg_{lk}e^{-uw_l/a}\Big)\,du.
\]

Each component is a finite signed mixture of exponentials whose rates lie in [1,b/a]. Define C_0=w/a and C_k=m*(w/a)*sum_l|g_lk|. Its absolute integrand is bounded by C_k exp(-u); its fourth u derivative is bounded by C_k*(b/a)^4 exp(-u). The integral above U therefore has absolute value at most C_k exp(-U).

For composite Simpson with panel width 2h and U=2Nh, the per-panel fourth-derivative error is bounded by h^5/90 times that derivative supremum. Summing the geometric envelope gives

\[
|\mathrm{error}_k|
\le C_k e^{-U}
 +\frac{C_k(b/a)^4h^5}{90}
    \sum_{r=0}^{N-1}e^{-2rh}
= C_k e^{-U}
 +\frac{C_k(b/a)^4h^5}{90}
    \frac{1-e^{-U}}{1-e^{-2h}}.
\]

The code chooses U so the tail bound is at most half the requested tolerance and doubles the panel count until the complete bound meets it for A and every gradient component. Evaluation factors out exp(-u), leaving conditional Laplace averages in [0,1], which avoids exponent underflow unnecessarily corrupting the derivative. Common weight scaling leaves these factors unchanged. Accepted roundoff-sized violations of the input probability/gradient sum constraints are normalized by the corresponding quotient rule.

The incoming load k=M*R*A, where R is the weighted query subdensity divided by donor subdensity. Its gradient error is bounded componentwise by M*(|grad R|*error_A+|R|*error_gradA). Absolute residuals, cutoff and chart tangents propagate this to the returned bound for j_z.

These are analytic quadrature bounds for the supplied finite empirical law in exact arithmetic. They exclude floating-point roundoff, kernel bias, sampling error and nuisance error. To apply the consistency theorem, choose tolerances tending to zero; the fixed default is a numerical reference setting, not an asymptotic sequence. Statistical validity additionally requires the reviewed own-(score,weight) centering, signed interior support, ordinary marked-density smoothness, fitted-chart conditions and the point estimator theorem. The returned `assumptions_verified` remains false.
