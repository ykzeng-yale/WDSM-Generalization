N = 1e+06
n = N * 0.005
n_covar = 7
n_strata = 10
n_cluster = 20
n_scenarios = 1
tau_strata = 0.35
min_prob = 0
tau_cluster = c(0.15)
s.m = 1
set.seed(20220930)
data_sets = list()
for (scenario in 1:n_scenarios) {
    mu_strata = matrix(rnorm(n_strata * n_covar, mean = 0, sd = tau_strata), 
        nrow = n_strata, ncol = n_covar)
    mu_cluster = matrix(rnorm(n_cluster * n_covar, mean = 0, 
        sd = tau_cluster[scenario]), nrow = n_cluster, ncol = n_covar)
    strata_row = rep(NA, n_strata)
    for (strata in 1:n_strata) {
        strata_row[strata] = paste("Strata", strata, sep = "_")
    }
    rownames(mu_strata) = strata_row
    strata_col = rep(NA, n_covar)
    for (covariate in 1:n_covar) {
        strata_col[covariate] = paste("Covariate", covariate, 
            sep = "_")
    }
    colnames(mu_strata) = strata_col
    cluster_row = rep(NA, n_cluster)
    for (cluster in 1:n_cluster) {
        cluster_row[cluster] = paste("Cluster", cluster, sep = "_")
    }
    rownames(mu_cluster) = cluster_row
    cluster_col = rep(NA, n_covar)
    for (covariate in 1:n_covar) {
        cluster_col[covariate] = paste("Covariate", covariate, 
            sep = "_")
    }
    colnames(mu_cluster) = cluster_col
    cluster_l = list()
    strata_l = list()
    variable_l = list()
    for (variable in 1:n_covar) {
        for (strata in 1:n_strata) {
            for (cluster in 1:n_cluster) {
                data = cbind(rep((strata - 1) * n_cluster + cluster, 
                  N/(n_strata * n_cluster)), rep(strata, N/(n_strata * 
                  n_cluster)), rnorm(N/(n_strata * n_cluster), 
                  mean = mu_cluster[cluster, variable] + mu_strata[strata, 
                    variable]))
                colnames(data) = c("Cluster", "Strata", paste("X", 
                  variable, sep = ""))
                cluster_l[[cluster]] = data
            }
            strata_l[[strata]] = do.call("rbind", cluster_l)
        }
        variable_l[[variable]] = data.frame(do.call("rbind", 
            strata_l))
    }
    data = variable_l[[1]]
    for (variable in 2:n_covar) {
        data = cbind(data, variable_l[[variable]][, 3])
    }
    var_col = rep(NA, n_covar)
    for (covariate in 1:n_covar) {
        var_col[covariate] = paste("x", covariate, sep = "")
    }
    colnames(data) = c("Cluster", "Strata", var_col)
    expit = function(x) {
        exp(x)/(1 + exp(x))
    }
    a0 = log(35/80)
    a1 = log(1.1)
    a2 = log(1.25)
    a3 = log(1.5)
    a4 = log(1.75)
    a5 = log(2)
    a6 = log(2.5)
    a7 = log(1.1)
    delta_z = 0.6
    formula = a0 + delta_z * (a1 * data$x1 + a2 * data$x2 + a3 * 
        data$x3 + a4 * data$x4 + a5 * data$x5 + a6 * data$x6 + 
        a7 * data$x1 * data$x2)
    prob_treat = expit(formula)
    data$prob_treat = prob_treat
    data$z = rbinom(N, 1, prob_treat)
    data$TreatmentGroup <- factor(data$z, levels = c(0, 1), labels = c("Control (z=0)", 
        "Treatment (z=1)"))
    b0 = 0
    b1 = 2.5
    b2 = -2
    b3 = 1.75
    b4 = -1.25
    b5 = 1.5
    b6 = 1.1
    b70 = 2.5
    b71 = 1.5
    delta0 = 0.3
    delta1 = 1
    delta2 = 0.2
    data$y0 = b0 + delta0 * (b1 * data$x1 + b2 * data$x2 + b3 * 
        data$x3 + b4 * data$x4 + b5 * data$x5 + b6 * data$x6 + 
        b70 * data$x1 * data$x2) + rnorm(N)
    data$y1 = data$y0 + delta1 + delta2 * (b1 * data$x1 + b2 * 
        data$x2 + b3 * data$x3 + b4 * data$x4 + b5 * data$x5 + 
        b6 * data$x6 + b71 * data$x1 * data$x2) + rnorm(N)
    data$y = data$z * data$y1 + (1 - data$z) * data$y0
    data_sets[[scenario]] = data
}
SuperPopGen_GoodOverlap_AllinOne_0118_2025 <- data_sets
