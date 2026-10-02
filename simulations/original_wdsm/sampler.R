draw_sample <- function(replication) {
    set.seed(config$sample_seed_base + replication)
    if (design == "retrospective") {
        c0 <- log(config$n_sample/nrow(population)/(1 - config$n_sample/nrow(population)))
        lp <- c0 + log(0.9) * population$A + log(1.05) * population$X1 + 
            log(1.1) * population$X2 + log(1.15) * population$X3 + 
            log(1.1) * population$X4 + log(1.05) * population$X5 + 
            log(1.1) * population$X6
        probability <- exp(lp)/(1 + exp(lp))
        selected <- rbinom(nrow(population), 1, probability) == 
            1L
        sample <- population[selected, ]
        sample$id <- which(selected)
        sample$survey_weight <- 1/probability[selected]
    }
    else {
        pop <- population
        pop$ones <- 1
        pop$id <- seq_len(nrow(pop))
        size3 <- c(rep(850/5, 5), rep(750/5, 5), rep(700/5, 5), 
            rep(650/5, 5), rep(600/5, 5), rep(400/5, 5), rep(350/5, 
                5), rep(300/5, 5), rep(250/5, 5), rep(150/5, 
                5))
        sampled_design <- sampling::mstage(pop, stage = list("stratified", 
            "cluster", ""), varnames = list("Strata", "Cluster", 
            "ones"), size = list(rep(1e+05, 10), rep(5, 10), 
            size3), method = list("", "srswor", "srswor"))
        sample <- sampling::getdata(pop, sampled_design)[[3L]]
        sample$survey_weight <- 1/sample$Prob
    }
    sample
}
