set.seed(280814)

source("got_readin.R")
source("fit_model12.R")
source("simulate_model12.R")
source("plotting_for_model12.R")
source("inference_test_model12.R")

cat("Performing main model fit\n")
fit <- Gibbs12(101000, M.smoothed, burnin=1000, thin=100)

plot_24(fit)
plot_prob0(fit, to.png=T)

cat("Performing validation on M.1\n")
fit.1 <- Gibbs12(101000, M.1, burnin=1000, thin=100)

plot_validation(fit.1, true=M[1:9, 3], to.png=T)

cat("Performing inference test\n")

OUT <- inference_test(N.sims=100, pars=pars)

plot_coverage_png(OUT)

put_matrix(fit$pred1, char.names)