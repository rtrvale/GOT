pars <- c(1.3, 0.75, 2, 1, 4, 1.5)

source("fit_model12.R")
source("simulate_model12.R")

# note: dumping the output partway through makes it a lot slower. OUT is an 
# object of size 18 MB
#

inference_test <- function(N.sims=100, pars, sim.len=1100, burnin=100, thin=10,
N=24, d=5, sd.loc=0.1, sd.log_scale=0.01, trace=1, omit.zero.rows=F
){

# This function simulates data from the model N.sims times with parameters which are
# slightly modified versions of pars.
#
# -- N.sims number of simulations
# -- pars vector of 6 parameters
# -- sim.len length of each simulation
# -- burnin length of burnin for each simulation
# -- thin thinning for each simulation (as in Gibbs12 function)
# -- N, d size of sim$M
# -- sd.loc, sd.log_scale, size of changes in location and log scale parameters
#    in the data produced for the simulations
#
#
# -- OUT: a list of N.sims lists, each one with a component sim, the simulated
#    data object and fit, the fitted model object

OUT <- list()

for (i in 1:N.sims){
  if (trace) print(i)
  pars.new <- pars
  pars.new[c(1,3,5)] <- pars[c(1,3,5)] + rnorm(3, 0, sd.loc)
  pars.new[c(2,4,6)] <- pars[c(2,4,6)]*exp(rnorm(3, 0, sd.log_scale))
  sim.dat <- sim(N, d, pars.new)

  if (omit.zero.rows) sim.dat$M <- sim.dat$M[rowSums(sim.dat$M)!=0, ]  

  fit <- Gibbs12(its=sim.len, sim.dat$M, burnin=burnin, thin=thin, trace=0)
  OUT[[i]] <- list(sim=sim.dat, fit=fit)
}

OUT
} # end inference_test function

### Helper function for quantiles of columns of a matrix ###

colQ <- function(X, alpha) apply(X, 2, function(x) quantile(x, alpha))

count_pars_hits <- function(OUT, alpha){

# Count the coverage in an object returned by inference_test
# returns a single number, the coverage of central alpha-level credible intervals
# for the parameters over all model fits

counter <- function(fitlist){
  sim <- fitlist$sim
  pars <- fitlist$fit$pars
  lower <- colQ(pars, (1-alpha)/2)
  upper <- colQ(pars, 1-(1-alpha)/2)
  true <- sim$pars
  count <- sum(true >= lower & true <= upper)
}

hits <- unlist(lapply(OUT, counter))
avg.hits <- mean(hits)/ncol(OUT[[1]]$fit$pars)
avg.hits
}

count_pred1_hits <- function(OUT, alpha){

# same as count_pars_hits but for pred1

counter <- function(fitlist){
  sim <- fitlist$sim
  pred1 <- fitlist$fit$pred1
  lower <- colQ(pred1, (1-alpha)/2)
  upper <- colQ(pred1, 1-(1-alpha)/2)
  true <- sim$M.extra[,1]
  count <- sum(true >= lower & true <= upper)
}

hits <- unlist(lapply(OUT, counter))
avg.hits <- mean(hits)/ncol(OUT[[1]]$fit$pred1)
avg.hits
}

count_pred1_hits_logistic <- function(OUT, alpha){

# same as count_pred1_hits but also records actual coverage of the calculated intervals

hits <- c()
actual <- c()

counter <- function(fitlist){
  sim <- fitlist$sim
  pred1 <- fitlist$fit$pred1
  lower <- colQ(pred1, (1-alpha)/2)
  upper <- colQ(pred1, 1-(1-alpha)/2)
  true <- sim$M.extra[,1]
  actual.hits <- (true >= lower & true <= upper)
  actual.coverage <- rep(0, ncol(pred1))
  for (k in 1:ncol(pred1)){
    actual.coverage[k] <- sum(pred1[,k] >= lower[k] & pred1[,k] <= upper[k])/
nrow(pred1)
  }
 list(actual.hits=actual.hits, actual.coverage=actual.coverage)
} 

for (i in 1:length(OUT)){
   counted <- counter(OUT[[i]])
   hits <- c(hits, counted$actual.hits)
   actual <- c(actual, counted$actual.coverage)
}

list(hits=hits, actual=actual)
}