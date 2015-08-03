sim <- function(N, d, pars, tau_min=0, tau_max=d+2, beta_min=0, beta_max=d+2){

# simulate a data set from model 12.
# N -- number of characters
# d -- number of time periods
# pars -- 6-vector of parameters
# tau_min, tau_max, beta_min, beta_max as in the model.
#
# RETURN VALUE:
# -- pars
# -- M, the simulated data
# -- M.extra for the nex 2 time periods.
# -- log_lambda, tau, beta vectors of parameters

if (!require(truncnorm, quietly=T)){
  cat("\ntruncnorm package not available. Creating homebrew rtruncnorm function in global environment.\n")
  cat("This may possibly cause numerical instabilities.\n\n")
  
  rtruncnorm <- function(N, a, b, mu, sigma, tol=1e-6){

    Phi <- function(x) pnorm(x, mu, sigma)
    Phi.inv <- function(x) qnorm(x, mu, sigma)

    if (abs(Phi(b)-Phi(a)) < tol){
      if (mu > b) return(b)
      if (mu < a) return(a)
    }

    Phi.inv(runif(N)*(Phi(b)-Phi(a))+Phi(a))
  }
  assign("rtruncnorm", rtruncnorm, env=.GlobalEnv)
}

mu_lambda <- pars[1]
sigma_lambda <- pars[2]
mu_tau <- pars[3]
sigma_tau <- pars[4]
mu_beta <- pars[5]
sigma_beta <- pars[6]

log_lambda <- rnorm(N, mu_lambda, sigma_lambda)
tau <- rtruncnorm(N, tau_min, tau_max, mu_tau, sigma_tau)
beta <- rtruncnorm(N, beta_min, beta_max, mu_beta, sigma_beta)

M.plus <- matrix(0, N, d+2)

for (i in 1:N){
  for (t in 1:(d+2)){
    g <- abs(beta[i] - t)
    if (g < tau[i]){
      M.plus[i,t] <- rpois(1, exp(log_lambda[i]))
    }else{
      M.plus[i,t] <- 0
    }
  }
}
M <- M.plus[,1:d]
M.extra <- M.plus[,(d+1):(d+2)]

list(pars=pars, M=M, M.extra=M.extra, log_lambda=log_lambda, tau=tau, beta=beta)
}