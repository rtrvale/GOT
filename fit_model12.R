
# Model:
# X_it=Pois(lambda[i]) if g(i,t) > 0 or 0 otherwise where
# g(i,t) = tau[i] - abs(beta[i] - t) 
# log(lambda[i]) ~ N(mu_lambda, sigma_lambda^2)
# tau[i] ~ N(mu_tau, sigma_tau^2)*[tau_min, tau_max]
# beta[i] ~ N(mu_beta, sigma_beta^2)*[beta_min, beta_max]


Gibbs12 <- function(
its, M, 
burnin = 0,
thin = 1,
len = 99,
step = rep(3, N),
min_prob = 1e-6,
trace = 0,
priors = NULL,
start = NULL,
beta_min=0, beta_max = ncol(M)+2,
tau_min=0, tau_max = ncol(M)+2
){

# Performs the Gibbs sampling
#
# arguments:
# ==========
# + not optional:
# -- its = number iterations to output
# -- M = input matrix - must have at least 2 columns and 1 row.
#
# + optional:
# -- burnin = size of burnin to discard
# -- thin = size of thin (must divide its-burnin)
# -- len = length of vector of proposed lambdas (betas) in lambda (beta) sampler (smaller is faster)
# -- step = size of step in samples for individual lambdas
# -- min_prob = smallest possible prob in parameter samplers (beta and lambda)
# -- trace = 0 for no output, 1 for some, 2 for some + plots, 3 for some +
#    plots + pause after every plot
# -- priors and start, lists 
# -- beta_min, beta_max, tau_min, tau_max, fixed model parameters which can be changed if desired
#
# return value:
# =============
# a list with 6 components:
# -- pars: sampled values of pars, matrix with 6 columns
# -- log_lambda, matrix with N columns
# -- tau ,       "
# -- beta,       "
# -- pred1, integer matrix with N columns; pred values for next period
# -- pred2, "					 pred values for period after next

# PACKAGES
#
# check for truncnorm package and use something else if it doesn't exist
#
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
  assign("rtuncnorm", rtruncnorm, env=.GlobalEnv)
}

M <- as.matrix(M)
N <- nrow(M)
d <- ncol(M)

if (N < 2) stop("Too few rows in data matrix.")
if (d < 2) stop("Too few columns in data matrix.")

# STARTING VALUES AND PRIORS
#
# Use random start values if none specified.
#
if (is.null(start)){
  start <- list()
  start$log_lambda  <- rnorm(N)
  start$tau <- runif(N)*(tau_max - tau_min) + tau_min
  start$beta <- runif(N)*(beta_max - beta_min) + beta_min
  start$pars <- runif(6)
}
# Use vague priors if none specified
#
if (is.null(priors)){
  priors <- list()
  priors$mu_lambda <- c(0,1000)
  priors$sigma_lambda <- c(0.001,0.001)
  priors$mu_tau <- c(0, 1000)
  priors$sigma_tau <- c(0.001, 0.001)
  priors$mu_beta <- c(0, 1000)
  priors$sigma_beta <- c(0.001, 0.001)
}

# SET UP THE OUTPUT
#
# (note: names of these matrices begin with "out_" so that they don't clash with parameter names)
#
out_pars <- matrix(0, its, 6)
out_log_lambda <- out_tau <- out_beta <- pred1 <- pred2 <- matrix(0, its, N)
colnames(out_pars) <- c("mu_lambda", "sigma_lambda", "mu_tau", "sigma_tau", "mu_beta", "sigma_beta")
colnames(out_log_lambda) <- paste("log_lambda", 1:N, sep="")
colnames(out_tau) <- paste("tau", 1:N, sep="")
colnames(out_beta) <- paste("beta", 1:N, sep="")

params <- start

# MAIN LOOP
#
for (it in 1:its){
  if (trace) cat("\niteration", it, ":\n\n")

  # SAMPLE HYPERPARAMETERS (pars)
  #  
  params$pars <- sample_pars(params, priors, trace) 

  # SAMPLE FOR EACH CHARACTER IN 1:N
  #
  for (i in 1:N){
    if (trace) cat("\n\ncharacter", i, ":\n\n")
    
    x <- M[i,]
    
    # SAMPLE log_lambda[i]
    #
    log_lambda_samp <- sample_log_lambda(x=x, 
                       tau = params$tau[i], 
                       beta = params$beta[i],
                       mu_lambda = params$pars[1],
                       sigma_lambda = params$pars[2],
                       min_prob=min_prob,
                       step = step[i], len=len, trace=trace, 
                       tau_min, tau_max, beta_min, beta_max)
    params$log_lambda[i] <- log_lambda_samp$samp
    step[i] <- log_lambda_samp$step

    # SAMPLE tau[i]
    #
    params$tau[i] <- sample_tau(x, log_lambda=params$log_lambda[i], 
                         beta=params$beta[i], mu_tau=params$pars[3],
                         sigma_tau=params$pars[4], trace=trace, 
                         tau_min, tau_max, beta_min, beta_max, min_prob)

    # SAMPLE beta[i]
    #
    params$beta[i] <- sample_beta(x, log_lambda=params$log_lambda[i], 
                         tau=params$tau[i], mu_beta=params$pars[5],
                         sigma_beta=params$pars[6], trace=trace,
                         tau_min, tau_max, beta_min, beta_max, len, min_prob)

    # MAKE PREDICTIONS
    #
    g_i1 <-  params$tau[i] - abs(d+1 - params$beta[i])
    pred1[it, i] <- ifelse(g_i1 > 0, rpois(1, exp(params$log_lambda[i])), 0)
    g_i2 <-  params$tau[i] - abs(d+2 - params$beta[i])
    pred2[it, i] <- ifelse(g_i2 > 0, rpois(1, exp(params$log_lambda[i])), 0)


  } # end for (i in 1:N)

out_pars[it,] <- params$pars
out_log_lambda[it,] <- params$log_lambda
out_tau[it,] <- params$tau
out_beta[it,] <- params$beta

} # end for (it in 1:its)

# DELETE BURNIN AND THIN
#
outputs <- c("out_pars", "out_log_lambda", "out_tau", "out_beta", "pred1", "pred2")
for (k in 1:6){
  output <- get(outputs[k])
  output <- output[(burnin+1):nrow(output),]
  output <- output[thin*(1:round(nrow(output)/thin)),]
  assign(outputs[k], output)
}

# RETURN
#
list(pars=out_pars, log_lambda=out_log_lambda, tau=out_tau, beta=out_beta, 
pred1=pred1, pred2=pred2)
} # end Gibbs12 function

## FUNCTIONS FOR SAMPLING INDIVIDUAL PARTS
##########################################

sample_pars <- function(params, priors, trace){

# function for sampling hyperparameters
#
# arguments:
# ----------
# params - list of current params
# priors - list of priors
# trace - controls output to console
#
# output:
# -------
# vector of parameters, length 6.

pars <- params$pars

samp_hyp_lambda <- N_sample(params$log_lambda, pars[1], pars[2], 
                   priors$mu_lambda, priors$sigma_lambda, trace)
samp_hyp_tau    <- N_sample(params$tau, pars[3], pars[4], 
                   priors$mu_tau, priors$sigma_tau, trace)
samp_hyp_beta   <- N_sample(params$beta, pars[5], pars[6], 
                   priors$mu_beta, priors$sigma_beta, trace)

if (trace){
  cat("sampling of hyperparameters complete.", "\n")
  cat("old: ", pars, "\n")
  cat("new: ", c(samp_hyp_lambda, samp_hyp_tau, samp_hyp_beta), "\n")
}

c(samp_hyp_lambda, samp_hyp_tau, samp_hyp_beta)
} # end sample_pars function

N_sample <- function(x, mu, sigma, prior.normal, prior.gamma, trace=0){

# sample mu and sigma from a normal distirbution using the usual theory
#
# arguments:
# ----------
# x - data vector
# mu - mean
# sigma - sd
# prior.normal, a 2-vector (prior.mu, prior.sigma)
# prior.gamma, a 2-vector (prior.a, prior.b)
#
# output:
# -------
# a 2-vector containing the sampled values of mu and sigma

N <- length(x)
prior.mu <- prior.normal[1]
prior.sigma <- prior.normal[2]
prior.a <- prior.gamma[1]
prior.b <- prior.gamma[2]

if (trace >= 2) cat("sampling from normal distribution:", "\n")

v <- 1/(prior.sigma^-2 + N*sigma^-2)
m <- v * (prior.mu/prior.sigma^2 + sum(x)/sigma^2)
sampled_mu <- rnorm(1, m, sqrt(v))

if (trace >= 2){
  cat("sd:", sqrt(v), "computed from\n")
  cat("- prior sd:", prior.sigma, "\n")
  cat("- current sd:", sigma, "\n") 
  cat("- N:", N, "\n")
  cat("mean:", m, "\n")
  cat("previous value of mean:", mu, "\n")
  cat("sampled value of mean:", sampled_mu, "\n")
}

if (trace >= 2) cat("sampling from gamma distribution:", "\n")

a <- prior.a + N/2
b <- prior.b + sum((x-mu)^2)/2
sampled_sigma <- sqrt(1/rgamma(1, a, b))

if (trace >= 2){
  cat("a:", a, "\n")
  cat("b:", b, "\n")
  cat("previous value of sigma:", sigma, "\n")
  cat("sampled value of sigma:", sampled_sigma, "\n")
}

c(sampled_mu, sampled_sigma)
} # end N_sample function

###### sample_tau ######
########################

sample_tau <- function(x, log_lambda, beta, mu_tau, sigma_tau, trace=0, tau_min, tau_max, beta_min, beta_max, min_prob){

if (trace >= 3){
  cat("\nSampling tau with log_lambda = ", log_lambda, "\nbeta = ", beta, "\nmu_tau = ",mu_tau,
"\nsigma_tau = ", sigma_tau, "\n")
}

d <- length(x)
breaks <- abs(beta - 1:d)
a0 <- ifelse(sum(x)==0, tau_min, max(breaks[x!=0]))
if (trace) cat( "Value of a0: ", a0, "\n")
a0 <- max(a0, tau_min)
if (a0 >= tau_max){
  stop("Error in sample_tau. a0 is ", a0, "but tau_max is", tau_max, ".\n")
}

if (all(x!=0)) return(rtruncnorm(1, a0, tau_max, mu_tau, sigma_tau))

a <- sort(breaks[x==0])

if (trace) cat("a = ", a, "\n")

if (max(a) <= a0){
  if (trace) cat(a0,"\n", tau_max,"\n", mu_tau,"\n", sigma_tau)
  return(rtruncnorm(1, a0, tau_max, mu_tau, sigma_tau))
}
lambda <- exp(log_lambda)
a <- c(a0, a[(a > a0) & (a < tau_max)], tau_max)
L <- length(a)
areas <- rep(0, L-1)
for (i in 1:(L-1)){
  areas[i] <- pnorm(a[i+1], mu_tau, sigma_tau) - pnorm(a[i], mu_tau, sigma_tau)
  areas[i] <- areas[i]*exp(-lambda*i)
}

if (sum(areas) <= 0){ # this can happen if mu_tau is very large
  areas <- rep(min_prob, length(areas))
  if (trace) cat("Warning. Uniform distribution in sample_tau with mu_tau =", mu_tau, "and sigma_tau=", sigma_tau, "\n")
}

wedge <- sample(1:(L-1), 1, prob=areas)

sampled_tau <- rtruncnorm(1, a[wedge], a[wedge+1], mu_tau, sigma_tau)

if (trace) cat("Sampled value of tau: ", sampled_tau, "\n")

sampled_tau
} # end sample_tau function 

###### sample_beta ######
#########################

sample_beta <- function(x, log_lambda, tau, mu_beta, sigma_beta, trace=0, tau_min, tau_max, beta_min, beta_max,
                        len, min_prob){

if (trace) cat("\nSampling beta\n")
if (trace >= 3){
  cat("log_lambda = ", log_lambda, "\ntau = ", tau, "\nmu_beta = ",mu_beta,
"\nsigma_beta = ", sigma_beta, "\n")
}

d <- length(x)
beta_vect <- seq(beta_min, beta_max, len=len)
lambda <- exp(log_lambda)

if (trace){
  cat("lambda = ",lambda,"\n")
  cat("tau = ",tau,"\n")
}

g <- abs(matrix(rep(1:d, len), nc=d, byrow=T) - beta_vect) < tau
L <- apply(g, 1, function(v) prod(v[x!=0])*prod(exp(-lambda)*v[x==0] + (1-v[x==0])) )
L <- L* exp(-0.5*(beta_vect - mu_beta)^2/sigma_beta^2)
L[L <= 0] <- min_prob

sampled_beta <- sample(beta_vect, 1, prob=L) + (runif(1)-0.5)*(beta_vect[2]-beta_vect[1])

if (trace) cat("\nSampled value of beta: ", sampled_beta, "\n")

sampled_beta
} # end sample_beta function

##### sample_log_lambda ######
##############################

sample_log_lambda <- function(x, tau, beta, mu_lambda, sigma_lambda,
                     min_prob, step, len, trace, tau_min, tau_max, beta_min, beta_max){

d <- length(x)
log_lambda_vect <- seq(-step, step, len=len)

if (trace) cat("\nSampling log_lambda with tau :", tau, " beta: ", beta, "\n")

g <- tau - abs(1:d - beta)
lambda_vect <- exp(log_lambda_vect)
LL <- -sum(g>0)*lambda_vect + sum(x*(g>0))*log_lambda_vect
LL <- LL - 0.5*(log_lambda_vect - mu_lambda)^2/sigma_lambda^2
probs <- exp(LL-max(LL))

if (trace >= 3) plot(log_lambda_vect, probs, "l")

if (probs[1] > min_prob|probs[len] > min_prob){
  step <- step+1
  if (trace) cat("\nStep for lambda_samp updated to", step, "\n")
}
samp <- sample(log_lambda_vect, 1, prob=probs) + 
        (runif(1)-0.5)*(log_lambda_vect[2]-log_lambda_vect[1])

if (trace) cat("Sampled value of log_lambda: ", samp, "\n", "(sampled value of lambda", exp(samp), ")\n")

list(samp=samp, step=step)
} # end sample_log_lambda function