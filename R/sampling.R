library(MASS)

# Given opinion poll data and modelling parameters, produce information for drawing random samples
# in the multi logit model

election_model <- function(data, covobs_weight, covsamp_weight, covdiag_weight, df) {
  
  # Multi logit tranform, each row now sums to zero
  tdata <- log(data) - apply(log(data),1,mean)
  np <- ncol(tdata)
  
  em <- list()
  em$df <- df
  
  em$m <- apply(tdata,2,mean)
  
  covobs <- cov(tdata)
  
  p <- exp(em$m)/sum(exp(em$m))
  q <- rep(1/np,times = np)
  pinv <- 1/p
  
  # Sampling covariance as it looks in the transformed space. I must have a derivation somewhere...
  covsamp <- diag(pinv) - q %*% t(pinv) - pinv %*% t(q) + sum(pinv) * (q %*% t(q)) 
  
  # Diagonal noise but with correction to enforce that sampled number sum to zero 
  covdiag <- diag(1, np) - matrix(1/np, nrow=np, ncol=np)
  
  # Add covariances, and put a little extra on the diag for numerical stability
  em$c <- covobs_weight*covobs +
    covsamp_weight*covsamp +
    covdiag_weight*covdiag + 
    diag(1e-12, np)
  
  em
}


# Draw samples from the multi logit model

r_election_model <- function(emodel, ns, seed = NULL) {
  
  if (is.null(seed)) set.seed(123)
  else if (is.integer(seed)) set.seed(seed)
  
  np <- length(emodel$m)
  
  # Draw normal variates for variables * replication
  rr <- matrix(rnorm(np*ns),nrow = np)
  
  # Modify normal variates by a random factor for each replication. Makes them each t distributed
  r <- rr %*% diag((rchisq(ns,emodel$df)/emodel$df)^-.5)
  
  # Transform using covariance matrix and mean
  mlres <- t(chol(emodel$c)) %*% r + emodel$m
  
  # Transform back to probabilities
  apply(mlres,2,function(x) exp(x)/sum(exp(x)))
}

# Compute the chi^2 of result given the model

chi2_election_model <- function(emodel, result) {
  
  tresult <- log(result) - mean(log(result))
  
  diffvec <- tresult - emodel$m
  
  c(t(diffvec) %*% MASS::ginv(emodel$c) %*% diffvec)
}

# Compute p value for election result given the model

p_election_model <- function(emodel, result) {
  
  chi2 <- chi2_election_model(emodel, result)
  np <- length(result)
  
  pf(chi2/(np-1),np-1,emodel$df,lower.tail = FALSE)
}