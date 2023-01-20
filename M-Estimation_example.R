##############################################
#
# M-estimation: a worked example with connections 
# to maximum likelihood estimation 
# XXX et al.
#
# R Code prepared by XXX (2023/01/20)
#
##############################################

library("tidyverse")
library("numDeriv")
library("rootSolve")
library("geex")

# Create dataset
dat <- tibble(anemia=c(rep(0,4),rep(1,4)),
              bp=rep(c(rep(0,2),rep(1,2)),2),
              ptb=rep(c(0,1),4),
              n=c(496,74,113,25,85,15,15,3)) %>%
  uncount(n)  

# Save nsize
n <- nrow(dat)

##############################
#
# Illustrative example      
#
##############################

##########################
### MAXIMUM LIKELIHOOD ###
##########################

### Using glm()
mle <- glm(ptb ~ anemia + bp, data=dat, family="binomial")
mle$coefficients #parameter estimates
solve(vcov(mle))/n #information matrix - hessian based estimator
vcov(mle) #covariance matrix


### Manually
# A function for the log likelihood
logl <- function(beta){
  p <- plogis(beta[1] + beta[2]*dat$anemia + beta[3]*dat$bp)
  logli <- ifelse(dat$ptb==1,log(p),log(1-p))
  logl <- sum(logli)
  return(logl)
}

# Get MLE using optim() and numerically approximate the Hessian
optfit <- optim(c(-2,0,0),logl,control=list(fnscale=-1), hessian = TRUE)
optfit$par #parameter estimates
-1*optfit$hessian/n #information matrix
solve(n*-1*optfit$hessian/n) #covariance matrix


##########################
###   M-ESTIMATION     ###
##########################

### Manually
# A function for the estimating functions, i.e., score functions
score <- function(beta){
  p <- plogis(beta[1] + beta[2]*dat$anemia + beta[3]*dat$bp)
  ef_1 <- (dat$ptb - p)
  ef_2 <- (dat$ptb - p)*dat$anemia
  ef_3 <- (dat$ptb - p)*dat$bp
  return(cbind(ef_1, ef_2, ef_3)) # Returns an n by 3 matrix
}

# A function to sum each estimating function 
sumstack <- function(parms,estfx){ 
  sums <- colSums(estfx(parms))
  return(sums) # Returns a 1 by 3 row vector
}

# Estimate parameters by finding root (where estimating equations = 0)
mest <- rootSolve::multiroot(f=sumstack, start=c(-2,0,0), estfx=score)
mest$root

# Empirical sandwich
bread <- -numDeriv::jacobian(sumstack, mest$root, estfx=score)/n #-Hessian/n = information matrix
score_i <- score(mest$root) 
meat <- t(score_i) %*% score_i/n #transpose order is different than text bc of matrix orientation
cov_S <- (solve(bread)%*%meat%*%t(solve(bread)))/n
cov_S

### Using geex package
# Define estimating function
geex_ef <- function(data){
  ptb <- data$ptb
  anemia <- data$anemia
  bp <- data$bp
  function(theta){
    p <- plogis(theta[1] + theta[2]*anemia + theta[3]*bp)
    c((ptb - p),
      (ptb - p)*anemia,
      (ptb - p)*bp)
  }
 }

geexresults <- m_estimate(estFUN=geex_ef, data=dat, root_control = setup_root_control(start = c(-2,0,0)))
roots(geexresults)
vcov(geexresults)

##################################################
###   Comparison of MLE and M-estimation       ###
##################################################

# Beta estimates - showing up to 12 decimal places
sprintf("%.12f",mle$coefficients)
sprintf("%.12f",mest$root)

# Compare information matrix with bread and meat
-1*optfit$hessian/n #information matrix - Hessian based estimator
bread

t(score_i) %*% score_i/n #information matrix - outerproduct of score functions estimator
meat

# Cov matrix
vcov(mle) #using built in glm function (Hessian based estimator)
solve(n*-1*optfit$hessian/n) #Using Hessian based estimator
solve(n*t(score_i) %*% score_i/n) #using outer product of score functions estimator
cov_S #empirical sandwich estimator


##############################
#
# Standardization example      
#
##############################

### Manually
# A function for the estimating functions
fx_ef <- function(theta){
  beta <- theta[1:3]
  mu <- theta[4:5]
  delta <- theta[6]
  
  p <- plogis(beta[1] + beta[2]*dat$anemia + beta[3]*dat$bp)
  ef_1 <- (dat$ptb - p)
  ef_2 <- (dat$ptb - p)*dat$anemia
  ef_3 <- (dat$ptb - p)*dat$bp
  ef_r1 <- (plogis(beta[1] + beta[2]*1 + beta[3]*dat$bp) - mu[1])
  ef_r0 <- (plogis(beta[1] + beta[2]*0 + beta[3]*dat$bp) - mu[2])
  ef_rd <- as.vector(rep(mu[1] - mu[2], n) - delta[1])
  
  return(cbind(ef_1,ef_2,ef_3,ef_r1,ef_r0,ef_rd))
}

# Find root, calculated bread, meat, and sandwich cov
mest1 <- multiroot(f=sumstack, start=c(-2,0,0,.1,.1,0), estfx=fx_ef)
bread1 <- -numDeriv::jacobian(sumstack, mest1$root, estfx=fx_ef)/n
meat1 <- t(fx_ef(mest1$root)) %*% fx_ef(mest1$root)/n
cov_S1 <- (solve(bread1)%*%meat1%*%t(solve(bread1)))/n

# Print results
mest1$root
cov_S1
# risk difference with 95% CI
c(mest1$root[6],
  mest1$root[6]-1.96*sqrt(cov_S1[6,6]),
  mest1$root[6]+1.96*sqrt(cov_S1[6,6]))

### Using geex package
# Define estimating function
geex_ef1 <- function(data){
  ptb <- data$ptb
  anemia <- data$anemia
  bp <- data$bp
  
  function(theta){
    beta <- theta[1:3]
    mu <- theta[4:5]
    delta <- theta[6]
    
    p <- plogis(beta[1] + beta[2]*anemia + beta[3]*bp)
    ef_1 <- (ptb - p)
    ef_2 <- (ptb - p)*anemia
    ef_3 <- (ptb - p)*bp
    ef_r1 <- plogis(beta[1] + beta[2]*1 + beta[3]*bp) - mu[1]
    ef_r0 <- plogis(beta[1] + beta[2]*0 + beta[3]*bp) - mu[2]
    ef_rd <- mu[1] - mu[2] - delta[1]
    return(c(ef_1,ef_2,ef_3,ef_r1,ef_r0,ef_rd))
  }
}

geexresults1 <- m_estimate(estFUN=geex_ef1, data=dat, root_control = setup_root_control(start = c(-2,0,0,0.1,0.1,0)))
roots(geexresults1)
vcov(geexresults1)
