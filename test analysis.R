#M-estimation example of logistic regression
library("haven")
library("tidyverse")
library("numDeriv")
library("rootSolve")

#Set working directory
orig <- paste0(Sys.getenv("HOME"), "/03 UNC/03 Research Project Materials/10 Zambia PTB/03 Projects/01 HS Interventions/04 Data/02 Clean")
setwd(orig)

#Load original data, filter for no missing data
dat <- read_sas("data.sas7bdat") %>% 
  select(c("outcome_ptb","bp_elev","anemia_pos")) %>%
  filter(!is.na(outcome_ptb)&!is.na(bp_elev)&!is.na(anemia_pos))

#Save dataset
new <- paste0(Sys.getenv("HOME"), "/03 UNC/03 Research Project Materials/06 Cole/10 M Est teaching paper/")
setwd(new)

write.csv(dat, paste0(new,"dat.csv"), row.names=FALSE)
dat <- read.csv(paste0(new,"dat.csv"))
n <- nrow(dat)

### ML estimates using glm
mle <- glm(outcome_ptb ~ anemia_pos + bp_elev, data=dat, family="binomial")
summary(mle)
vcov(mle)

### ML estimates manually
# Fx for log-likelihood
logl <- function(beta){
  p <- plogis(beta[1] + beta[2]*dat$anemia_pos + beta[3]*dat$bp_elev)
  logli <- ifelse(dat$outcome_ptb==1,log(p),log(1-p))
  logl <- sum(logli)
  return(logl)
}

# Fit using optim
optfit <- optim(c(-2,0,0),logl,control=list(fnscale=-1), hessian = TRUE)
optfit$par
hessian <- optfit$hessian
info <- -1*hessian
cov <- solve(info)

# log-likelihood at ML estimates
logl(mle$coefficients)
# first-order partial derivatives (score function)
numDeriv::grad(logl,x=mle$coefficients)
# second-order partial derivatives (Hessian)
numDeriv::hessian(logl,x=mle$coefficients)
# information matrix (first way to define)
infom <- -(numDeriv::hessian(logl,x=mle$coefficients))
# cov matrix
solve(-(numDeriv::hessian(logl,x=mle$coefficients)))

### This is not coming out to the be the same as the information matrix
# information matrix (2nd way to define) - outer product of gradient
grad(logl,x=mle$coefficients) %*% t(grad(logl,x=mle$coefficients))

hessfx <- function(beta){
  p <- plogis(beta[1] + beta[2]*dat$anemia_pos + beta[3]*dat$bp_elev)
  m11 <- -p*(1-p)
  m12 <- -p*(1-p)*dat$anemia_pos
  m13 <- -p*(1-p)*dat$bp_elev
  m21 <- m12
  m22 <- -p*(1-p)*dat$anemia_pos^2
  m23 <- -p*(1-p)*dat$anemia_pos*dat$bp_elev
  m31 <- m13
  m32 <- m23
  m33 <- -p*(1-p)*dat$bp_elev^2
  cbind(m11,m12,m13,m21,m22,m23,m31,m32,m33)
}
matrix(colSums(hessfx(mle$coefficients)),nrow = 3)
matrix(colSums(-hessfx(mle$coefficients)),nrow = 3)
solve(matrix(colSums(-hessfx(mle$coefficients)),nrow = 3))


### M-est

# Estimating functions
score <- function(beta){
  p <- plogis(beta[1] + beta[2]*dat$anemia_pos + beta[3]*dat$bp_elev)
  ef_1 <- (dat$outcome_ptb - p)
  ef_2 <- (dat$outcome_ptb - p)*dat$anemia_pos
  ef_3 <- (dat$outcome_ptb - p)*dat$bp_elev
  return(cbind(ef_1, ef_2, ef_3))
}

# Function sums each column of estimating fx 
sumstack <- function(parms, estfx){ 
  sums <- colSums(estfx(parms))
  return(sums)
}

# Solve for root
gest <- multiroot(f=sumstack, start=c(-2,0,0), estfx=score, atol=1e-12)
gest$root

# Empirical sandwich
gradi <- score(gest$root)
colSums(gradi)/n
meat <- t(gradi) %*% gradi/n
bread <- -numDeriv::jacobian(sumstack, gest$root, estfx=score)/n
# solve(bread)
# n*solve(-numDeriv::jacobian(sumstack, gest$root, estfx=score))
sandwich <- (solve(bread)%*%meat%*%t(solve(bread)))/n
stderr <- sqrt(diag(sandwich))

solve(-numDeriv::jacobian(sumstack, gest$root, estfx=score))%*%(t(gradi) %*% gradi)%*%t(solve(-numDeriv::jacobian(sumstack, gest$root, estfx=score)))
solve(infom) %*% (infom) %*% t(solve(infom))


# Standarization
score <- function(beta){
  p <- plogis(beta[1] + beta[2]*dat$anemia_pos + beta[3]*dat$bp_elev)
  ef_1 <- (dat$outcome_ptb - p)
  ef_2 <- (dat$outcome_ptb - p)*dat$anemia_pos
  ef_3 <- (dat$outcome_ptb - p)*dat$bp_elev
  ef_a1 <- (plogis(beta[1] + beta[2]*1 + beta[3]*dat$bp_elev) - beta[4])
  ef_a0 <- (plogis(beta[1] + beta[2]*0 + beta[3]*dat$bp_elev) - beta[5])
  ef_rd <- as.vector(rep(beta[4] - beta[5], n) - beta[6])
  
  return(cbind(ef_1, ef_2, ef_3,ef_a1,ef_a0,ef_rd))
}

gest <- multiroot(f=sumstack, start=c(-2,0,0,.1,.1,0), estfx=score, atol=1e-12)
gest$root

gradi <- score(gest$root)
meat <- t(gradi) %*% gradi/n
bread <- -numDeriv::jacobian(sumstack, gest$root, estfx=score)/n
(solve(bread)%*%meat%*%t(solve(bread)))/n

solve(-numDeriv::jacobian(sumstack, gest$root, estfx=score))%*%(t(gradi) %*% gradi)%*%t(solve(-numDeriv::jacobian(sumstack, gest$root, estfx=score)))

