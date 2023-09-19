###################################################################################################################
# M-estimation: introduction and applied examples
# 
# Rachael Ross (2023/09/15)
###################################################################################################################

############################################
# Loading Libraries 

library("tidyverse")
library("numDeriv")
library("rootSolve")
library("geex")

############################################
# EXAMPLE 1: LOGISTIC REGRESSION

############################################
# Data 

dat <- tibble(anemia=c(rep(0,4),rep(1,4)),
              bp=rep(c(rep(0,2),rep(1,2)),2),
              ptb=rep(c(0,1),4),
              n=c(496,74,113,25,85,15,15,3)) %>%
  uncount(n) 

n <- nrow(dat)     # Number of observations

############################################
# Regression by MLE 

reg_mle <- glm(ptb ~ anemia + bp, data=dat, family="binomial")   

ests_mle <-as.data.frame(cbind("beta"=reg_mle$coefficients,         # Formatting results to display in table
                               "se"=sqrt(diag(vcov(reg_mle))))) %>%
  mutate(lcl = beta - 1.96*se,                                      # Construct confidence intervals
         ucl = beta + 1.95*se)
intlist <- row.names(ests_mle)                                      # Label rows of table

print("Estimated logistic regression") 
print("MLE")                                             # Display results
print(round(ests_mle,3))           

                        

############################################
# Defining estimating equation

estimating_function <- function(beta){
  p <- plogis(beta[1] + beta[2]*dat$anemia + beta[3]*dat$bp)  # Predicted probability of the outcome
  ef_1 <- (dat$ptb - p)*1                                     # Estimating function for beta_0 (intercept)
  ef_2 <- (dat$ptb - p)*dat$anemia                            # Estimating function for beta_1
  ef_3 <- (dat$ptb - p)*dat$bp                                # Estimating function for beta_2
  return(cbind(ef_1, ef_2, ef_3)) 
}

estimating_equation <- function(beta){
  estf = estimating_function(beta)        # Return estimating function
  este = colSums(estf)                    # Estimating equations are sum
  return(este)
}

############################################
# Root-finding

proc <- rootSolve::multiroot(f = estimating_equation,     # Function to find root(s) of
                             start = c(-2,0,0))           # Starting values for root-finding procedure

# Starting values need to be provided. A good starting value is within the plausible 
# range and not close to the bounds. For example, if the parameter is a risk then a 
# starting value of 0.5 would be a good choice. For regression, one can generally 
# provide starting values of 0. To increase computational efficiency of M-estimation, 
# subsets of estimating equations can be solved separately and then used as the starting 
# values for the overall estimating equations. For example, in Example 2, we can obtain 
# the point estimates for propensity score model parameters using built-in 
# functions/procedures for logistic regression use those as starting values.

beta_root <- proc$root

############################################
# Baking the bread (approximate derivative)

deriv <- numDeriv::jacobian(func = estimating_equation,   # Function to find derivative of
                            x = beta_root)                # Array of values to compute derivative at (root of estimating equation)

bread <- -1*deriv / n

############################################
# Cooking the filling (matrix algebra)

outerprod <- t(estimating_function(beta_root)) %*% estimating_function(beta_root) # Outer product of the residuals
filling <- outerprod/n 

############################################
# Assembling the sandwich (matrix algebra)

sandwich <- solve(bread) %*% filling %*% t(solve(bread))
se <- sqrt(diag(sandwich / n))                            # Extract diagonal and take sqrt to get se

ests_mest <-as.data.frame(cbind("beta"=beta_root,         # Formatting results to display in table
                                "se"=se)) %>%
  mutate(lcl = beta - 1.96*se,                            # Construct confidence intervals
         ucl = beta + 1.96*se)

row.names(ests_mest) <- intlist                           # Label rows of table
            
print("M-Estimation, by hand")                            # Display results
print(round(ests_mest,3))


####################################################################################################################
# Using geex 

geex_ef <- function(data){               # Function of estimating functions (to be used in geex::m_estimate)
  ptb <- data$ptb
  anemia <- data$anemia
  bp <- data$bp
  function(theta){
    p <- plogis(theta[1] + theta[2]*anemia + theta[3]*bp)   # Predicted probability of the outcome
    c((ptb - p)*1,                                          # Estimating function for beta_0 (intercept)
      (ptb - p)*anemia,                                     # Estimating function for beta_1
      (ptb - p)*bp)                                         # Estimating function for beta_2
  }
}

# Compute M-estimator using geex package  
mestr <- geex::m_estimate(estFUN = geex_ef,                                       # Function of estimating functions
                          data = dat,                                             # Data to be used (must be data frame)
                          root_control = setup_root_control(start = c(-2,0,0)))   # Set starting values

beta_geex <- roots(mestr)             # Extract point estimates from geex
se_geex <- sqrt(diag(vcov(mestr)))    # Extract finite sample variance and take sqrt to get se

ests_geex <-as.data.frame(cbind("beta"=beta_geex,   # Formatting results to display in table
                                "se"=se_geex)) %>%
  mutate(lcl = beta - 1.96*se,                      # Construct confidence intervals
         ucl = beta + 1.96*se)

row.names(ests_geex) <- intlist                     # Label rows of table

print("M-Estimation, by geex")                      # Display results
print(round(ests_geex,3))

############################################
# Results for logistic regression

print("Estimated logistic regression")
print("MLE")
print(round(ests_mle,3))

print("M-Estimation, by hand")
print(round(ests_mest,3))

print("M-Estimation, using geex")
print(round(ests_geex,3))


################################################
# EXAMPLE 2: STANDARDIZATION BY IPW

################################################
# Using geex 

geex_ef2 <- function(data){             # Function of estimating functions (to be used in geex::m_estimate)  
  ptb <- data$ptb
  anemia <- data$anemia
  bp <- data$bp
  
  function(theta){
    alpha <- theta[1:2]                            # Parameters for propensity score model
    mu <- theta[3:4]                               # Risks
    delta <- theta[5]                              # Risk difference
    
    pscore <- plogis(alpha[1] + alpha[2]*bp)       # Predicted propensity score
    ef_1 <- (anemia - pscore)                      # Estimating function for alpha_0
    ef_2 <- (anemia - pscore)*bp                   # Estimating function for alpha_1
    
    wt <- anemia/pscore + (1-anemia)/(1-pscore)    # IPW
    ef_r0 <- (1 - anemia)*wt*ptb - mu[1]           # Estimating function for mu_0
    ef_r1 <- anemia*wt*ptb - mu[2]                 # Estimating function for mu_1

    ef_rd <- mu[2] - mu[1] - delta[1]              # Estimating function for delta
    return(c(ef_1,ef_2,ef_r0,ef_r1,ef_rd))
  }
}

mest_2 <- geex::m_estimate(estFUN = geex_ef2,   # Compute M-estimator using geex package                                      
                            data = dat,                                             
                            root_control = setup_root_control(start = c(-2,0,.1,.1,0)))   

theta2 <- roots(mest_2)                         # Extract point estimates from geex              
se2 <- sqrt(diag(vcov(mest_2)))                 # Extract finite sample variance and take sqrt to get se 

ests_2 <-as.data.frame(cbind("beta"=theta2,     # Formatting results to display in table
                              "se"=se2)) %>%
  mutate(lcl = beta - 1.96*se,                  # Construct confidence intervals
         ucl = beta + 1.96*se)

row.names(ests_2) <- c(intlist[1],intlist[3],"risk0","risk1","rd") # Label rows of table

print("IPW by M-Estimation, using geex")        # Display results
print(round(ests_2,2))


################################################
# EXAMPLE 3: STANDARDIZATION BY G-COMPUTATION

################################################
# Using geex 

geex_ef3 <- function(data){                       # Function of estimating functions (to be used in geex::m_estimate)            
  ptb <- data$ptb
  anemia <- data$anemia
  bp <- data$bp
  
  function(theta){
    beta <- theta[1:3]                                # Parameters of outcome model
    mu <- theta[4:5]                                  # Risks
    delta <- theta[6]                                 # Risk difference
    
    p <- plogis(beta[1] + beta[2]*anemia + beta[3]*bp)  # Predicted probability of the outcome
    
    ef_1 <- (ptb - p)                                 # Estimating function for beta_0  
    ef_2 <- (ptb - p)*anemia                          # Estimating function for beta_1
    ef_3 <- (ptb - p)*bp                              # Estimating function for beta_2 
    
    ef_r0 <- plogis(beta[1] + beta[2]*0 + beta[3]*bp) - mu[1]  # Estimating function for mu_0
    ef_r1 <- plogis(beta[1] + beta[2]*1 + beta[3]*bp) - mu[2]  # Estimating function for mu_1
    
    ef_rd <- mu[2] - mu[1] - delta[1]                 # Estimating function for delta
    return(c(ef_1,ef_2,ef_3,ef_r0,ef_r1,ef_rd))
  }
}

mest_3 <- geex::m_estimate(estFUN = geex_ef3,       # Compute M-estimator using geex package                                  
                            data = dat,                                             
                            root_control = setup_root_control(start = c(-2,0,0,.1,.1,0)))   

theta3 <- roots(mest_3)                             # Extract point estimates from geex                 
se3 <- sqrt(diag(vcov(mest_3)))                     # Extract finite sample variance and take sqrt to get se

ests_3 <-as.data.frame(cbind("beta"=theta3,         # Formatting results to display in table
                              "se"=se3)) %>%
  mutate(lcl = beta - 1.96*se,                      # Construct confidence intervals
         ucl = beta + 1.96*se)

row.names(ests_3) <- c(intlist,"risk0","risk1","rd")  # Label rows of table

print("G-computation by M-Estimation, using geex")   # Display results
print(round(ests_3,3))


############################################
# EXAMPLE 4: OUTCOME MISCLASSIFICATION

############################################
# Data 

datfusion <- tibble(r=c(rep(1,950),rep(0,331)),
              y=c(rep(0,950),                           # Missing Y in R=1 set to 0
                  rep(1,242),rep(0,89)),
              w=c(rep(1,680),rep(0,950-680),
                  rep(1,204),rep(0,38),rep(1,18),rep(0,71)))



################################################
# Using geex 

geex_ef4 <- function(data){                       # Function of estimating functions (to be used in geex::m_estimate)               
  r <- data$r
  y <- data$y
  w <- data$w
  
  function(theta){
    ef_1 <- r*(w - theta[1])                      # Estimating function for misclassified prevalence
    ef_2 <- (1 - r)*y*(w - theta[2])              # Estimating function for sensitivity
    ef_3 <- (1 - r)*(1 - y)*((1 - w) - theta[3])  # Estimating function for specificity
    ef_4 <- theta[4]*(theta[2] - (1 - theta[3])) - (theta[1] - (1 - theta[3]))  # Estimating function for prevalence (Rogan Gladen)

    return(c(ef_1,ef_2,ef_3,ef_4))
  }
}

mest_4 <- geex::m_estimate(estFUN = geex_ef4,        # Compute M-estimator using geex package                               
                           data = datfusion,                                             
                           root_control = setup_root_control(start = c(.7,1,1,0.7)))   

theta4 <- roots(mest_4)                              # Extract point estimates from geex  
se4 <- sqrt(diag(vcov(mest_4)))                      # Extract finite sample variance and take sqrt to get se

ests_4 <-as.data.frame(cbind("beta"=theta4,          # Formatting results to display in table
                              "se"=se4)) %>%
  mutate(lcl = beta - 1.96*se,                       # Construct confidence intervals
         ucl = beta + 1.96*se)

row.names(ests_4) <- c("pr_w","se","sp","pr_y")      # Label rows of table

print("ME correction, using geex")                   # Display results
print(round(ests_4,3))

