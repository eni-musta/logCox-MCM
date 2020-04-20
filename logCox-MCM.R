###################################################################################
# A presmoothing approach for estimation in the logistic/Cox mixture cure model   #
###################################################################################

# logistic/Cox model
# Incidence:  P('cured'|X)=1-phi(gamma_0,X),   phi(gamma,X)=1/(1+exp(-gamma*X)),
#             where gamma_0 is the vector of unknown parameters and X is the vector of covariates
# Latency:  for the uncured subjects, S_u(t|Z)=S_0(t)^exp(beta*Z),
#           where beta_0 is a vectorof unknown parameters and S_0 is the baseline survival function


# Libraries

library(smcure)
library(np)
library(survival)

# Predefined functions

source("functions_logCox_MCM")

# Load the data into the dataframe "Data" and remove the observations with missing values. 
# The dataframe consists of the following columns: 
# Y - follow-up time
# Delta - censoring status:  censored (Delta=0), noncensored (Delta=1)
# Treatment - a binary covariate (0=control, 1=treatment)
# Age - a continuous covariate centered to the mean
# Gender - a binary covariate (0=male, 1=female)

data("e1684")
Data = e1684   
Data = na.omit(Data) 
colnames(Data) = c("Treatment","Y","Delta","Age","Sex")


# Plot the Kaplan Meier estimator of the survival function

plot(survfit(Surv(Data$Y,Data$Delta)~1))


# Compute the maximum likelihood estimators with the package smcure

MCM = smcure(Surv(Y,Delta)~Age+Treatment+Sex,cureform=~Age+Treatment+Sex,data=Data,model="ph",Var=TRUE,nboot=500)

#smcure estimators
gamma_smcure = MCM$b
beta_smcure = MCM$beta
s_smcure = MCM$s


################
# New approach #
################


# Step 1. Estimate the parameters of the incidence 
    

# data_log is a dataframe containing the follow-up times, censoring status and the covariates used to model the incidence
# In this case all covariates are used for both incidence and latency 
# The continuous covariates are standardized and come before the discrete covariates (as columns of data_log)
# The observations (rows) are ordered according to the follow-up time.

  data_log = data.frame(cbind(Data$Y,Data$Delta,Data$Age/sd(Data$Age),Data$Treatment,Data$Sex))
  colnames(data_log) = c('Y','Delta','Age','Treatment','Sex')
  ord = order(data_log[,1],1-data_log[,2])
  data_log = data_log[ord,]
  

  g_init = glm(Delta~Age+Treatment+Sex,data_log,family=quasibinomial)$coefficients  # initail value of gamma computed by logistic regression   
  
  Y_r = max(data_log[which(data_log[,2]==1),1])    # the largest observed event time
  hopt = npcdistbw(Y~Age,data=data_log,gydat=data_log[which(data_log[,1]<=Y_r),1],cxkertype="epanechnikov")$xbw  # compute bandwidth for estimation of  H(t|X)
  cmax = 2
  hopt = min(hopt,cmax)   # truncate bandwidth from above
  
  
  w = weight(data_log[,3:5],hopt) # compute the Nadaray-Watson weights for the Beran estimator
  pi_hat = beran(data_log,w)      # nonparametric estimator of cure probabilities
  g = optim(par=g_init,logLik_logit,V=cbind(pi_hat,data_log[,3:5]),method='BFGS',control=list(fnscale=-1,maxit=50))$par  #estimator using logistic-likelihood 
  
  gamma_new = g/c(1,sd(Data$Age),1,1)  # rescale the estimator so that it corresponds to the initial covariates (not standardized ones)
  
  
  
# Step 2. Estimate the regression coefficients and the baseline hazard of the Cox component 
  
  
  n = dim(Data)[1]  #sample size
  X = as.matrix(cbind(rep(1, n), Data$Age, Data$Treatment, Data$Sex)) # covariates for the incidence including the intercept 
  Z = as.matrix(cbind(Data$Age, Data$Treatment, Data$Sex)) #covariates for the latency
  Time = Data$Y          # Follow-up time
  Status = Data$Delta    # Censoring status
  
  eps = 1e-07     # tolerance criteria for the convergence of the EM algorithm
  emmax = 50      # maximum number of iterations for the EM algorithm
  
  beta_init = coxph(Surv(Data$Y, Data$Delta) ~ Age+Treatment+Sex,data=Data, method = "breslow")$coef # initial estimator for beta as in the classical Cox model without cure
  
  EM_fit = EM(Time, Status, X, Z, gamma_new, beta_init, emmax, eps)   # perform the EM algorithm
  
  beta_new = EM_fit$beta
  s_new = EM_fit$survival
  
   
# Estimate the variance via bootstrap 
  
  nboot = 500        # number of bootstrap samples
  nbeta = dim(Z)[2]  # number of covariates for the latency
  ngamma = dim(X)[2] # number of covariates for the incidence
  log_var=c('Intercept','Age','Treatment','Sex')
  cox_var=c('Age','Treatment','Sex')
    
  gamma_boot = matrix(rep(0, nboot * ngamma), nrow = nboot)
  beta_boot = matrix(rep(0, nboot * nbeta), nrow = nboot)
  iter = matrix(rep(0, nboot), ncol = 1)
   
  tempdata = cbind(Time, Status, X, Z)
  data1 = subset(tempdata, Status == 1)
  data0 = subset(tempdata, Status == 0)
  n1 = nrow(data1)
  n0 = nrow(data0)
    
  i = 1
  while (i <= nboot) {
    
    # generate the bootstrap sample
    
      id1 = sample(1:n1, n1, replace = TRUE)
      id0 = sample(1:n0, n0, replace = TRUE)
      bootdata = rbind(data1[id1, ], data0[id0, ])
      bootZ = bootdata[,(3+ngamma):(2+ngamma+nbeta)]
      bootX = bootdata[,3:(2+ngamma)]
        
    # estimate the parametersr of the incidence  
      
      data_log.boot = data.frame(cbind(bootdata[,1],bootdata[,2],bootX[,-1]))
      colnames(data_log.boot) = c('Y','Delta',log_var)
      data_log.boot$Age = data_log.boot$Age/sd(data_log.boot$Age)
      ord.boot = order(data_log.boot[,1],1-data_log.boot[,2])
      data_log.boot = data_log.boot[ord.boot,]
      
      
      w.boot = weight(data_log.boot[,3:(1+ngamma)],hopt) 
      pi_hat.boot = beran(data_log.star,w.star)  #nonparametric estimaytor of cure probabilities
      g.boot = optim(par=g,logLik_logit,V=cbind(pi_hat.boot,data_log.boot[,3:(1+ngamma)]),method='BFGS',control=list(fnscale=-1,maxit=50))$par  #estimator using logistic-likelihood, 
      
      gamma.boot = g.boot/c(1,sd(bootX[,2]),rep(1,ngamma-2))
      
    # estimate the Cox component 
      
      bootfit <- EM(bootdata[, 1], bootdata[, 2], bootX, bootZ, g.star, beta, emmax, eps)
      
      gamma_boot[i, ] = gamma.boot
      beta_boot[i, ] = bootfit$beta
      
      if (bootfit$tau < eps)
        i = i + 1
    }
    
    gamma_var = apply(gamma_boot, 2, var)
    beta_var = apply(beta_boot, 2, var)
    gamma_sd = sqrt(gamma_var)
    beta_sd = sqrt(beta_var)
  
  fit = list()
  fit$beta = beta_new
  fit$gamma = gamma_new
  fit$gamma_var = gamma_var
  fit$gamma_sd = gamma_sd
  fit$gamma_zvalue = gamma_new/gamma_sd
  fit$gamma_pvalue = (1 - pnorm(abs(fit$gamma_zvalue)))*2
  fit$beta_var = beta_var
  fit$beta_sd = beta_sd
  fit$beta_zvalue = beta_new/beta_sd
  fit$beta_pvalue = (1 - pnorm(abs(fit$beta_zvalue)))*2
  fit$gammanm =  log_var
  fit$betanm = cox_var
  fit$s = s
  print_logCox_MCM(fit)  


