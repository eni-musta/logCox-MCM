# compute the Nadaraya-Watson weights for the Beran estimator 
# input:  xx = matrix containing the covariates X as columns (the continuous covariates come first) ordeed with respect to the follow-up times, 
#          h = vector of bandwidths for each continuous covariate in X 
# output:  w = matrix with elements w_ij=K_h(x_j-x_i)/sum_l K_h(x_l-x_i), the Epanechnikov kernel is used

weight = function(xx,h)
{
  n = dim(xx)[1]       # sample size
  m = dim(xx)[2]       # number of covariates 
  d = length(h)        # number of continuous covariates 
  w = matrix(1,n,n)
  
  for(i in 1:d){
     xxi = xx[,i]
     xxi = matrix(xxi,n,n,byrow=TRUE)
      zi = (xxi-t(xxi))/h[i]
      zi = replace(zi,abs(zi)>1,-1)  
      wi = 3/4*(1-zi^2) 
       w = w*wi
  }
  if(m>d){
      for(i in (d+1):m){
        xxi = xx[,i]
        xxi = matrix(xxi,n,n,byrow=TRUE)
        zi = (xxi-t(xxi))
        wi = (zi==0)
        mode(wi) = 'numeric'
        w = w*wi
      }
  }
  
  sumw = apply(w,1,sum)  
  w = w/matrix(sumw,n,n,byrow=FALSE)  
  return(w)  
}


# Compute the nonparametric estimator of the cure probability pi_hat, i.e. the Beran estimator at the last observed event time Y_{(n)}
# input:  data = dataframe containing the follow-up time, censoring status and covariates of the logistic component ordered with respect to the follow-up times
#         w = Nadaraya-Watson weights returned by the previous function

beran = function(data,w)
{
  m = nrow(w)
  n = ncol(w)
  k = cbind(rep(0,m),w[,1:(n-1)])
  cumsumk = t(apply(k,1,cumsum))   
  cumw = 1-cumsumk               # H([t,infinit)|x]), the  j col is t=Y.(j) 
  w = replace(w,cumw==0,0)  
  cumw = replace(cumw,cumw==0,1)  
  q = 1 - w/cumw 
  delta = matrix(data[,2],m,n,byrow=TRUE)
  q = replace(q,delta==0,1)  
  q[which(q<10^(-8))] = 0 # to avoid numerical problems when dealing with extremely small numbers
  q = apply(q,1,prod)  # P(T=infinity|x) = S(Y_(n)|x) 
  return(q)
}


# Logarithm of the logistic likelihood 
# L(gamma)=product_i {phi(gamma,X_i)^(1-pi_hat(X_i))}{[1-phi(gamma,X_i)]^pi_hat(X_i)}
# input: gamma parameter and a matrix V containing the nonparamtric estimator pi_hat as first column and covariates X in the other columns

logLik_logit = function(par,V) 
{ 
  m = length(par)
  phi = exp(-par[1]-c(par[2:m]%*%t(V[,2:m])))   
  phi = 1/(1+phi)
  likelihood = sum((1-V[,1])*log(phi)+V[,1]*log(1-phi))
  return(likelihood)
}


# Estimates baseline survival function
# input: Time = follow-up times, Status = censoring indicator, Z = covariates for the latency
#        beta = estimator of beta coefficients, w = weights equal to P('not cured'| observations, estimators from the previous step)
# zero-tail constraint: S_0(t)=0 for t>Y_(r)

surv_cox = function(Time,Status,Z,beta,w){    
    death_point = sort(unique(subset(Time, Status==1)))
    coxexp = exp((beta)%*%t(Z))  
    lambda = numeric()
    event = numeric()
    for(i in 1: length(death_point)){
      event[i] = sum(Status*as.numeric(Time==death_point[i]))
      temp = sum(as.numeric(Time>=death_point[i])*w*drop(coxexp))
      lambda[i] = event[i]/temp
    }
    CumHazard = numeric()
    for(i in 1:length(Time)){
      CumHazard[i] = sum(as.numeric(Time[i]>=death_point)*lambda)
      if(Time[i]>max(death_point))CumHazard[i] = Inf
      if(Time[i]<min(death_point))CumHazard[i] = 0
    }
    survival = exp(-CumHazard)
    list(survival=survival) 
  }



# EM algorithm for the cox component, does not update the logstic parameters
# input:  Time = follow-up times, Status = censoring indicator, X = covariates for the incidence, Z = covariates for the latency
#         gamma = estimator of the parameters of the incidence from Step 1, beta_init = initial value for beta
#         enmax = maximum number of iterations, eps = tolerance for convergence
# output:  list containing the estimator of beta, of the baseline survival function, convergence (distance between the old and the new estimators in the last iteration)

EM = function(Time,Status,X,Z,gamma,beta_init,emmax,eps){ 
  
    uncureprob = matrix(1/(1+exp(-(gamma)%*%t(X))),ncol=1)  # estimate probability of being noncured  using the estimator of gamma from Step 1
   
    w = Status	 # initial value of the weights
    beta=beta_init
    s <- surv_cox(Time,Status,Z,beta,w)$survival  #initial estimate for the baseline survival function ignoring the cure proportion and using the initial estimate of beta
    
    convergence = 1000
    i = 1
    
    while (convergence > eps & i < emmax){  
 
      ## E step: update the weights 
      
      survival = drop(s^(exp((beta)%*%t(Z))))
      w = Status+(1-Status)*(uncureprob*survival)/((1-uncureprob)+uncureprob*survival)
      
      ## M step: update the estimators
      
      update_beta = coxph(Surv(Time, Status)~Z+offset(log(w)), subset=w!=0, method="breslow")$coef
      update_s = surv_cox(Time,Status,Z,beta,w)$survival
        
      convergence = sum(c(update_beta-beta)^2)+sum((s-update_s)^2)
      beta = update_beta 
      s = update_s
      i = i+1
    }
    
    EM = list(beta=beta,survival=s,tau=convergence)
  }


# Function printing the results (estimator, estimated standard deviation and p-values)
# input: x = list containing the fitted parameters, their estimated standad deviations and p-values

print_logCox_MCM = function (x){
  
  cat("\nCure probability model:\n")
  
    gamma = array(x$gamma, c(length(x$gamma), 4))
    rownames(gamma) = x$gammanm
    colnames(gamma) = c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
    gamma[, 2] = x$gamma_sd
    gamma[, 3] = x$gamma_zvalue
    gamma[, 4] = x$gamma_pvalue
  
  print(gamma)
  cat("\n")
  cat("\nFailure time distribution model:\n")
  
    beta = array(x$beta, c(length(x$beta), 4))
    rownames(beta) = x$betanm
    colnames(beta) = c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
    beta[, 2] <- x$beta_sd
    beta[, 3] <- x$beta_zvalue
    beta[, 4] <- x$beta_pvalue
  
  print(beta)
  invisible(x)
}
