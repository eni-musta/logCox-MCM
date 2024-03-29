# logCox-MCM


The logistic/Cox mixture cure model is commonly used to handle survival data in the presence of a cure fraction. It assumes that the population consists of two subpopulations: the cured and the noncured ones. Logistic regression is used to model the incidence and a 
Cox proportional hazards model is used for the latency. Maximum likelihood estimators are computed via the Expectation-Maximization algorithm and are  implemented in the package smcure (Cai et al. (2012)). 

We provide a new two-step estimation method for the logistic/Cox mixture cure model. First, a prelimminary smooth estimator of the 
cure probability is computed independently of the model assumed for the latency. Then, the parameters of the incidence are estimated
by 'projecting' the previous nonparametric estimator on the class of logistic functions. In the second step, we estimate the survival distribution of the susceptible subjects by maximizing the Cox component of the likelihood. In this step, the expectation maximization algorithm is used to compute the estimators of the baseline cumulative hazard and the regression parameters. A detailed description of the method and its theoretical properties can be found in Musta, Patilea and Van Keilegom (2020). 

We demonstrate the method for the melanoma data (ECOG phase III clinical trial e1684) from the smcure package (Cai et al. (2012)). 
The purpose of this study was to evaluate the e_ect of treatment (high dose interferon alpha-2b regimen) as the postoperative adjuvant therapy. The event time is the time from initial treatment to recurrence of melanoma and three covariates have been considered: age (continuous variable centered to the mean), gender (0=male and 1=female) and treatment (0=control and 1=treatment). The data consists 
of 284 observations (after deleting missing data) out of which 196 had recurrence of the melanoma cancer (around 30% censoring). The Kaplan Meier estimator of the survival function can be found in 'KM_estimator.png'.

One can run the example via 'logCox-MCM.R'. This produces estimators for the parameters of the incidence, the regression coefficients and the baseline survival of the Cox component (latency). The bootstrap procedure is used to estimate the variance of the parameter estimates and the resulting p-values. The R file uses the packages 'smcure' for the data and 'np' for bandwidth selection. Additional functions are defined in the file 'functions_logCox-MCM.R'.


References

Cai, C., Y. Zou, Y. Peng, and J. Zhang (2012). smcure: An r-package for estimating semiparametric
mixture cure models. Comput. Meth. Prog. Bio. 108 (3), 1255-1260.

Musta, E., V. Patilea and I. Van Keilegom (2020). A presmoothing approach for estimation in mixture cure models. 
