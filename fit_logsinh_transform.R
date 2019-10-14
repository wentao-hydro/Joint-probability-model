fit_logsinh_transform <- function(data, zero_threshold, method='MLE' )
{
#coded by Wentao Li, wentaoli@student.bnu.edu.cn
  c <- 5/max(data)
  data_scaled <- c * data
  threshold_scaled<- c * zero_threshold

  #1.Log-sinh par fitting
  if (method=='MLE')
  {logsinh.para <- logsinh_MLE (data_scaled ,threshold_censored=threshold_scaled)#!!use scaled threshold, for scaled data
  }else
    {logsinh.para <- logsinh_MAP (data_scaled , threshold_scaled) 
  }
  #cat( logsinh.para$par[2],logsinh.para$par[1],logsinh.para$par[3],logsinh.para$par[4],'\n') #log(epsilon),log(lambda),mean/sd,log(sd)
  
  
  
  #2.apply log-sinh to data
  epsilon <- exp(logsinh.para$par[1])
  lambda <- exp(logsinh.para$par[2])
  sd <- exp(logsinh.para$par[4])
  mean <- sd * (logsinh.para$par[3])
  par_logsinh=c(epsilon,lambda,c)#output: par in untransformed space
  
  
  zero_threshold_tranformed<-(1/lambda) * log( sinh(epsilon + lambda * (threshold_scaled) ) ) # eq.1
  data_logsinhed <- logsinh_function(data , par_logsinh, zero_threshold_tranformed) # eq.1
 
  
  
  
  #output:
  list( data_logsinhed=data_logsinhed,zero_threshold_tranformed=zero_threshold_tranformed,par_logsinh=par_logsinh,par_dist=c(mean,sd) )
}




#logsinh fitting using MLE 
logsinh_MLE <- function(x , threshold_censored )
{
  #Log-sinh transformation
  #input: x samples
  #output: fitted parameters: epsilon,lambda,mean, standard deviation
  
  #MLE
  par0 <- c( -2, 0.1, -1, 1)  #4 parameter, re-parametrised, (transformed parameters)
  # par0 <- c( -2, 0.1, -1, 1)  #4 parameter, re-parametrised, (transformed parameters)
  
  #-------- global var-----------------
  obs <<- x 
  threshold_censored_ya <<- threshold_censored  
  #------------------------------------
  
  
  #optim: use simplex
  opt.res <- fminsearch(MLE.bjp,
                        x0=par0)
  res<-list()
  res$par<-opt.res$optbase$xopt
  res 
  
  
  
}




#The likelihood function for logsinh fitting using MLE
MLE.bjp<-function(par )
{
  
  #get original parameters:
  epsilon <- exp(par[1])
  lambda <- exp(par[2])
  sd <- exp(par[4]) #standard deviation
  mean <- par[3] * sd
  
  miss <- is.na(obs)
  ya <- obs[!miss]  #data in original space
  yb <- (1/lambda) * log( sinh( epsilon + lambda * ya ) ) # eq.1, data in Normal space
  # cat(summary(yb))
  # browser()
  
  #1. Likelihood function for non-censored data: use density * Jacob
  Jacob_value <- 1 / (tanh(epsilon + lambda * ya) )   #eq.8
  density_norm <- dnorm( yb, mean = mean, sd = sd)
  density_yb <- density_norm * Jacob_value
  
  
  #2 CDF value for censored data (zero precipitation): need to use CDF instead of density
  threshold_censored_yb <- (1/lambda) * log( sinh ( epsilon + lambda * threshold_censored_ya ) )
  CDF_censored <- pnorm( threshold_censored_yb, mean = mean, sd = sd)
  
  
  #3. negative log likelihood for all samples together
  ign_CDF<- -log(CDF_censored ) *  (ya <= threshold_censored_ya)
  ign_CDF<- ign_CDF[(!is.infinite(ign_CDF))  ]
  ign_CDF_sum<-sum(ign_CDF, na.rm=T)
  
  ign_PDF<- -log(density_yb) * (ya > threshold_censored_ya)
  ign_PDF<- ign_PDF[(!is.infinite(ign_PDF))  ]
  ign_PDF_sum<-sum(ign_PDF, na.rm=T)
  neg_loglikelihood <- ign_PDF_sum +  ign_CDF_sum
 

  # 
  # cat(epsilon,lambda,  neg_loglikelihood,'\n')
  # cat(par, neg_loglikelihood,'\n')
  return( neg_loglikelihood )#return negative log likelihood 
  
}





logsinh_function <- function( x ,par_logsinh , zero_threshold_tranformed)
  # input: threshold in transformed space!!!
{ 
  epsilon<-par_logsinh[1]
  lambda<-par_logsinh[2]
  c<-par_logsinh[3]
  
  data_scaled <- c * x #scaled
  data_logsinhed<-(1/lambda) * log( sinh(epsilon + lambda * data_scaled ) ) # eq.1
  data_logsinhed[data_logsinhed<=zero_threshold_tranformed ]<- zero_threshold_tranformed
  return(data_logsinhed)
}




logsinh_inv_function <- function( x, par_logsinh, zero_threshold )
  # input: threshold in original space!!!
{
  epsilon<-par_logsinh[1]
  lambda<-par_logsinh[2]
  c<-par_logsinh[3]
  
  data_inv_transformed <-(1/lambda) * ( asinh(exp(lambda * x)) - epsilon ) # eq.2
  data_inv_transformed <- data_inv_transformed * (1/c)   #scaled back
  # data_inv_transformed[data_inv_transformed<=zero_threshold ]<-zero_threshold #set as the censoring threshold, not zero
  data_inv_transformed[data_inv_transformed<=zero_threshold ]<-0 #set as exactly zero
  return(data_inv_transformed )
  
}










