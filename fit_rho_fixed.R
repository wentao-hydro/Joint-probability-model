fit_rho_fixed <-function(fcst,obs, threshold_fcst,threshold_obs, obs_par,fcst_par)
{
#coded by Wentao Li, wentaoli@student.bnu.edu.cn

  # fit the parameter of correlation coefficient in Normal space
  # MLE: minimize the negative log likelihood function
  
  par0 <- c(0.25)
  opt.res <- optim(par0,
                   likelihood_fixed_corr,
                   fcst = fcst,
                   obs = obs,
                   threshold_fcst = threshold_fcst,
                   threshold_obs = threshold_obs,
                   obs_par = obs_par,
                   fcst_par = fcst_par,
                   method = "L-BFGS-B",
                   lower = 0.001,
                   upper = 0.95 )   


  neg_loglikelihood=opt.res$value
  N_sample=length(fcst)
  N_par=1  #one parameter, rho
  AIC=2*N_par + 2*neg_loglikelihood
  BIC=log(N_sample)*N_par + 2*neg_loglikelihood
  
  
  list(par=opt.res$par, neg_loglikelihood=c(neg_loglikelihood,AIC,BIC) ) 


}
  
  
  likelihood_fixed_corr <- function(par_corr, fcst,obs, threshold_fcst,threshold_obs,obs_par,fcst_par )  
  {
    ind_pos_fcst <- fcst >threshold_fcst
    ind_pos_obs <- obs >threshold_obs
    mu1<-fcst_par[1]
    sigma1<-fcst_par[2]    
    mu2<-obs_par[1]
    sigma2<-obs_par[2]

 
    #four cases of likelihood functions    
    func11 <- dmvnorm(x=cbind(fcst,obs), mean=c(mu1,mu2),sigma=matrix(c(sigma1^2, par_corr*sigma1*sigma2, par_corr*sigma1*sigma2, sigma2^2), 2 ) )
    
    func10 <- dnorm(fcst, mean=mu1, sd=sigma1) * pnorm(threshold_obs,mean=(mu2 + par_corr*sigma2/sigma1*(fcst-mu1)), sd=sqrt(1-par_corr^2)*sigma2 )#fcst>0,obs=0
    
    func01 <- dnorm(obs, mean=mu2, sd=sigma2) * pnorm(threshold_fcst,mean=(mu1 + par_corr*sigma1/sigma2*(obs-mu2)),  sd=sqrt(1-par_corr^2)*sigma1 )#fcst=0,obs>0
    
    func00 <- pmvnorm(lower=-Inf, upper=c(threshold_fcst,threshold_obs), mean=c(mu1,mu2),sigma=matrix(c(sigma1^2, par_corr*sigma1*sigma2, par_corr*sigma1*sigma2, sigma2^2), 2 ) )
    
    

    f1<-  ( func11 * ((ind_pos_fcst==1) & (ind_pos_obs ==1)) )#if not satisfy the condition in each case, the f value equals zero
    f2<-  ( func10 * ((ind_pos_fcst==1) & (ind_pos_obs ==0)) ) 
    f3<-  ( func01 * ((ind_pos_fcst==0) & (ind_pos_obs ==1)) ) 
    f4<-  ( func00 * ((ind_pos_fcst==0) & (ind_pos_obs ==0)) ) 

    
    

    #negative log likelihood function value
    neg_loglikelihood <- -sum(log(f1[f1>0]), na.rm=T) - sum(log(f2[f2>0]), na.rm=T) - sum(log(f3[f3>0]), na.rm=T) - sum(log(f4[f4>0]), na.rm=T)




    return( neg_loglikelihood )#return negative log likelihood function value
    
  }
  
  