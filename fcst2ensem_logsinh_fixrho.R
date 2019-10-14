fcst2ensem_logsinh_fixrho <-function(x_fcstnew,parameter_list,n_mem_output )
{  
  #coded by Wentao Li, wentaoli@student.bnu.edu.cn
  # make prediction given one new forecast 
  # 1.transform
  # 2.compute the conditianal distribution in Normal space: two situations
  # 3.back-transform

  
   #-----------------------------------------------
  # get fitted parameters  
  logsinh_fcst=parameter_list$logsinh_fcst
  logsinh_obs=parameter_list$logsinh_obs 
  threshold_fcst=parameter_list$threshold_fcst 
  threshold_obs=parameter_list$threshold_obs
  zx_fcst_threshold=parameter_list$z_fcst_thre 
  zx_obs_threshold=parameter_list$z_obs_thre
  obs_par=parameter_list$obs_par
  fcst_par=parameter_list$fcst_par
  mu1<-fcst_par[1]
  sigma1<-fcst_par[2]
  mu2<-obs_par[1]
  sigma2<-obs_par[2]
  #-----------------------------------------------

  
  #1. log-sinh transformation for new fcsts:
  zx_fcstnew<-logsinh_function(x_fcstnew,logsinh_fcst ,zx_fcst_threshold)


  
  
  #2. compute conditional distribution in Normal space: two situations 
  parameter_rho=parameter_list$rho0  
  #2.1 F(y|x=0)
  if (x_fcstnew<=threshold_fcst)
  {
    
    threshold_CDF<-pnorm(zx_fcst_threshold, mean=mu1, sd= sigma1)
    rand_CDF<-runif(n_mem_output, 0, threshold_CDF)#generate N CDF value follows uniform distribution [0, CDF(Threshold) ]
    x_samples<-qnorm(rand_CDF, mean=mu1, sd= sigma1)#get N random samples less than the threshold
    
    
    #the conditional distribution of B(v|u<=u0)
    mu.cond <- mu2 + parameter_rho*sigma2/sigma1*(x_samples-mu1) #use the generated "x_samples" here to substitute x in the cond. dist. function
    sigma.cond <- sqrt(1-parameter_rho*parameter_rho) * sigma2
    ens_y <- rnorm(n_mem_output,mu.cond,sigma.cond)#!!   draw random samples here
  }
  #2.2 F(y|x>0)
  else# (x>threshold)
  {
    #conditional distribution of B(v|u>u0)
    mu.cond <- mu2 + parameter_rho*sigma2/sigma1*(zx_fcstnew-mu1)
    sigma.cond <- sqrt(1-parameter_rho*parameter_rho) * sigma2
    ens_y <- rnorm(n_mem_output,mu.cond,sigma.cond)   #!!   draw random samples here
  }
  

  
  
  
  #3. back-transform: use marginal parameters for obs 
  ens_y_ori<-logsinh_inv_function(ens_y, logsinh_obs  ,threshold_obs)

 
  return(ens_y_ori) #return value: y in original space
  
}
# x_samples<-c()
# while (length(x_samples) < n_mem_output)
#   {x_samples_new<-rnorm(n_mem_output, mu1, sigma1)
#   x_samples_used<-x_samples_new[x_samples_new<=zx_fcst_threshold]#only use  samples less than the threshold
#   x_samples<-c(x_samples,x_samples_used )
# }
# x_samples<-x_samples[1:n_mem_output]#get N random samples less than the threshold
