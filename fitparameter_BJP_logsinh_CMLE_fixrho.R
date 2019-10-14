
fitparameter_BJP_logsinh_CMLE_fixrho<-function(obs_train,fcst_train_mean,thre_fix )
{#coded by Wentao Li, wentaoli@student.bnu.edu.cn
  # parameter fitting by MLE, using log-sinh transformation
  # input: obs, fcst mean, censoring threshold (here I just set a same censoring threshold for obs and fcst)
  # output: fitted parameters, see comments in Line 38-43 for detail
  
  threshold_fcst<-  threshold_obs <- thre_fix #here set same censoring threshold for obs and fcst
  
  
  #1. log-sinh transformation for obs and fcst
  res.fcst <-fit_logsinh_transform(fcst_train_mean, threshold_fcst  )
  res.obs <-fit_logsinh_transform(obs_train, threshold_obs  )
  
  
  z_fcst <-res.fcst$data_logsinhed     # transformed fcst
  z_obs <- res.obs$data_logsinhed
  z_fcst_thre<-res.fcst$zero_threshold_tranformed  #transformed zero threshold for fcst
  z_obs_thre<-res.obs$zero_threshold_tranformed
  logsinh_fcst<- res.fcst$par_logsinh  #  log-sinh parameters for fcst: epsilon, lambda, c
  logsinh_obs<- res.obs$par_logsinh  
  fcst_par<- res.fcst$par_dist         #  Normal distribution parameters for transformed fcst: mean, standard deviation
  obs_par<-  res.obs$par_dist
  
  
  
  
  #2. fit the parameter of correlation coefficient in Normal space
  result_rho <- fit_rho_fixed(z_fcst,z_obs,z_fcst_thre,z_obs_thre,obs_par,fcst_par)
  rho0=result_rho$par 
  neg_loglikelihood_fixedrho=result_rho$neg_loglikelihood 
  
  
  
  
  
  #output:
  parameter_list<-list(threshold_fcst=threshold_fcst,  threshold_obs=threshold_obs , #threshold for fcst and obs in untransformed space
                       logsinh_fcst=logsinh_fcst,logsinh_obs=logsinh_obs,            #log-sinh parameters for fcst and obs 
                       z_obs_thre=z_obs_thre,z_fcst_thre=z_fcst_thre,               #threshold for fcst and obs in transformed space
                       obs_par=obs_par,fcst_par=fcst_par,                          #parameters (mean, standard deviation) for fcst and obs in transformed Normal space
                       rho0=rho0  )                                               # correlation coefficient in Normal space
  
  
  
  
  list(parameter_list=parameter_list,neg_loglikelihood_rho=neg_loglikelihood_fixedrho)
  
  
  
  
}



