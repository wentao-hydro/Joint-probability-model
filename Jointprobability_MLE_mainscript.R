#coded by Wentao Li, wentaoli@student.bnu.edu.cn

rm(list=ls())

library(MASS)
library(neldermead)#optimization method of downhill simplex by Nelder and Mead
library(mvtnorm) #multivariate Normal

#load functions: fit the joint probability model using MLE, then make prediction given new forecasts
BJP_function_path<-'D:\\Auxiliary_BJP_MLE_logsinh_fixrho_commented\\'  #change the R function diretory
source(paste(BJP_function_path,'fitparameter_BJP_logsinh_CMLE_fixrho.R',sep=''))#
source(paste(BJP_function_path,'fit_rho_fixed.R',sep='')) 
source(paste(BJP_function_path,'fit_logsinh_transform.R',sep=''))#
source(paste(BJP_function_path,'fcst2ensem_logsinh_fixrho.R',sep=''))#

#load example data
load( file='D:\\Auxiliary_BJP_MLE_logsinh_fixrho_commented\\fcst_obs_example.Rdata')#change the example data diretory



#fit the joint probability model
threshold_zero=0.1 #threshold for zero precipitation 
result_list<-fitparameter_BJP_logsinh_CMLE_fixrho(obs_train,fcst_oneevent_train_mean, threshold_zero  )
parameter_list=result_list$parameter_list



#predict given a new forecast; still post-processing for the training data here, just as an example
n_mem_output<-1000  #ensemble size, need how many ensemble members; 
#if using random sampling like the current code, the ensemble size should be large enough, like 500 or 1000; 
#if want less ensemble members, I think can compute the quantiles from the 1000 members 

n_fcst<-length(fcst_oneevent_train_mean)
ens_y<- array(dim=c(n_fcst,n_mem_output))
for (ifcst in 1:n_fcst)
{
  x_fcstnew=fcst_oneevent_train_mean[ifcst] #still post-processing for the training data here, just as an example
   ens_y [ifcst, ]<-fcst2ensem_logsinh_fixrho(x_fcstnew,parameter_list,n_mem_output) #predict for the forecast on each day
}

ens_y_mean<-apply(ens_y, 1, mean)
cat('raw forecast mean corr: ',cor(fcst_oneevent_train_mean,obs_train ),'\n')
cat('post-processed forecast mean, corr: ',cor(ens_y_mean,obs_train ), '\n')
 


