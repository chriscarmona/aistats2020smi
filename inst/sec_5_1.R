rm(list = ls())
options(scipen=999, stringsAsFactors=FALSE)
set.seed(0)

# pkgbuild::compile_dll()
# Rcpp::compileAttributes()
# devtools::document()
# devtools::install()

# loading required packages #
req.pck <- c( "aistats2020smi",
              "tidyverse","foreach","doParallel","doRNG","cowplot",
              "abind" )
req.pck_bool <- sapply(X=req.pck,FUN=require,character.only=T)
if(!all(req.pck_bool)) {
  sapply(X=req.pck[!req.pck_bool],FUN=install.packages,character.only=T);
  sapply(X=req.pck,FUN=require,character.only=T)
}

# Parallel processing
parallel_comp = TRUE
if(parallel_comp){
  n_cores = 25
  options(cores=n_cores)
  doParallel::registerDoParallel()
  getDoParWorkers()
}


# auxilliar functions used in foreach loops #
comb <- function(...) { mapply('rbind', ..., SIMPLIFY=FALSE) }
acomb <- function(...) { mapply('abind', ..., MoreArgs=list(along=3),SIMPLIFY=FALSE) }

### Generative (true) parameters ###

phi = 0
theta = 1
n=25 # Sample size for Z
m=50 # Sample size for Y
sigma_z = 2 # Likelihood variance for Z
sigma_y = 1 # Likelihood variance for Y
sigma_phi=Inf # Prior variance phi
sigma_theta=0.5 # Prior variance eta

param_names = c('phi','theta','theta_tilde')
param_true = c(phi,theta,theta)

# sequence of eta values in (0,1)
eta_all = seq(0,1,0.01)


# Illustrate convenience of SMI with one synthetic dataset #
if(F) {
  set.seed(123)
  Z = rnorm( n=n, mean=phi, sd=sigma_z)
  Y = rnorm( n=m, mean=phi+theta, sd=sigma_y )
  cat('Z_mean=',mean(Z),'; Y_mean=',mean(Y))
  
  post_eta_all = foreach::foreach(eta = eta_all,.combine='comb', .multicombine=TRUE) %do% {
    # eta = 0.01
    # Compute posterior mean and variance
    posterior = aistats2020smi::SMI_post_biased_data( Z=Z, Y=Y, sigma_z=sigma_z, sigma_y=sigma_y, sigma_phi=sigma_phi, sigma_theta=sigma_theta, sigma_theta_tilde=sigma_theta, eta=eta )
    list( t(posterior[[1]]), diag(posterior[[2]]) )
  }
  # Compute MSE
  mse_eta_all = post_eta_all[[2]] + ( post_eta_all[[1]] - t( matrix(param_true,3,length(eta_all)) ) )^2
  
  # Plot posterior vs true value
  aistats2020smi::set_ggtheme()
  post_plot_all = foreach::foreach( par_i = seq_along(param_names) ) %do% {
    
    p = data.frame( eta=eta_all,
                    post_mean=post_eta_all[[1]][,par_i],
                    post_sd=post_eta_all[[2]][,par_i]^0.5 ) %>%
      ggplot() + 
      geom_line( aes(x=eta,y=post_mean), col='red' ) +
      geom_line( aes(x=eta,y=post_mean+post_sd), col='blue', lty=3 ) +
      geom_line( aes(x=eta,y=post_mean-post_sd), col='blue', lty=3 ) +
      geom_hline(yintercept=param_true[par_i], lty=2) +
      labs(y=paste(param_names[par_i]," posterior",sep=""))
    p
  }
  cowplot::plot_grid(  post_plot_all[[1]], post_plot_all[[2]], post_plot_all[[3]], ncol=1, align='v' )
  
  # Plot MSE
  mse_plot_all = foreach::foreach( par_i = seq_along(param_names) ) %do% {
    # par_i=3
    p = data.frame( eta=eta_all,
                    mse=mse_eta_all[,par_i] ) %>%
      ggplot() + 
      geom_line( aes(x=eta,y=mse), col='red' ) +
      labs(y=paste("MSE ",param_names[par_i],sep=""))
    p
  }
  cowplot::plot_grid(  mse_plot_all[[1]], mse_plot_all[[2]], mse_plot_all[[3]], ncol=1, align='v' )
  
  # Compare theta vs theta_tilde
  p = mse_eta_all %>%
    `colnames<-`(param_names) %>%
    as.data.frame() %>%
    mutate(eta=eta_all) %>%
    tidyr::pivot_longer(cols=all_of(param_names),names_to='parameter') %>%
    dplyr::filter(parameter %in% c('theta','theta_tilde')) %>%
    ggplot() +
    geom_line( aes(x=eta,y=value,col=parameter) ) +
    labs(y="MSE theta")
  print(p)
  
}

# Illustrate convenience of SMI in expectation #
if(T) {
  set.seed(0)
  
  # Average Mean Square Error #
  
  # Compute Posterior mean and sd for each iteration
  n_iter = 1000
  post_eta_all_iter = foreach::foreach(iter_i = 1:n_iter,.combine='acomb', .multicombine=TRUE)  %dorng% {
    # iter_i=1
    Z = rnorm( n=n, mean=phi, sd=sigma_z)
    Y = rnorm( n=m, mean=phi+theta, sd=sigma_y )
    post_eta_all = foreach::foreach(eta_i = seq_along(eta_all), .combine='comb', .multicombine=TRUE) %do% {
      # eta_i=1
      posterior = aistats2020smi::SMI_post_biased_data( Z=Z, Y=Y, sigma_z=sigma_z, sigma_y=sigma_y, sigma_phi=sigma_phi, sigma_theta=sigma_theta, sigma_theta_tilde=sigma_theta, eta=eta_all[eta_i] )
      list( t(posterior[[1]]), diag(posterior[[2]]) )
    }
    # mse_eta_all = post_eta_all[[2]] + ( post_eta_all[[1]] - t( matrix(param_true,3,length(eta_all)) ) )^2
    
    post_eta_all
  }
  
  # Compute MSE for each iteration
  param_true_array = aperm( array(param_true,dim=c(3,length(eta_all),n_iter)) , c(2,1,3) )
  MSE_all_iter = post_eta_all_iter[[2]] + (post_eta_all_iter[[1]]-param_true_array)^2
  
  # # Verifying array computation of MSE
  # eta_i=2; iter_i=230
  # MSE_all_iter[eta_i,,iter_i]
  # post_eta_all_iter[[2]][eta_i,,iter_i] + (post_eta_all_iter[[1]][eta_i,,iter_i]-param_true)^2
  
  # Average across iterations
  post_eta_all_average = list( apply(post_eta_all_iter[[1]],c(1,2),mean),
                               apply(post_eta_all_iter[[2]],c(1,2),mean) + apply(post_eta_all_iter[[1]],c(1,2),var))
  MSE_average = apply(MSE_all_iter,c(1,2),mean)
  
  
  # Plot posterior vs true value
  aistats2020smi::set_ggtheme()
  post_plot_all = foreach::foreach( par_i = seq_along(param_names) ) %do% {
    
    p = data.frame( eta=eta_all,
                    post_mean=post_eta_all_average[[1]][,par_i],
                    post_sd=post_eta_all_average[[2]][,par_i]^0.5 ) %>%
      ggplot() + 
      geom_line( aes(x=eta,y=post_mean), col='red' ) +
      geom_line( aes(x=eta,y=post_mean+post_sd), col='blue', lty=3 ) +
      geom_line( aes(x=eta,y=post_mean-post_sd), col='blue', lty=3 ) +
      geom_hline(yintercept=param_true[par_i], lty=2) +
      labs(y=paste(param_names[par_i]," posterior",sep=""))
    p
  }
  p = cowplot::plot_grid(  post_plot_all[[1]], post_plot_all[[2]], post_plot_all[[3]], ncol=1, align='v' )
  # print(p)
  
  # Plot MSE
  aistats2020smi::set_ggtheme()
  mse_plot_all = foreach::foreach( par_i = seq_along(param_names) ) %do% {
    # par_i=3
    p = data.frame( eta=eta_all,
                    mse=MSE_average[,par_i] ) %>%
      ggplot() + 
      geom_line( aes(x=eta,y=mse), col='red' ) +
      labs(y=paste(param_names[par_i]," MSE",sep=""))
    p
  }
  p = cowplot::plot_grid(  mse_plot_all[[1]], mse_plot_all[[2]], mse_plot_all[[3]], ncol=1, align='v' )
  # print(p)
  
  # Compare theta vs theta_tilde
  aistats2020smi::set_ggtheme()
  p = MSE_average %>%
    `colnames<-`(param_names) %>%
    as.data.frame() %>%
    mutate(eta=eta_all) %>%
    tidyr::pivot_longer(cols=all_of(param_names),names_to='parameter') %>%
    dplyr::filter(parameter %in% c('theta','theta_tilde')) %>%
    ggplot() +
    geom_line( aes(x=eta,y=value,col=parameter) ) +
    labs(y="MSE average theta")
  # print(p)
}
