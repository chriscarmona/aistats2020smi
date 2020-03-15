rm(list = ls())
options(scipen=999, stringsAsFactors=FALSE)
set.seed(0)

# pkgbuild::compile_dll()
# Rcpp::compileAttributes()
# roxygen2::roxygenise()
# devtools::document()
# devtools::install()
# devtools::build_vignettes()
# devtools::check()
# devtools::build()

# loading required packages #
req.pck <- c( "aistats2020smi",
              "tidyverse","foreach","doParallel","doRNG","cowplot",
              "abind","mvtnorm" )
req.pck_bool <- sapply(X=req.pck,FUN=require,character.only=T)
if(!all(req.pck_bool)) {
  sapply( X=req.pck[!req.pck_bool], FUN=install.packages, character.only=T);
  sapply( X=req.pck, FUN=require, character.only=T)
}

# Parallel processing
parallel_comp = TRUE
if(parallel_comp){
  n_cores = 6
  options(cores=n_cores)
  doParallel::registerDoParallel()
  getDoParWorkers()
}


# auxilliar functions used in foreach loops #
lrcomb <- function(...) { mapply('rbind', ..., SIMPLIFY=FALSE) }
lacomb <- function(...) { mapply('abind', ..., MoreArgs=list(along=3),SIMPLIFY=FALSE) }
acomb <- function(...) {abind(..., along=3)}

### Generative (true) parameters ###

phi = 0
theta = 1
n = 25 # Sample size for Z
m = 50 # Sample size for Y
sigma_z = 2 # Likelihood variance for Z
sigma_y = 1 # Likelihood variance for Y
sigma_phi = Inf # Prior variance phi
sigma_theta = 0.5 # Prior variance eta

param_names = c('phi','theta','theta_tilde')
param_true = c(phi,theta,theta)

# sequence of eta values in (0,1)
eta_all = seq(0,1,0.025)


# Illustrate convenience of SMI with one synthetic dataset #
if(T) {
  set.seed(123)
  Z = rnorm( n=n, mean=phi, sd=sigma_z)
  Y = rnorm( n=m, mean=phi+theta, sd=sigma_y )
  cat('Z_mean=',mean(Z),'; Y_mean=',mean(Y))

  post_eta_all = foreach::foreach(eta_i = seq_along(eta_all),.combine='lrcomb', .multicombine=TRUE) %dopar% {
    # eta_i = 75
    # Compute posterior mean and variance
    posterior = aistats2020smi::SMI_post_biased_data( Z=Z, Y=Y, sigma_z=sigma_z, sigma_y=sigma_y, sigma_phi=sigma_phi, sigma_theta=sigma_theta, sigma_theta_tilde=sigma_theta, eta=eta_all[eta_i] )
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
  p = cowplot::plot_grid(  post_plot_all[[1]], post_plot_all[[2]], post_plot_all[[3]], ncol=1, align='v' )

  if(F) {
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
    p = cowplot::plot_grid(  mse_plot_all[[1]], mse_plot_all[[2]], mse_plot_all[[3]], ncol=1, align='v' )

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
}

# Illustrate convenience of SMI in expectation #
if(T) {
  set.seed(123)

  # Average Mean Square Error #
  set.seed(123)
  # Compute Posterior mean and sd for each iteration
  n_iter = 1000
  Z = matrix(rnorm( n=n*n_iter, mean=phi, sd=sigma_z),n_iter,n)
  Y = matrix(rnorm( n=m*n_iter, mean=phi+theta, sd=sigma_y ),n_iter,m)
  post_eta_all_iter = foreach(iter_i = 1:n_iter,.combine='lacomb', .multicombine=TRUE)  %:%
    foreach(eta_i = seq_along(eta_all), .combine='lrcomb', .multicombine=TRUE) %dopar% {
      # eta_i=1
      posterior = aistats2020smi::SMI_post_biased_data( Z=Z[iter_i,], Y=Y[iter_i,], sigma_z=sigma_z, sigma_y=sigma_y, sigma_phi=sigma_phi, sigma_theta=sigma_theta, sigma_theta_tilde=sigma_theta, eta=eta_all[eta_i] )
      list( t(posterior[[1]]), diag(posterior[[2]]) )
    }
  # Compute MSE for each iteration
  param_true_array = aperm( array(param_true,dim=c(3,length(eta_all),n_iter)) , c(2,1,3) )
  MSE_all_iter = post_eta_all_iter[[2]] + (post_eta_all_iter[[1]]-param_true_array)^2

  # Average across iterations
  post_eta_all_average = list( apply(post_eta_all_iter[[1]],c(1,2),mean),
                               apply(post_eta_all_iter[[2]],c(1,2),mean) + apply(post_eta_all_iter[[1]],c(1,2),var))
  MSE_average = apply(MSE_all_iter,c(1,2),mean)

  # Plot posterior vs true value
  if(F){
    aistats2020smi::set_ggtheme()
    post_plot_all = foreach::foreach( par_i = seq_along(param_names) ) %do% {
      # par_i=2
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
  }

  # Plot MSE
  aistats2020smi::set_ggtheme()
  mse_plot_all = foreach::foreach( par_i = seq_along(param_names) ) %do% {
    # par_i=1
    p = data.frame( eta=eta_all,
                    mse=MSE_average[,par_i] ) %>%
      ggplot() +
      geom_line( aes(x=eta,y=mse), col='red' ) +
      labs(y=paste("MSE( ",param_names[par_i]," )",sep=""))
    p
  }
  # p = cowplot::plot_grid(  mse_plot_all[[1]], mse_plot_all[[2]], mse_plot_all[[3]], ncol=1, align='v' )
  # print(p)

  # Compare theta vs theta_tilde
  aistats2020smi::set_ggtheme()
  p_mse_theta = MSE_average %>%
    `colnames<-`(param_names) %>%
    as.data.frame() %>%
    mutate(eta=eta_all) %>%
    tidyr::pivot_longer(cols=all_of(param_names),names_to='parameter') %>%
    dplyr::filter(parameter %in% c('theta','theta_tilde')) %>%
    ggplot() +
    geom_line( aes(x=eta,y=value,col=parameter) ) +
    labs(y="MSE average theta")
  # print(p_mse_theta)

  # ELPD approximation via Monte Carlo
  set.seed(123)
  # Compute Posterior mean and sd for each iteration
  n_iter = 1000
  Z = matrix(rnorm( n=n*n_iter, mean=phi, sd=sigma_z),n_iter,n)
  Y = matrix(rnorm( n=m*n_iter, mean=phi+theta, sd=sigma_y ),n_iter,m)

  # generate data from the ground-truth distribution
  n_new = 1000
  Z_new = matrix( rnorm( n=n_iter*n_new, mean=phi, sd=sigma_z), n_iter, n_new )
  Y_new = matrix( rnorm( n=n_iter*n_new, mean=phi+theta, sd=sigma_y), n_iter, n_new )

  log_pred_eta_all_iter = foreach(iter_i = 1:n_iter, .combine='acomb', .multicombine=TRUE) %:%
    foreach( eta_i = seq_along(eta_all), .combine=rbind ) %dopar% {
      # iter_i=1
      # new_i = 1
      # eta_i=1
      # posterior = aistats2020smi::SMI_post_biased_data( Z=Z[iter_i,], Y=Y[iter_i,], sigma_z=sigma_z, sigma_y=sigma_y, sigma_phi=sigma_phi, sigma_theta=sigma_theta, sigma_theta_tilde=sigma_theta, eta=eta_all[eta_i] )
      predictive = aistats2020smi::SMI_pred_biased_data( Z=Z[iter_i,], Y=Y[iter_i,], sigma_z=sigma_z, sigma_y=sigma_y, sigma_phi=sigma_phi, sigma_theta=sigma_theta, sigma_theta_tilde=sigma_theta, eta=eta_all[eta_i] )
      as.numeric( aistats2020smi::dmvnorm_arma( x=cbind(Z_new[iter_i,],Y_new[iter_i,]), mean=as.numeric(predictive[[1]]) , sigma=predictive[[2]], logd=TRUE ) )
    }

  # Average elpd
  elpd_eta_all = apply(log_pred_eta_all_iter,1,mean)
  plot(x=eta_all,y=elpd_eta_all)

  # Optimal eta for each dataset
  elpd_eta_all_datasets = t( apply(log_pred_eta_all_iter,c(1,3),mean) )
  infer_best_data = data.frame( eta_star = eta_all[apply(elpd_eta_all_datasets,1,which.max)] ) %>%
    dplyr::mutate( infer_best = dplyr::case_when( (eta_star==0)~"cut",
                                           (eta_star>0)&(eta_star<1)~"smi",
                                           (eta_star==1)~"full" ) ) %>%
    dplyr::mutate( infer_best = factor(infer_best,levels=c("smi","cut","full")) )

  # Histogram of best eta (minimizing elpd) across datasets
  aistats2020smi::set_ggtheme()
  p_eta_star_hist = infer_best_data %>%
    ggplot(aes(x=eta_star, fill=infer_best)) +
    geom_histogram( bins=30, alpha=0.75 ) +
    theme(legend.position="none") + # Remove legend
    geom_text( aes(x=x,y=y,label=prop), hjust="inward", vjust="inward",
               data = infer_best_data %>%
                 dplyr::group_by(infer_best) %>%
                 dplyr::tally() %>%
                 dplyr::mutate( prop = n/sum(n),
                                x=dplyr::case_when( (infer_best=="cut")~0,
                                                    (infer_best=="smi")~0.3,
                                                    (infer_best=="full")~1 ),
                                y=dplyr::case_when( (infer_best=="cut")~400,
                                                    (infer_best=="smi")~50,
                                                    (infer_best=="full")~100 ) ) )

  # Comparing Distribution of elpd: Cut vs SMI vs Full
  p_elpd_dist = data.frame( cut=elpd_eta_all_datasets[,1],
                            smi=elpd_eta_all_datasets[cbind(1:nrow(elpd_eta_all_datasets),apply(elpd_eta_all_datasets,1,which.max))],
                            full=elpd_eta_all_datasets[,ncol(elpd_eta_all_datasets)] ) %>%
    dplyr::mutate( smi_minus_cut=smi-cut, cut_minus_full=cut-full) %>%
    dplyr::select( c("smi_minus_cut","cut_minus_full") ) %>%
    tidyr::pivot_longer( cols = c("smi_minus_cut","cut_minus_full"), names_to='parameter' ) %>%

    dplyr::mutate( parameter = gsub("_"," ",parameter) ) %>%
    dplyr::mutate( parameter = paste( "elpd", parameter ) ) %>%
    dplyr::mutate( parameter = gsub("minus","- elpd",parameter) ) %>%
    dplyr::rename( elpd=value )

  aistats2020smi::set_ggtheme()
  p = p_elpd_dist %>%
    # filter(elpd!=0) %>%
    ggplot( ) +
    geom_histogram( aes(x=elpd),bins=50,alpha=0.7) +
    facet_wrap( vars(parameter), ncol=1, scales = "free_y" )
  p

  # Comparing Distribution of MSE: Cut vs SMI vs Full
  param_i = match('phi',param_names)
  p_mse_diff_data = data.frame( cut=MSE_all_iter[1,param_i,],
                            smi=MSE_all_iter[ cbind( apply(elpd_eta_all_datasets,1,which.max),
                                                     1,
                                                     1:nrow(elpd_eta_all_datasets) ) ],
                            full=MSE_all_iter[dim(MSE_all_iter)[1],param_i,] ) %>%
    dplyr::mutate( cut_minus_smi=cut-smi, full_minus_cut=full-cut) %>%
    dplyr::select( c("cut_minus_smi","full_minus_cut") ) %>%
    tidyr::pivot_longer( cols = c("cut_minus_smi","full_minus_cut"), names_to='parameter' ) %>%
    dplyr::mutate( parameter = gsub("_minus_"," - ",parameter) ) %>%
    dplyr::mutate( parameter = paste( "MSE(",param_names[param_i],") :", parameter ) )

  aistats2020smi::set_ggtheme()
  p_biased_data_mse_diff = p_mse_diff_data %>%
    # filter(elpd!=0) %>%
    ggplot( ) +
    stat_bin( aes(x=value),bins=25,alpha=0.7, breaks=seq(-1.5,1.5,0.1)) +
    geom_vline(xintercept=0,lty=2)+
    facet_wrap( vars(parameter), ncol=1, scales = "free_y" )
  ggsave( plot=p_biased_data_mse_diff,
          filename="SMI_biased_data_mse_diff.pdf",
          device="pdf", width=15,height=12, units="cm")

  # Plot ELPD and MSE
  aistats2020smi::set_ggtheme()
  elpd_eta_star = data.frame( eta = eta_all[which.max(elpd_eta_all)],
                              elpd = max(elpd_eta_all) )
  curves_data = data.frame( eta=rep(eta_all,3),
                            value=c(-elpd_eta_all,MSE_average[,1],MSE_average[,2]),
                            stat=c(rep("-ELPD( Z, Y )",length(eta_all)),rep("MSE( phi )",length(eta_all)),rep("MSE( theta )",length(eta_all))) )
  labels_data = data.frame( label=c("A","B","C",
                                    "D","E","F",
                                    "G","H","I"),

                            eta=c( eta_all[1], eta_all[length(eta_all)], elpd_eta_star[1,"eta"],
                                   eta_all[1], eta_all[length(eta_all)], eta_all[which.min(MSE_average[,1])],
                                   eta_all[1], eta_all[length(eta_all)], eta_all[which.min(MSE_average[,2])] ),

                            value=c( -elpd_eta_all[1], -elpd_eta_all[length(elpd_eta_all)], -elpd_eta_star[1,"elpd"],
                                     MSE_average[1,1], MSE_average[length(eta_all),1], min(MSE_average[,1]),
                                     MSE_average[1,2], MSE_average[length(eta_all),2], min(MSE_average[,2]) ),

                            stat=c( "-ELPD( Z, Y )","-ELPD( Z, Y )","-ELPD( Z, Y )",
                                    "MSE( phi )","MSE( phi )","MSE( phi )",
                                    "MSE( theta )","MSE( theta )","MSE( theta )" )  )
  p_biased_data = curves_data %>%
    ggplot( aes(x=eta,y=value) ) +
    geom_line( col='red' ) +
    facet_wrap( vars(stat), ncol=1, scales = "free_y" ) +
    geom_vline( aes(xintercept=eta), col="purple", lty=2, data=elpd_eta_star ) +
    geom_point( col="blue", pch=20, size=5,
                data=labels_data ) +
    geom_label( aes(label=label), size=5, hjust="inward", vjust="inward", label.padding=unit(0.1, "lines"),
                data=labels_data ) +
    theme( axis.title.y=element_blank() )
  ggsave( plot=p_biased_data,
          filename="SMI_biased_elpd_MSE.pdf",
          device="pdf", width=15,height=12, units="cm")
}
