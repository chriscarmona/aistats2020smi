#' @keywords internal
lrcomb <- function(...) { mapply('rbind', ..., SIMPLIFY=FALSE) }

#' @importFrom abind abind
#' @keywords internal
lacomb <- function(...) { mapply('abind', ..., MoreArgs=list(along=3),SIMPLIFY=FALSE) }

#' @importFrom abind abind
#' @keywords internal
acomb <- function(...) {abind::abind(..., along=3)}

#' @title SMI on synthetic data, single dataset
#' @description Illustrate convenience of SMI on a randomly generated dataset
#' @param out_dir Directory where outputs are stored
#' @param phi Generative (true) value assumed for phi
#' @param theta Generative (true) value assumed for theta
#' @param n Sample size for Z
#' @param m Sample size for Y
#' @param sigma_z Likelihood variance for Z
#' @param sigma_y Likelihood variance for Y
#' @param sigma_phi Prior std dev phi
#' @param sigma_theta Prior std dev theta
#' @param eta_all Degree of influence of the unbiased module into phi
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom dplyr all_of
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats rnorm
#' @export
smi_sec_5_1_single_dataset <- function(
  out_dir,
  phi = 0,
  theta = 1,
  n = 25,
  m = 50,
  sigma_z = 2,
  sigma_y = 1,
  sigma_phi = Inf,
  sigma_theta = 0.5,
  eta_all = seq(0,1,0.025)
) {

  param_names = c('phi','theta','theta_tilde')
  param_true = c(phi,theta,theta)

  Z = stats::rnorm( n=n, mean=phi, sd=sigma_z)
  Y = stats::rnorm( n=m, mean=phi+theta, sd=sigma_y )
  cat('Z_mean=',mean(Z),'; Y_mean=',mean(Y))

  eta_i = 1
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
  par_i=1
  post_plot_all = foreach::foreach( par_i = seq_along(param_names) ) %do% {

    p = data.frame( eta=eta_all,
                    post_mean=post_eta_all[[1]][,par_i],
                    post_sd=post_eta_all[[2]][,par_i]^0.5 ) %>%
      ggplot() +
      geom_line( aes(x=.data$eta,y=.data$post_mean), col='red' ) +
      geom_line( aes(x=.data$eta,y=.data$post_mean+.data$post_sd), col='blue', lty=3 ) +
      geom_line( aes(x=.data$eta,y=.data$post_mean-.data$post_sd), col='blue', lty=3 ) +
      geom_hline(yintercept=param_true[par_i], lty=2) +
      labs(y=paste(param_names[par_i]," posterior",sep=""))
    p
  }
  p = cowplot::plot_grid(  post_plot_all[[1]], post_plot_all[[2]], post_plot_all[[3]], ncol=1, align='v' )

  ggsave( plot = p,
          filename = paste(out_dir, "/biased_data_posterior_single_dataset.pdf", sep=""),
          device="pdf", width=10,height=15, units="cm")

  if(F) {
    # Plot MSE
    par_i=1
    mse_plot_all = foreach::foreach( par_i = seq_along(param_names) ) %do% {
      # par_i=3
      p = data.frame( eta=eta_all,
                      mse=mse_eta_all[,par_i] ) %>%
        ggplot() +
        geom_line( aes(x=.data$eta,y=.data$mse), col='red' ) +
        labs(y=paste("MSE ",param_names[par_i],sep=""))
      p
    }
    p = cowplot::plot_grid(  mse_plot_all[[1]], mse_plot_all[[2]], mse_plot_all[[3]], ncol=1, align='v' )

    # Compare theta vs theta_tilde
    p = mse_eta_all %>%
      `colnames<-`(param_names) %>%
      as.data.frame() %>%
      mutate(eta=eta_all) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(param_names), names_to='parameter') %>%
      dplyr::filter(.data$parameter %in% c('theta','theta_tilde')) %>%
      ggplot() +
      geom_line( aes(x=.data$eta,y=.data$value,col=.data$parameter) ) +
      labs(y="MSE theta")
    print(p)
  }

  return(TRUE)
}



#' @title SMI on synthetic data, on expectation
#' @description Illustrate convenience of SMI using Monte Carlo integration over the data generation
#' @param out_dir Directory where outputs are stored
#' @param n_iter number of iterations
#' @param n_new number of iterations
#' @param phi Generative (true) value assumed for phi
#' @param theta Generative (true) value assumed for theta
#' @param n Sample size for Z
#' @param m Sample size for Y
#' @param sigma_z Likelihood variance for Z
#' @param sigma_y Likelihood variance for Y
#' @param sigma_phi Prior std dev phi
#' @param sigma_theta Prior std dev eta
#' @param sigma_theta Prior std dev eta
#' @param eta_all Degree of influence of the unbiased module into phi
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom dplyr all_of case_when filter group_by mutate rename select tally
#' @importFrom foreach foreach %do% %dopar% %:%
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats rnorm
#' @importFrom tidyr pivot_longer
#' @export
smi_sec_5_1_expectation <- function(
  out_dir,
  n_iter = 1000,
  n_new = 1000,
  phi = 0,
  theta = 1,
  n = 25,
  m = 50,
  sigma_z = 2,
  sigma_y = 1,
  sigma_phi = Inf,
  sigma_theta = 0.5,
  eta_all = seq(0,1,0.025)
) {

  param_names = c('phi','theta','theta_tilde')
  param_true = c(phi,theta,theta)


  # Average Mean Square Error #

  # Compute Posterior mean and sd for each iteration

  Z = matrix( stats::rnorm( n=n*n_iter, mean=phi, sd=sigma_z), n_iter, n)
  Y = matrix( stats::rnorm( n=m*n_iter, mean=phi+theta, sd=sigma_y ), n_iter, m)

  iter_i = 1; eta_i=1
  post_eta_all_iter = foreach::foreach(iter_i = 1:n_iter,.combine='lacomb', .multicombine=TRUE)  %:%
    foreach::foreach(eta_i = seq_along(eta_all), .combine='lrcomb', .multicombine=TRUE) %dopar% {
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
  aistats2020smi::set_ggtheme()
  par_i=1
  post_plot_all = foreach::foreach( par_i = seq_along(param_names) ) %do% {
    # par_i=2
    p = data.frame( eta=eta_all,
                    post_mean=post_eta_all_average[[1]][,par_i],
                    post_sd=post_eta_all_average[[2]][,par_i]^0.5 ) %>%
      ggplot() +
      geom_line( aes(x=.data$eta,y=.data$post_mean), col='red' ) +
      geom_line( aes(x=.data$eta,y=.data$post_mean+.data$post_sd), col='blue', lty=3 ) +
      geom_line( aes(x=.data$eta,y=.data$post_mean-.data$post_sd), col='blue', lty=3 ) +
      geom_hline(yintercept=param_true[par_i], lty=2) +
      labs(y=paste(param_names[par_i]," posterior",sep=""))
    p
  }
  p = cowplot::plot_grid(  post_plot_all[[1]], post_plot_all[[2]], post_plot_all[[3]], ncol=1, align='v' )
  # print(p)

  # Plot MSE
  aistats2020smi::set_ggtheme()
  par_i=1
  mse_plot_all = foreach::foreach( par_i = seq_along(param_names) ) %do% {
    # par_i=1
    p = data.frame( eta=eta_all,
                    mse=MSE_average[,par_i] ) %>%
      ggplot() +
      geom_line( aes(x=.data$eta,y=.data$mse), col='red' ) +
      labs(y=paste("MSE( ",param_names[par_i]," )",sep=""))
    p
  }
  # p = cowplot::plot_grid(  mse_plot_all[[1]], mse_plot_all[[2]], mse_plot_all[[3]], ncol=1, align='v' )
  p = cowplot::plot_grid(  mse_plot_all[[1]], mse_plot_all[[2]], ncol=1, align='v' )
  ggsave( plot=p,
          filename="biased_data_mse_average.pdf",
          device="pdf", width=10,height=10, units="cm")

  # Compare theta vs theta_tilde
  aistats2020smi::set_ggtheme()
  mse_average_theta = MSE_average %>%
    `colnames<-`(param_names) %>%
    as.data.frame() %>%
    mutate(eta=eta_all) %>%

    tidyr::pivot_longer(cols = dplyr::all_of(param_names),names_to='parameter') %>%
    dplyr::filter(.data$parameter %in% c('theta','theta_tilde')) %>%
    ggplot() +
    geom_line( aes(x=.data$eta, y=.data$value, col=.data$parameter) ) +
    # coord_cartesian(ylim=c(0.28,0.50)) +
    labs(y="MSE( theta )")
  ggsave( plot=mse_average_theta,
          filename="biased_data_mse_average_theta.pdf",
          device="pdf", width=10,height=10, units="cm")


  # ELPD approximation via Monte Carlo
  # generate data from the ground-truth distribution
  Z_new = matrix( stats::rnorm( n=n_iter*n_new, mean=phi, sd=sigma_z), n_iter, n_new )
  Y_new = matrix( stats::rnorm( n=n_iter*n_new, mean=phi+theta, sd=sigma_y), n_iter, n_new )

  iter_i = 1; eta_i=1
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
  # plot(x=eta_all,y=elpd_eta_all)

  p = data.frame(eta=eta_all, elpd=elpd_eta_all) %>%
    ggplot() +
    geom_line( aes(x=.data$eta,y=-.data$elpd),col='red' ) +
    labs(y="- elpd(z,y)")
  ggsave( plot=p,
          filename=paste(out_dir, "/biased_data_elpd.pdf", sep=""),
          device="pdf", width=10,height=7, units="cm")

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
  biased_data_curves = curves_data %>%
    ggplot( aes(x=.data$eta, y=.data$value) ) +
    geom_line( col='red' ) +
    facet_wrap( vars(.data$stat), ncol=1, scales = "free_y" ) +
    # geom_vline( aes(xintercept=eta), col="purple", lty=2, data=elpd_eta_star ) +
    # geom_point( col="blue", pch=20, size=5,
    #             data=labels_data ) +
    geom_label( aes(label=.data$label), size=3, hjust="inward", vjust="inward", label.padding=unit(0.1, "lines"),
                data=labels_data ) +
    theme( axis.title.y=element_blank() )
  ggsave( plot=biased_data_curves,
          filename=paste(out_dir, "/biased_elpd_MSE.pdf", sep=""),
          device="pdf", width=15,height=12, units="cm")

  # Optimal eta for each dataset
  elpd_eta_all_datasets = t( apply(log_pred_eta_all_iter,c(1,3),mean) )
  infer_best_data = data.frame( eta_star = eta_all[apply(elpd_eta_all_datasets,1,which.max)] ) %>%
    dplyr::mutate( infer_best = dplyr::case_when( (.data$eta_star==0)~"cut",
                                           (.data$eta_star>0)&(.data$eta_star<1)~"smi",
                                           (.data$eta_star==1)~"full" ) ) %>%
    dplyr::mutate( infer_best = factor(.data$infer_best,levels=c("smi","cut","full")) )

  # Histogram of best eta (minimizing elpd) across datasets
  aistats2020smi::set_ggtheme()
  biased_data_eta_star_hist = infer_best_data %>%
    ggplot(aes(x=.data$eta_star, fill=.data$infer_best)) +
    geom_histogram( bins=30, alpha=0.75 ) +
    theme(legend.position="none") + # Remove legend
    geom_label( aes( x=.data$x, y=.data$y, label=.data$prop, fill=.data$infer_best),alpha=0.3, hjust="inward", vjust="inward", label.padding=unit(0.1, "lines"),
                data = infer_best_data %>%
                  dplyr::group_by(.data$infer_best) %>%
                  dplyr::tally() %>%
                  dplyr::mutate( prop = paste(round(100*n/sum(n)),"%"),
                                 x=dplyr::case_when( (.data$infer_best=="cut")~0,
                                                     (.data$infer_best=="smi")~0.3,
                                                     (.data$infer_best=="full")~1 ),
                                 y=dplyr::case_when( (.data$infer_best=="cut")~400,
                                                     (.data$infer_best=="smi")~50,
                                                     (.data$infer_best=="full")~100 ) ) )
  ggsave( plot=biased_data_eta_star_hist,
          filename=paste(out_dir, "/biased_data_eta_star_hist.pdf", sep=""),
          device="pdf", width=15,height=12, units="cm")

  # Comparing Distribution of elpd: Cut vs SMI vs Full
  p_elpd_dist = data.frame( cut=elpd_eta_all_datasets[,1],
                            smi=elpd_eta_all_datasets[cbind(1:nrow(elpd_eta_all_datasets),apply(elpd_eta_all_datasets,1,which.max))],
                            full=elpd_eta_all_datasets[,ncol(elpd_eta_all_datasets)] ) %>%
    dplyr::mutate( smi_minus_cut=.data$smi-.data$cut, cut_minus_full=.data$cut-.data$full) %>%
    dplyr::select( c("smi_minus_cut","cut_minus_full") ) %>%
    tidyr::pivot_longer( cols = c("smi_minus_cut","cut_minus_full"), names_to='parameter' ) %>%
    dplyr::mutate( parameter = gsub("_"," ",.data$parameter) ) %>%
    dplyr::mutate( parameter = paste( "elpd", .data$parameter ) ) %>%
    dplyr::mutate( parameter = gsub("minus","- elpd", .data$parameter) ) %>%
    dplyr::rename( elpd=.data$value )

  aistats2020smi::set_ggtheme()
  biased_data_elpd_diff = p_elpd_dist %>%
    # filter(elpd!=0) %>%
    ggplot( ) +
    geom_histogram( aes(x=.data$elpd),bins=50,alpha=0.7) +
    facet_wrap( vars(.data$parameter), ncol=1, scales = "free_y" )
  ggsave( plot=biased_data_elpd_diff,
          filename="biased_data_elpd_diff.pdf",
          device="pdf", width=15,height=12, units="cm")

  # Comparing Distribution of MSE: Cut vs SMI vs Full
  param_i = match('phi',param_names)
  mse_diff_data = data.frame( cut=MSE_all_iter[1,param_i,],
                            smi=MSE_all_iter[ cbind( apply(elpd_eta_all_datasets,1,which.max),
                                                     1,
                                                     1:nrow(elpd_eta_all_datasets) ) ],
                            full=MSE_all_iter[dim(MSE_all_iter)[1],param_i,] ) %>%
    dplyr::mutate( cut_minus_smi=.data$cut-.data$smi, full_minus_cut=.data$full-.data$cut) %>%
    dplyr::select( c("cut_minus_smi","full_minus_cut") ) %>%
    tidyr::pivot_longer( cols = c("cut_minus_smi","full_minus_cut"), names_to='parameter' ) %>%
    dplyr::mutate( parameter = gsub("_minus_"," - ",.data$parameter) ) %>%
    dplyr::mutate( parameter = paste( "MSE(",param_names[param_i],") :", .data$parameter ) )

  aistats2020smi::set_ggtheme()
  biased_data_mse_diff = mse_diff_data %>%
    # filter(elpd!=0) %>%
    ggplot( ) +
    stat_bin( aes(x=.data$value),bins=25,alpha=0.7, breaks=seq(-1.5,1.5,0.1)) +
    geom_vline(xintercept=0, lty=2)+
    facet_wrap( vars(.data$parameter), ncol=1, scales = "free_y" )
  ggsave( plot=biased_data_mse_diff,
          filename=paste(out_dir, "/biased_data_mse_diff.pdf", sep=""),
          device="pdf", width=15,height=12, units="cm")

  biased_data_mse_comparison = mse_diff_data %>%
    dplyr::group_by(.data$parameter) %>%
    dplyr::summarise(mean(.data$value<0),mean(.data$value==0),mean(.data$value>0)) %>%
    as.data.frame()

  ### Combining results for biased data ###
  p1 = cowplot::plot_grid( biased_data_eta_star_hist, biased_data_mse_diff, rel_heights =c(2,3), ncol=1 )
  p_biased_data = cowplot::plot_grid( biased_data_curves, p1, ncol=2 )
  ggsave( plot=p_biased_data,
          filename="SMI_biased_data.pdf",
          device="pdf", width=12,height=10, units="cm")

  return(biased_data_mse_comparison)

}
