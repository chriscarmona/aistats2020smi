
#' @title Compute MCMC in HPV model
#' @description Compute MCMC used in section 5.3
#' @param out_dir Directory where outputs are stored
#' @param eta_all Levels of influence of the Poisson module
#' @param force_compute_mcmc To avoid computation time, completed MCMC is not recomputed, unless force_compute_mcmc is True
#' @param n_iter Number of iterations in the main chain
#' @param n_iter_sub Number of iterations in the sub-chain for the two-stage MCMC
#' @param n_warmup Number of iterations used for warm-up
#' @param n_cores Number of cores used by stan
#' @param seed Set random seed before MCMC starts for reproducibility
#' @param pkg_dir Internal directory where the package is installed. Used to retrieve stan scripts for the HPV model.
#' @importFrom foreach foreach %dopar%
#' @importFrom rstan stan_model
#' @export
smi_sec_5_4_compute_mcmc <- function(
  out_dir,
  eta_all = seq(0.10,1.00,by=0.10),
  force_compute_mcmc = FALSE,
  n_iter = 5000,
  n_iter_sub = 500,
  n_warmup = 1000,
  n_cores = 8,
  seed = 123,
  pkg_dir = system.file(package = "aistats2020smi")
) {

  cat("Compiling HPV model in stan...")
  hpv_model_whole_stan = rstan::stan_model(file = paste(pkg_dir, "/stan/hpv_model_whole.stan", sep="") )
  hpv_model_poisson_stan = rstan::stan_model(file = paste(pkg_dir, "/stan/hpv_model_poisson.stan", sep="") )

  cat("MCMC sampling for HPV model starts...")
  i_eta=1
  foreach::foreach( i_eta = seq_along(eta_all) ) %do% {
    # i_eta <- 2

    eta_pois <- eta_all[i_eta]
    cat(eta_pois,", ")
    file_i <- paste(out_dir,"/hpv_smi_",formatC(eta_pois,digits=3,format="f",flag="0"),".rda",sep="")

    if( !file.exists(file_i) | force_compute_mcmc ) {
      set.seed(seed)
      hpv_smi_mcmc_i <- aistats2020smi::mcmc_hpv( HPV=aistats2020smi::HPV,

                                                  hpv_model_whole_stan = hpv_model_whole_stan,
                                                  hpv_model_poisson_stan = hpv_model_poisson_stan,

                                                  # Number of iterations
                                                  n_iter = n_iter, # main chain
                                                  n_warmup = n_warmup,
                                                  n_chains_mcmc = 4, # Number of chains
                                                  n_iter_sub = n_iter_sub, # Subchain

                                                  # Cut rate
                                                  eta_pois = eta_pois,
                                                  eta_binom = 1,

                                                  out_file_rda = file_i,
                                                  n_cores=n_cores )
    }
    NULL
  }

  return(TRUE)
}


#' @title Compute summary of MCMC results in HPV model
#' @description Compute summary of MCMC used in section 5.4
#' @param mcmc_dir Directory where the MCMC results are read from
#' @param out_dir Directory where outputs are saved
#' @param eta_all Levels of influence of the Poisson module
#' @importFrom foreach foreach %do%
#' @export
smi_sec_5_4_compute_summary <- function(
  mcmc_dir,
  out_dir,
  eta_all
){

  i_eta=1
  hpv_smi_mcmc_all = foreach::foreach( i_eta = seq_along(eta_all), .combine="rbind" ) %do% {
    # i_eta <- 1
    eta_pois <- eta_all[i_eta]
    cat(eta_pois,", ")
    file_i <- paste(out_dir,"/hpv_smi_",formatC(eta_pois,digits=3,format="f",flag="0"),".rda",sep="")
    hpv_mcmc_smi = NULL
    if( file.exists(file_i) ){
      load(file = file_i)
      hpv_mcmc_smi$eta = eta_pois
    }
    hpv_mcmc_smi
  }

  save(hpv_smi_mcmc_all, file=paste(out_dir,"/hpv_smi_mcmc_all.rda",sep=""))

  return(TRUE)

}


#' @title Plot MCMC results in HPV model
#' @description Plot MCMC results from section 5.4
#' @param mcmc_dir Directory with MCMC summaries produced by smi_sec_5_4_compute_summary
#' @param out_dir Directory where outputs are saved
#' @import ggplot2
#' @importFrom dplyr mutate filter select
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom grDevices colorRampPalette
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @export
smi_sec_5_4_plot_mcmc <- function(
  mcmc_dir,
  out_dir
) {

  hpv_smi_mcmc_all <- NULL
  load(file = paste(mcmc_dir, "/hpv_smi_mcmc_all.rda", sep=""))

  # Yellow and black
  col_aux = grDevices::colorRampPalette(c("#000000","#FFD500"))( 2 )

  # Joint posterior of theta1 and theta2 under the full and cut model
  p_theta_joint_cut_full = hpv_smi_mcmc_all %>%
    dplyr::select( one_of(c("eta",'theta1','theta2'))) %>%
    dplyr::filter( .data$eta %in% c(0,1) ) %>%
    dplyr::mutate( eta=as.factor(.data$eta)) %>%
    ggplot()+
    geom_point( aes(x=.data$theta1, y=.data$theta2, col=.data$eta), alpha=0.1 ) +
    coord_cartesian(xlim=c(-3,-1), ylim=c(0,40)) +
    scale_color_manual(values=col_aux, name="eta")+
    guides(colour=guide_legend(override.aes = list(alpha=1)))+
    # labs(title="Posterior distribution of theta",subtitle="Partial Cut method") +
    theme_bw() + theme(legend.position = "bottom")
  ggsave( plot=p_theta_joint_cut_full,
          filename=paste(out_dir,"/hpv_smi_theta_joint_post.pdf",sep=""),
          device="pdf",width=10,height=8,units="cm")

  # Marginal posterior of theta1 and theta2 under SMI for all eta
  aistats2020smi::set_ggtheme()
  p_theta_post = hpv_smi_mcmc_all %>%
    dplyr::select(one_of(c("eta",'theta1','theta2'))) %>%
    tidyr::pivot_longer(c('theta1','theta2'), names_to = "param") %>%
    ggplot( aes(x=.data$value, y=.data$eta, group=.data$eta) ) +
    ggridges::geom_density_ridges(scale=3, alpha=0.5 ) +
    scale_y_continuous(expand = c(0, 0)) +     # will generally have to set the `expand` option
    scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
    facet_wrap( ~.data$param, scales="free_x" ) +
    ggridges::theme_ridges()
  ggsave( plot=p_theta_post,
          filename=paste(out_dir,"/hpv_smi_theta_post.pdf",sep=""),
          device="pdf", width=12,height=8, units="cm")

  return(TRUE)
}


#' @title Compute eta selection in HPV model
#' @description Compute eta selection in HPV model for section 5.4
#' @param mcmc_dir Directory where the MCMC results are read from
#' @param out_dir Directory where outputs are saved
#' @param eta_all Levels of influence of the poisson module
#' @importFrom doRNG %dorng%
#' @importFrom dplyr filter mutate
#' @importFrom foreach foreach
#' @importFrom loo loo waic relative_eff
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyr gather
#' @export
smi_sec_5_4_eta_selection <- function(
  mcmc_dir,
  out_dir,
  eta_all = seq(0.10,1.00,by=0.10)
) {

  hpv_smi_mcmc_all <- NULL
  load(file = paste(mcmc_dir, "/hpv_smi_mcmc_all.rda", sep=""))

  n_obs_hpv = nrow(aistats2020smi::HPV)
  i_eta = 1
  hpv_smi_model_eval <- foreach::foreach( i_eta = seq_along(eta_all), .combine=rbind ) %dorng% {
    loglik_aux <- apply( hpv_smi_mcmc_all %>% dplyr::filter(.data$eta==eta_all[i_eta]), 1,
                         FUN = function(x,HPV) {
                           aistats2020smi::hpv_loglik( data=HPV,
                                       theta=x[paste("theta",1:2,sep="")],
                                       phi=x[paste("phi_",1:n_obs_hpv,sep="")] ) },
                         HPV = aistats2020smi::HPV )
    loglik_aux <- t(loglik_aux); rownames(loglik_aux)<-NULL
    ll_pois <- loglik_aux[,1:n_obs_hpv]
    ll_binom <- loglik_aux[,n_obs_hpv+1:n_obs_hpv]
    waic_aux <- data.frame( poisson = c( loo::waic( ll_pois )$estimates[,1],
                                         loo::loo( ll_pois, r_eff=loo::relative_eff(exp(ll_pois),chain_id=rep(1,nrow(ll_pois))) )$estimates[,1] ),
                            binomial = c( loo::waic( ll_binom )$estimates[,1],
                                          loo::loo( ll_binom, r_eff=loo::relative_eff(exp(ll_binom),chain_id=rep(1,nrow(ll_binom))) )$estimates[,1] )
                            )
    waic_aux <- waic_aux %>%
      dplyr::mutate( score_id = gsub("loo","psis",rownames(waic_aux))) %>%
      tidyr::gather( key="module", value="score", -.data$score_id) %>%
      dplyr::mutate( eta_pois = eta_all[i_eta] )
    waic_aux
  }

  save( hpv_smi_model_eval, file=paste(out_dir,"/hpv_smi_model_eval.rda",sep=""))

  return(TRUE)

}


#' @title Plotting model evaluation results in HPV model
#' @description Plot model evaluation results in HPV model for section 5.4
#' @param mcmc_dir Directory where the MCMC results are read from
#' @param out_dir Directory where outputs are stored
#' @import ggplot2
#' @importFrom dplyr filter mutate
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
smi_sec_5_4_plot_eta_selection <- function(
  mcmc_dir,
  out_dir
) {

  hpv_smi_model_eval <- NULL
  load(paste(mcmc_dir, "/hpv_smi_model_eval.rda", sep=""))

  aistats2020smi::set_ggtheme()
  p_model_eval = hpv_smi_model_eval %>%
    dplyr::mutate( module = dplyr::recode(.data$module, 'poisson'='-elpd poisson', 'binomial'='-elpd binomial') ) %>%
    dplyr::filter( .data$score_id=='elpd_waic', .data$score>-5000) %>%
    ggplot( aes(x=.data$eta_pois, y=-.data$score) ) +
    geom_line( col="red" ) +
    geom_point( col="red" ) +
    facet_wrap( vars(.data$module), ncol=2, scales="free_y" ) +
    theme( axis.title.y=element_blank() ) +
    labs( x="eta (over poisson module)" )
  ggsave( plot=p_model_eval,
          filename=paste(out_dir,"/hpv_smi_model_eval.pdf",sep=""),
          device="pdf", width=12, height=6, units="cm")

  return(TRUE)
}
