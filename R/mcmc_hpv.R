#' @title
#'     MCMC procedure for ecological study of HPV and cervical cancer
#'
#' @description
#'    MCMC sampling for the ecological study of HPV and cervical cancer.
#'
#' @param HPV Data frame. Data for correlation between human papillomavirus (HPV) prevalence and cervical cancer incidence.
#' @param hpv_model_whole_stan stanmodel. The compiled stan model with the powered version of the HPV model.
#' @param hpv_model_poisson_stan stanmodel. The compiled stan model with the poisson module of the HPV model.
#' @param eta_pois Float. Degree of influence of the poisson module.
#' @param eta_binom Float. Degree of influence of the poisson module.
#' @param n_iter Integer. Number of iterations in the main MCMC chain.
#' @param n_warmup Integer. Number of updates discarded when "warming-up" the MCMC.
#' @param n_chains_mcmc Integer. Number of parallel chains produced by the sampler
#' @param n_iter_sub Integer. Number of updates in the subchain of the two-stage MCMC
#' @param out_file_rda Indicates a file (.rda) where the output should be saved
#' @param n_cores Number of cores used by stan in the sampling process.
#'
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach
#' @importFrom rstan sampling extract
#'
#' @export
mcmc_hpv <- function( HPV, # HPV data

                      hpv_model_whole_stan,
                      hpv_model_poisson_stan,

                      # SMI degree of influence for each module
                      eta_pois = 1,
                      eta_binom = 1,

                      # Number of iterations
                      n_iter = 10000, # main chain
                      n_warmup = 1000,
                      n_chains_mcmc = 4, # Number of chains
                      n_iter_sub = 200, # Subchain

                      out_file_rda=NULL,
                      n_cores=1 ) {

  # Data in stan format #
  HPV_data_stan <- list( n_obs=nrow(HPV),
                         nhpv = HPV$nhpv,
                         Npart = HPV$Npart,
                         ncases = HPV$ncases,
                         Npop = HPV$Npop,
                         eta_pois=eta_pois,
                         eta_binom=eta_binom,
                         phi=NULL )

  ### Power likelihood ###
  stan_fit <- rstan::sampling( hpv_model_whole_stan,
                               data = HPV_data_stan,
                               iter = n_iter,
                               warmup = n_warmup,
                               chains = n_chains_mcmc,
                               cores = n_cores,
                               # pars=c("phi","theta1","theta2"),
                               show_messages=FALSE )
  hpv_mcmc_pow <- rstan::extract(stan_fit)

  ### Multiple imputation ###
  # Conventional Bayesian inference on poisson model, conditional on imputed phi #
  HPV_data_stan$eta_pois = 1
  HPV_data_stan$eta_binom = 0

  imp_i=1
  hpv_mcmc_stage2 <- foreach::foreach( imp_i = 1:nrow(hpv_mcmc_pow$phi),
                            .combine = rbind ) %dorng% {
                              # imp_i=1

                              # Select one imputed value of phi
                              HPV_data_stan$phi = hpv_mcmc_pow$phi[imp_i, ]

                              # Sample the poisson module conditional on such imputed phi
                              stan_fit <- rstan::sampling( hpv_model_poisson_stan,
                                                           data = HPV_data_stan,
                                                           iter = n_iter_sub,
                                                           chains = 1,
                                                           pars=c("theta1","theta2"),
                                                           show_messages=FALSE )
                              stan_fit_mcmc <- rstan::extract(stan_fit)

                              # Return only the last sample
                              n_aux <- length(stan_fit_mcmc$theta1)
                              c(stan_fit_mcmc$theta1[n_aux], stan_fit_mcmc$theta2[n_aux])
                            }

  colnames(hpv_mcmc_pow$phi) <- paste("phi_", 1:nrow(HPV), sep = "")
  colnames(hpv_mcmc_stage2) <- paste("theta", 1:2, sep = "")

  hpv_mcmc_smi <- data.frame( hpv_mcmc_pow$phi,
                              hpv_mcmc_stage2,
                              theta_tilde_1 = hpv_mcmc_pow$theta1,
                              theta_tilde_2 = hpv_mcmc_pow$theta2 )
  rownames(hpv_mcmc_smi) = NULL

  # Save results #
  if( !is.null(out_file_rda) ){
    save( hpv_mcmc_smi, file=out_file_rda )
  }

  return( hpv_mcmc_smi )

}
