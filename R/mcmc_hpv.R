#' @title
#'     MCMC procedure for ecological study of HPV and cervical cancer
#'
#' @description
#'     This function perform MCMC sampling for the ecological study of HPV and cervical cancer.
#'     There are two modules in this model: Poisson and Binomial.
#'
#' @param HPV Data frame. Data for correlation between human papillomavirus (HPV) prevalence and cervical cancer incidence.
#'
#' @param impute_spec Character. Specification of imputed values of ManureLevel: "full", "cut", "smi"
#'
#' @param power_w_PO Numeric. Raise the likelihood in the PO module to a power when performing M-H steps.
#' @param power_w_HM Numeric. Raise the likelihood in the HM module to a power when performing M-H steps.
#'
#' @param prior_spec_PO specification of prior distributions for the parameters in the PO model: "flat" or "proper"
#' @param prior_spec_HM specification of prior distributions for the parameters in the HM model: "flat" or "proper"
#'
#' @param PO_site_rnd_eff Boolean. Shall the PO module use random effects by Site
#' @param HM_site_rnd_eff Boolean. Shall the HM module use random effects by Site
#'
#' @param n_iter Integer. Number of iterations in the main MCMC chain.
#' @param n_iter_sub Integer. Number of updates in the subchain for parameters in the PO module.
#'
#' @param theta_ini Numeric vector. Initial values for the parameters
#' @param theta_min_max matrix with two columns, minimum and maximum values for each parameter
#' @param theta_prop_int Used to control the width of the proposal distribution for parameters in PO and HM modules
#' @param theta_prop_kernel Shape of the proposal distribution: "uniform" or "normal"
#'
#' @param n_epoch_adapt Integer. Number of epochs that adaptation runs of the MCMC, to addapt the proposal distribution, is performed before the real chain.
#' @param n_iter_adapt Integer. Number of iterations in the adaptation runs of the MCMC.
#'
#' @param ManureLevel_arc_imp_ini Initial values for imputing the missing values of Rainfall
#' @param Rainfall_arc_imp_ini Initial values for imputing the missing values of Rainfall
#'
#' @param gibbs_hm Boolean. Shall we use Gibss to update parameters in the HM module. If FALSE, M_H is used.
#'
#' @param PO_expand Indicates if the Proportional Odds models should be expanded by considering an additional mixture component
#' @param POmixda Indicates if the inference for the mixture model uses data augmentation
#' @param lambda_prop_int Double. Width of the proposal distribution for the mixing weight, when PO_expand=TRUE
#' @param lambda_mix_PO_prior Numeric vector. indicates the two parameters of the beta prior for the mixture weight
#'
#' @param keep_imp Indicates if the imputed values for missing data should be returned
#'
#' @param out_file_rds Indicates a file (.rds) where the output should be saved
#' @param log_file Indicates a file (.txt) where the log of the process should be saved
#' @param devel Development mode
#'
#'
#' @details
#' The hierarchical model consists of two parts: The Proportional Odds (PO1) component, and the Gaussian Linear model (HM1).
#'   HM1:
#'      normd15N ∼ 1 + Rainfall + ManureLevel + (1|Site)
#'      weights = varIdent(form=~1|Category)
#'   PO1:
#'      ManureLevel ∼ 1 + Size + (1|Site)
#'      link = "logistic"
#'
#' MCMC details:
#'     Coefficients in the HM1 modulel can be updated using Gibbs sampling (gibbs_hm=TRUE) for the joint conditional posterior, M-H otherwise
#'     Parameters in the PO module are updated using M-H, updating one by one.
#'     ManureLevel missing values is updated using M-H one value at a time.
#'     Rainfall missing values are updated all together using M-H.
#'
#'     This function allows to perform adaptations for the proposal of the PO1 parameters. Change n_epoch_adapt and n_iter_adapt.
#'
#'     In this implementation, we allow to perform several types of inference.
#'     1) Conventional Bayes (default): impute_spec="full", set power_w_HM=1, power_w_PO=1
#'     2) Powered likelihood: impute_spec="full", gibbs_hm="FALSE", set "power_w_HM" and "power_w_PO" to control the influence of each module in the update of parameters and missing data.
#'     3) Cut model: impute_spec="cut" (deprecate set power_w_HM and power_w_PO). Bayesian multiple imputation for ManureLevel.
#'     4) smi imputation: impute_spec="smi", set "power_w_HM" and "power_w_PO" to control the influence of each module in the imputation of ManureLevel.
#'
#'     Parameters are initialized in the MLE using a single imputation of missing values.
#'
#' @return
#'
#' A list with three main elements, described below, and some details about the run.
#' \describe{
#'     \item{\code{theta_mcmc}}{Matrix with the chains of the parameters in the model.}
#'     \item{\code{ManureLevel_imp_mcmc}}{if keep_imp=TRUE, matrix with the chains of the imputed values of the missing ManureLevel.}
#'     \item{\code{Rainfall_imp_mcmc}}{if keep_imp=TRUE, matrix with the chains of the imputed values of the missing Rainfall.}
#' }
#'
#'
#' @examples
#'
#' \dontrun{
#'
#' ##### Cut model for NMeso data #####
#'
#' }
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom stats runif
#' @import nlme
#' @import ordinal
#'
#' @export
mcmc_hpv <- function( HPV, # HPV data

                      # SMI degree of influence for each module
                      eta_pois = 1,
                      eta_binom = 1,

                      # Number of iterations
                      n_iter_mcmc = 10000, # main chain
                      n_iter_warmup = 1000,
                      n_chains_mcmc = 4, # Number of chains
                      n_iter_mcmc_stage_2 = 200, # Subchain

                      # Other
                      mcmc_file=NULL,
                      n_cores ) {
  browser()
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
  stan_fit <- rstan::sampling( stanmodels$hpv_model_full_pow,
                               data = HPV_data_stan,
                               iter = n_iter_mcmc,
                               warmup = n_iter_warmup,
                               chains = n_chains_mcmc,
                               cores = n_cores,
                               # pars=c("phi","theta1","theta2"),
                               show_messages=FALSE )
  hpv_mcmc_pow <- rstan::extract(stan_fit)

  ### Multiple imputation ###
  # Conventional Bayesian inference on poisson model, conditional on imputed phi #
  HPV_data_stan$eta_pois = 1
  HPV_data_stan$eta_binom = 0

  hpv_mcmc_smi_stage2 <- foreach::foreach( imp_i = 1:nrow(hpv_mcmc_pow$phi),
                            .combine = rbind ) %dorng% {
                              # imp_i=1

                              # Select one imputed value of phi
                              HPV_data_stan$phi = hpv_mcmc_pow$phi[imp_i, ]

                              # Sample the poisson module conditional on such imputed phi
                              stan_fit <- rstan::sampling( stanmodels$hpv_model_pois_module,
                                                           data = HPV_data_stan,
                                                           iter = n_iter_mcmc_stage_2,
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
  if( !is.null(mcmc_file) ){
    saveRDS( hpv_mcmc_smi, file=mcmc_file )
  }

  return( hpv_mcmc_smi )

}
