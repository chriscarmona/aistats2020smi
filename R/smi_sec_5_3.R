
#' @title Compute MCMC in agricultural model
#' @description Compute MCMC used in section 5.3
#' @param out_dir Directory where outputs are stored
#' @param eta_all Levels of influence of the PO module
#' @param force_compute_mcmc To avoid computation time, completed MCMC is not recomputed, unless force_compute_mcmc is True
#' @param n_iter Number of iterations in the main chain
#' @param n_iter_sub Number of iterations in the sub-chain for the two-stage MCMC
#' @param n_warmup Number of iterations used for warm-up
#' @param n_thin Number of iterations used for thinning
#' @param n_epoch_adapt Number of epochs for adaptation of proposal distributions
#' @param n_iter_adapt Number of iterations in each adaptation epoch
#' @param seed Set random seed before MCMC starts for reproducibility
#' @importFrom foreach foreach %dopar%
#' @importFrom dplyr filter select
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr str_detect
#' @export
smi_sec_5_3_compute_mcmc <- function(
    out_dir,
    eta_all = seq(0.10,1.00,by=0.10),
    force_compute_mcmc = FALSE,
    n_iter = 100e3,
    n_iter_sub = 100,
    n_warmup = 10e3,
    n_thin = 3,
    n_epoch_adapt = 5,
    n_iter_adapt = 1000,
    seed = 123
  ) {

  i_eta=1
  out = foreach::foreach( i_eta = seq_along(eta_all), .errorhandling="pass" ) %dopar% {
    # i_eta=match(0.8,eta_all)

    cat("eta =", eta_all[i_eta] ,"\n")

    ##### LOAD AGRICURB DATA #####
    # Variables that will be transformed in log-scale
    vars_log_i = c("Rainfall")

    # Load data #
    Agricurb_data <- aistats2020smi::get_agricurb_data( arc_datasets="NMeso",
                                                        vars_log=vars_log_i,
                                                        vars_scale_mean_0_var_1=c("Size","Rainfall") )

    # Separates data into modern and archaeological
    data_arc = Agricurb_data %>%
      dplyr::filter( .data$dataset!="modern") %>%
      dplyr::select(c( "Size",
                       "Site",
                       "Rainfall_min",
                       "Rainfall_max",
                       "Category",
                       "normd15N"))
    data_mod = Agricurb_data %>%
      dplyr::filter( .data$dataset=="modern") %>%
      dplyr::select(c( "Site",
                       "ManureLevel",
                       "Rainfall",
                       "Category",
                       "normd15N"))

    # File with MCMC results #
    file_i = paste(out_dir, "/NMeso_smi_bayes_",formatC(eta_all[i_eta],digits=3,format="f",flag="0"),".rda",sep="")

    # Compute MCMC If the file doesn't exist or is forced to run
    if( !file.exists(file_i) | force_compute_mcmc ) {
      # Set seed for reproducilibity
      set.seed(seed)

      # Set limits for prior distributions
      Sites_arc = sort(unique(data_arc$Site))
      theta_min_max <- matrix(NA,4+length(Sites_arc),2)
      rownames(theta_min_max) <- c("alpha_po_1","alpha_po_2","gamma_po_1","sigma_eta_PO",paste("eta_po_",Sites_arc,sep=""))
      theta_min_max[,1] = -5
      theta_min_max[,2] = 5
      theta_min_max[rownames(theta_min_max)=="sigma_eta_PO",1]=0

      # MCMC #
      mcmc_res <- aistats2020smi::mcmc_agric_model( data_arc = data_arc,
                                                    data_mod = data_mod,

                                                    # Type of imputation: SMI
                                                    impute_spec = "smi",
                                                    power_w_PO = eta_all[i_eta],

                                                    # Prior specification: FLat improper priors
                                                    prior_spec_PO="flat",
                                                    prior_spec_HM="flat",


                                                    theta_min_max=theta_min_max,

                                                    # Include Random effects?
                                                    PO_site_rnd_eff = TRUE,
                                                    HM_site_rnd_eff = TRUE,

                                                    # Number of iterations in the main chain
                                                    n_iter = n_iter,
                                                    # Number of iterations in the cut secondary chain
                                                    n_iter_sub = n_iter_sub,
                                                    # Number of iterations used for warm-up and thinning
                                                    n_warmup = n_warmup,
                                                    n_thin = n_thin,

                                                    # number of cycles for adapting the proposal distribution of parameters in the PO
                                                    n_epoch_adapt = n_epoch_adapt,
                                                    n_iter_adapt = n_iter_adapt,

                                                    # Preserve imputed values as an output
                                                    keep_imp = TRUE,

                                                    # Preserve evaluation of the posterior log likelihood
                                                    keep_ll = TRUE,

                                                    # Output file where results are saved
                                                    out_file_rda = file_i )

    }
    TRUE
  }

  any_ok = FALSE
  for( out_i in out){
    if( class(out_i)=='logical' & out_i ) {
      any_ok = TRUE
    }
  }
  if(!any_ok){
    stop("All MCMC runs failed on smi_sec_5_3_compute_mcmc.")
  }

  return(TRUE)
}



#' @title Compute summary of MCMC results in agricultural model
#' @description Compute summary of MCMC used in section 5.3
#' @param mcmc_dir Directory where the MCMC results are read from
#' @param out_dir Directory where outputs are saved
#' @param eta_all Levels of influence of the PO module
#' @importFrom foreach foreach %dopar%
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
smi_sec_5_3_compute_summary <- function(
  mcmc_dir,
  out_dir,
  eta_all = seq(0.10,1.00,by=0.10)
  ) {

  i_eta=1
  agric_smi_summary <- foreach::foreach( i_eta = seq_along(eta_all), .errorhandling="pass", .combine="rbind" ) %dopar% {
      # i_eta=match(0.82,eta_all)

      file_i = paste(mcmc_dir, "/NMeso_smi_bayes_",formatC(eta_all[i_eta],digits=3,format="f",flag="0"),".rda",sep="")
      # file.exists(file_i)

      mcmc_res = NULL
      if( file.exists(file_i) ) {
        load( file=file_i )
      }

      # Summary of posteriors
      summary_theta_table <- data.frame(NULL)
      if( !is.null(mcmc_res) & !any(is.na(mcmc_res$theta_mcmc)) ){

        cat("power_w_PO =", eta_all[i_eta] ,"\n")

        # Summary of posteriors
        aux <- as.data.frame(mcmc_res$theta_mcmc)
        aux <-as.matrix(aux)

        summary_theta <- summary(mcmc_res$theta_mcmc)
        aux <- colnames(summary_theta[["quantiles"]])
        colnames(summary_theta[["quantiles"]]) <- paste("q",substr(aux,1,nchar(aux)-1),sep="")
        rm(aux)
        summary_theta_table <- data.frame( param=rownames(summary_theta[["statistics"]]),
                                           summary_theta[["statistics"]], summary_theta[["quantiles"]])
        rownames(summary_theta_table)<-NULL
        summary_theta_table <- summary_theta_table %>%
          tidyr::gather("statistic","value",-.data$param)

        # Probability of negative values for gamma_po #
        aux <- apply( mcmc_res$theta_mcmc<=0, 2, mean )
        aux <- data.frame(param=names(aux), statistic="prob_leq_0", value=as.numeric(aux))
        summary_theta_table <- rbind( summary_theta_table,
                                      aux %>% dplyr::filter( stringr::str_detect(.data$param,"^gamma_po_1" ) ) )
        rm(aux)

        # binding posterior summary with predictive performance measures
        if(!is.null(mcmc_res$predict_summary)){
          summary_theta_table <- rbind(summary_theta_table,mcmc_res$predict_summary)
        }

        # Meta data identifying this experiment #
        summary_theta_table$impute_spec = "smi"
        summary_theta_table$imp_playpen = 0
        summary_theta_table$eta_PO = eta_all[i_eta]
        summary_theta_table$priorPO = "flat"
        summary_theta_table$priorHM = "flat"
        summary_theta_table$arc_dataset = "NMeso"

      }

      summary_theta_table
    }

  save( agric_smi_summary, file=paste(out_dir,"/agric_smi_summary.rda",sep="") )

  ### GAMMA credibility intervals ###
  # Joy Division style #
  i_eta=1
  agric_smi_post_all = foreach( i_eta = seq_along(eta_all), .errorhandling="pass", .combine="rbind" ) %dopar% {
    # i_eta=match(0.840,eta_all)

    file_i = paste(mcmc_dir, "/NMeso_smi_bayes_",formatC(eta_all[i_eta],digits=3,format="f",flag="0"),".rda",sep="")
    # file.exists(file_i)

    PO_params = c("gamma_po_1", "alpha_po_1", "alpha_po_2", "sigma_eta_PO")

    mcmc_res = NULL
    post_i = NULL
    if( file.exists(file_i) ) {
      load( file=file_i ) # loads mcmc_res
      post_i = data.frame( arc_dataset="NMeso",
                           eta=eta_all[i_eta],
                           mcmc_res[["theta_mcmc"]][,PO_params,drop=F] )
    }
    post_i
  }

  agric_smi_post_all$arc_dataset = factor(agric_smi_post_all$arc_dataset,levels=c("NMeso"))
  save( agric_smi_post_all, file=paste(out_dir,"/agric_smi_post_all.rda",sep="") )

  return(TRUE)
}



#' @title Plot MCMC results in agricultural model
#' @description Plot MCMC results from section 5.3
#' @param mcmc_dir Directory with MCMC summaries produced by smi_sec_5_3_compute_summary
#' @param out_dir Directory where outputs are saved
#' @import ggplot2
#' @importFrom dplyr mutate filter one_of select
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyr spread
#' @export
smi_sec_5_3_plot_mcmc <- function(
  mcmc_dir,
  out_dir
) {

  # Loading summary of SMI posteriors
  agric_smi_summary <- agric_smi_post_all <- NULL
  load( file=paste(mcmc_dir,"/agric_smi_summary.rda",sep="") )
  load( file=paste(mcmc_dir,"/agric_smi_post_all.rda",sep="") )

  # Setting nice ggplot theme settings
  aistats2020smi::set_ggtheme()

  ### GAMMA credibility intervals ###
  p_gamma_interv <- agric_smi_summary %>%
    dplyr::mutate( arc_dataset=factor(.data$arc_dataset, levels=c("NMeso")) ) %>%
    dplyr::filter( .data$param=="gamma_po_1" ) %>%
    tidyr::spread(.data$statistic, .data$value) %>%
    ggplot( aes(x=.data$eta_PO, fill=.data$arc_dataset) ) +
    geom_ribbon( aes(ymax=.data$q97.5, ymin=.data$q2.5), alpha=0.25 ) +
    geom_ribbon( aes(ymax=.data$q75, ymin=.data$q25), alpha=0.50 ) +
    geom_line( aes(y=.data$Mean, col=.data$arc_dataset) ) +
    geom_hline( yintercept=0, lty=2) +
    geom_vline( xintercept=0.82, lty=4, col="purple" ) +
    coord_cartesian( xlim=c(0,1) ) +
    theme(legend.position="none") + # Remove legend
    labs( #title="Posterior mean and intervals of gamma",
      x="eta: PO module influence on ManureLevel", y="gamma" )
  # print(p_gamma_interv)
  ggsave( plot=p_gamma_interv,
          filename=paste(out_dir,"/agric_smi_gamma_post_interval.pdf",sep=""),
          device="pdf", width=20, height=8, units="cm" )


  ### SMI posterior ###
  # Joy Division style #

  # gamma credibility intervals
  param_i = "gamma_po_1"

  # Joy Division style #
  p_post_joy = agric_smi_post_all %>%
    dplyr::select(c("arc_dataset","eta",all_of(param_i))) %>% #head()
    'colnames<-'(value=c("arc_dataset","eta","value")) %>%
    ggplot( aes(x=.data$value,y=.data$eta, group=.data$eta, fill=.data$arc_dataset) ) +
    geom_vline(xintercept=0)+
    ggridges::geom_density_ridges(scale=4, alpha=0.9) +
    scale_y_continuous(expand = c(0, 0)) +     # will generally have to set the `expand` option
    scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
    # coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
    coord_cartesian(xlim = c(-4,1)) +
    # coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
    # facet_wrap( ~arc_dataset, scales = "free_x" ) +
    ggridges::theme_ridges() +
    labs( y="eta", x=param_i ) +
    theme(legend.position="none") # Remove legend
  # print(p_post_joy)
  ggsave( plot=p_post_joy,
          filename=paste("agric_smi_post_",param_i,".pdf",sep=""),
          device="pdf", width=13,height=8, units="cm")

  ### GAMMA BF for negative values ###
  p_gamma_neg <- agric_smi_summary %>%
    dplyr::mutate( arc_dataset=factor(.data$arc_dataset,levels=c("NMeso")) ) %>%
    dplyr::filter( .data$param=="gamma_po_1",
                   .data$statistic=="prob_leq_0") %>%
    dplyr::mutate( BF=.data$value/(1-.data$value) ) %>%
    dplyr::select( dplyr::one_of(c("arc_dataset","eta_PO","BF"))) %>%
    ggplot() +
    geom_line( aes(x=.data$eta_PO, y=.data$BF, col=.data$arc_dataset) ) +
    geom_vline( xintercept=0.82, lty=4, col="purple" ) +
    geom_point( aes(x=.data$eta_PO, y=.data$BF,col=.data$arc_dataset) ) +
    coord_cartesian( xlim=c(0,1) ) +
    labs( x="eta\n(degree of influence of PO module on imputation)", y="Bayes Factor for gamma < 0" )+
    theme(legend.position="none") # Remove legend
  # print(p_gamma_neg)
  ggsave( plot=p_gamma_neg,
          filename=paste(out_dir,"/agric_smi_gamma_leq_0_BF.pdf",sep=""),
          device="pdf", width=12,height=8, units="cm")

  return( TRUE )
}



#' @title Compute eta selection in agricultural model
#' @description Compute eta selection in agricultural model for section 5.3
#' @param mcmc_dir Directory where the MCMC results are read from
#' @param out_dir Directory where outputs are saved
#' @param eta_all Levels of influence of the PO module
#' @importFrom dplyr filter group_by
#' @importFrom foreach foreach %dopar%
#' @importFrom GPfit GP_fit predict.GP
#' @importFrom loo loo waic
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
smi_sec_5_3_eta_selection <- function(
  mcmc_dir,
  out_dir,
  eta_all = seq(0.10,1.00,by=0.10)
) {

  i_eta=1
  agric_smi_model_eval = foreach::foreach( i_eta = seq_along(eta_all), .errorhandling="pass",.combine = rbind ) %dopar% {
      # i_eta=match(0.6,eta_all)

      cat("eta_PO =", eta_all[i_eta] ,"\n")

      file_i = paste(mcmc_dir, "/NMeso_smi_bayes_",formatC(eta_all[i_eta],digits=3,format="f",flag="0"),".rda",sep="")
      # file.exists(file_i)

      mcmc_res = NULL
      model_eval_i = data.frame(NULL)
      if( file.exists(file_i) ) {
        load( file=file_i )
      }

      if( !is.null(mcmc_res) ){
        # elpd for predicting archaelogical data in HM model
        waic_i = loo::waic( mcmc_res$loglik_mcmc[,1:mcmc_res$n_obs_arc,2] )
        loo_i = loo::loo( mcmc_res$loglik_mcmc[,1:mcmc_res$n_obs_arc,2] )

        model_eval_i = data.frame(rbind(waic_i$estimates,loo_i$estimates))
        model_eval_i$stat = rownames(model_eval_i)
        rownames(model_eval_i) = NULL

        model_eval_i$arc_dataset = "NMeso"
        model_eval_i$eta_PO = eta_all[i_eta]
        if( (length(dim(model_eval_i))!=2) | ncol(model_eval_i)!=5 ) {
          cat("Problem with power_w_PO =", eta_all[i_eta] ,"\n")
          cat(dim(model_eval_i))
          model_eval_i = NULL
        }
      }

      model_eval_i
    }

  agric_smi_model_eval$arc_dataset = factor(agric_smi_model_eval$arc_dataset, levels="NMeso")
  save( agric_smi_model_eval, file=paste(out_dir,"/agric_smi_model_eval.rda",sep="") )

  # Select best eta #
  elpd_approx = c('elpd_waic','elpd_loo')[1]

  aux = agric_smi_model_eval %>%
    dplyr::filter( .data$stat==elpd_approx ) %>%
    dplyr::filter( .data$arc_dataset=="NMeso")
  gp_elpd = GPfit::GP_fit(X=aux$eta_PO,Y=aux$Estimate)
  elpd_gp_estimate = data.frame( eta_PO=seq(0,1,0.02),
                                 arc_dataset="NMeso",
                                 elpd_hat= GPfit::predict.GP(.data$gp_elpd, xnew=seq(0,1,0.02))$Y_hat )
  rm(aux,gp_elpd)

  elpd_gp_estimate$arc_dataset = factor(elpd_gp_estimate$arc_dataset, levels="NMeso")

  save( elpd_gp_estimate, file=paste(out_dir,"/elpd_gp_estimate.rda",sep="") )

  elpd_gp_best = elpd_gp_estimate %>%
    dplyr::group_by( .data$arc_dataset) %>%
    dplyr::filter( .data$elpd_hat == max(.data$elpd_hat)) %>%
    as.data.frame()

  # Duplicate best MCMC and save it as NMeso_mcmc_smi_best
  file_best = paste(mcmc_dir, "/NMeso_smi_bayes_",formatC(elpd_gp_best[1,'eta_PO'],digits=3,format="f",flag="0"),".rda",sep="")
  file.copy( from=file_best, to=paste(out_dir,"/NMeso_mcmc_smi_best.rda",sep="") )

  return(TRUE)
}


#' @title Plotting model evaluation results in agricultural model
#' @description Plot model evaluation results in agricultural model for section 5.3
#' @param mcmc_dir Directory where the MCMC results are read from
#' @param out_dir Directory where outputs are stored
#' @param extra_plots IF TRUE, Additional plots for MCMC evaluation are generated
#' @import ggplot2
#' @importFrom coda effectiveSize
#' @importFrom dplyr filter mutate one_of
#' @importFrom ggcorrplot ggcorrplot
#' @importFrom ggmcmc ggs ggs_traceplot ggs_autocorrelation
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats cor
#' @export
smi_sec_5_3_plot_eta_selection <- function(
  mcmc_dir,
  out_dir,
  extra_plots = FALSE
) {

  agric_smi_summary <- agric_smi_model_eval <- elpd_gp_estimate <- NMeso_mcmc_smi_best <- NULL
  load( file=paste(mcmc_dir,"/agric_smi_summary.rda",sep="") )
  load( file=paste(mcmc_dir,"/agric_smi_model_eval.rda",sep="") )
  load( file=paste(mcmc_dir,"/elpd_gp_estimate.rda",sep="") )
  load( file=paste(mcmc_dir,"/NMeso_mcmc_smi_best.rda",sep="") )

  elpd_approx = c('elpd_waic','elpd_loo')[1]

  elpd_gp_best = elpd_gp_estimate %>%
    dplyr::group_by(.data$arc_dataset) %>%
    dplyr::filter( .data$elpd_hat == max(.data$elpd_hat)) %>%
    as.data.frame()

  aistats2020smi::set_ggtheme()
  p_elpd <- agric_smi_model_eval %>%
    dplyr::filter( .data$stat==elpd_approx) %>%
    ggplot() +
    geom_point( aes(x=.data$eta_PO,y=-.data$Estimate, col=.data$arc_dataset) ) +
    geom_line( aes(x=.data$eta_PO,y=-.data$elpd_hat,col=.data$arc_dataset), lty=3, data=elpd_gp_estimate ) +
    geom_vline( aes(xintercept=.data$eta_PO), col="purple", lty=2, data=elpd_gp_best ) +
    geom_text( aes(x=elpd_gp_best[,'eta_PO'], y=-min(.data$Estimate), label=elpd_gp_best[,'eta_PO']))+
    theme(legend.position="none") + # Remove legend
    labs(y="-ELPD",x="eta")
  # print(p_elpd)
  ggsave( plot=p_elpd,
          filename=paste(out_dir,"/agric_smi_model_eval_",elpd_approx,".pdf",sep=""),
          device="pdf", width=12,height=8, units="cm")


  # MCMC for the best SMI posterior #

  # File with MCMC results #
  eta_star_i = elpd_gp_best[elpd_gp_best$arc_dataset=="NMeso","eta_PO"]

  # Plot convergence analysis for parameters of interest
  all_par = colnames(NMeso_mcmc_smi_best$theta_mcmc)
  interest_par = c("alpha_po_1","alpha_po_2","gamma_po_1","sigma_eta_PO",all_par[substr(all_par,1,7)=='beta_hm'],"sigma_hm","sigma_hm_eta")
  interest_par_tempered = c("alpha_po_1","alpha_po_2","gamma_po_1")

  mcmc_1 = NMeso_mcmc_smi_best$theta_mcmc[,interest_par]
  mcmc_2 = NMeso_mcmc_smi_best$theta_mcmc_tempered[,interest_par_tempered]; colnames(mcmc_2)=paste(colnames(mcmc_2),'_tilde',sep="")
  mcmc_interest <- coda::as.mcmc( cbind(mcmc_1,mcmc_2) )

  if(extra_plots) {
    print(coda::effectiveSize(mcmc_interest))

    model_S <- ggmcmc::ggs( mcmc_interest )
    # Plot MCMC traces for the main parameters in the model
    p_mcmc_trace = ggmcmc::ggs_traceplot(model_S) + labs( title=paste("NMeso eta=",formatC(eta_star_i,digits=3,format="f",flag="0"),sep="") )
    ggsave( plot=p_mcmc_trace,
            filename=paste(out_dir,"/agric_model_mcmc_trace.pdf",sep=""),
            device="pdf", width=30,height=30, units="cm")

    # Autocorrelation
    p_mcmc_autocor = ggmcmc::ggs_autocorrelation(model_S) +
      labs( title=paste("NMeso eta=",formatC(eta_star_i,digits=3,format="f",flag="0"),sep="") )
    ggsave( plot=p_mcmc_autocor,
            filename=paste(out_dir,"/agric_model_mcmc_autocorrelation.pdf",sep=""),
            device="pdf", width=30,height=20, units="cm")

    # Plot MCMC Cross correlation for the main parameters in the model
    cor_mat = round(stats::cor( mcmc_interest ),2); diag(cor_mat)=NA;
    p_mcmc_crosscor = cor_mat[colnames(mcmc_interest),rev(colnames(mcmc_interest))] %>%
      ggcorrplot::ggcorrplot( method = "circle", type='full', lab = TRUE, lab_size=2) +
      labs( title=paste("NMeso eta=",formatC(eta_star_i, digits=3, format="f",flag="0"), sep="") )
    ggsave( plot=p_mcmc_crosscor,
            filename=paste(out_dir,"/agric_model_mcmc_crosscorrelation.pdf",sep=""),
            device="pdf", width=40,height=20, units="cm")


    # Posterior distribution of gamma
    aistats2020smi::set_ggtheme()
    p <- NMeso_mcmc_smi_best[["theta_mcmc"]][,"gamma_po_1",drop=F] %>%
      as.data.frame() %>%
      dplyr::mutate(param='gamma_po_1') %>%
      ggplot() +
      geom_density( aes(x=.data$gamma_po_1, fill=.data$param), alpha=0.5 ) +
      geom_density( aes(x=.data$gamma_po_1, col=.data$param) ) +
      geom_vline(xintercept=0, lty=2) +
      theme(legend.position="none") + # Remove legend
      labs(x="gamma")
    # print(p)
    ggsave( plot=p,
            filename=paste(out_dir,"/gamma_po_1_smi.pdf",sep=""),
            device="pdf", width=15,height=10, units="cm", dpi="print")

  }

  # Posterior probability of negative gamma #
  # p_gamma_leq_0_smi <- mean( NMeso_mcmc_smi_best[["theta_mcmc"]][,"gamma_po_1"] < 0 )
  # print(p_gamma_leq_0_smi)

  # Bayes factor for negative gamma #
  # BF_gamma_leq_0_smi = p_gamma_leq_0_smi / (1-p_gamma_leq_0_smi)
  # print(BF_gamma_leq_0_smi)

  ### GAMMA BF for negative values ###
  BF_data = agric_smi_summary %>%
    dplyr::mutate( arc_dataset=factor(.data$arc_dataset, levels=c("NMeso")) ) %>%
    dplyr::filter( .data$param=="gamma_po_1",
                   .data$statistic=="prob_leq_0") %>%
    dplyr::mutate( value=.data$value/(1-.data$value) ) %>%
    dplyr::select( dplyr::one_of(c("arc_dataset", "eta_PO", "value")) ) %>%
    dplyr::mutate( stat='BF( gamma<0 )', value_hat=.data$value )

  ELPD_data = agric_smi_model_eval %>%
    dplyr::filter( .data$stat=='elpd_waic') %>%
    dplyr::mutate( value_hat=-.data$Estimate, stat="-elpd( Z )") %>%
    dplyr::select( dplyr::one_of(c('arc_dataset', 'eta_PO','value_hat', 'stat')) ) %>%
    merge( elpd_gp_estimate %>% dplyr::rename(value=.data$elpd_hat) %>% dplyr::mutate(value=-.data$value) )

  p_elpd <- rbind(ELPD_data,BF_data) %>%
    ggplot() +
    geom_line( aes( x=.data$eta_PO, y=.data$value, col=.data$arc_dataset) ) +
    geom_point( aes( x=.data$eta_PO, y=.data$value_hat, col=.data$arc_dataset) ) +
    facet_wrap( vars(.data$stat), ncol=1, scales = "free_y" ) +
    geom_vline( aes(xintercept=.data$eta_PO), col="purple", lty=2, data=elpd_gp_best ) +
    geom_text( aes(x=.data$eta_PO,y=-.data$elpd_hat,label=.data$eta_PO), hjust="outward", vjust="inward",
               data=elpd_gp_best%>%dplyr::mutate(stat="-elpd( Z )")) +
    theme( axis.title.y=element_blank(), # Remove y label
           legend.position="none" ) # Remove legend
  # print(p_elpd)
  ggsave( plot=p_elpd,
          filename=paste(out_dir,"/agric_smi_model_eval.pdf",sep=""),
          device="pdf", width=10,height=8, units="cm")

  return(TRUE)

}
