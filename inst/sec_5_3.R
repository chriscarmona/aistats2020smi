rm(list = ls())
options(scipen=999, stringsAsFactors=FALSE)
set.seed(0)

# pkgbuild::compile_dll()
# Rcpp::compileAttributes()
# devtools::document()
# devtools::install()

# loading required packages #
req.pck <- c( "aistats2020smi",
              "coda",
              "magrittr","tidyr","dplyr",
              "ggplot2","ggridges","GGally","ggmcmc","ggcorrplot",
              "foreach","doParallel",
              "GPfit" )
req.pck_bool <- sapply(X=req.pck,FUN=require,character.only=T)
if(!all(req.pck_bool)) {
  sapply(X=req.pck[!req.pck_bool],FUN=install.packages,character.only=T);
  sapply(X=req.pck,FUN=require,character.only=T)
}

out_dir <- './inst'

# setwd(out_dir)

# Indicates if the MCMC will be computed (TRUE), or loaded from previously computed results (FALSE)
compute_mcmc = FALSE
n_iter_mcmc = 500e3
n_iter_mcmc_stage2 = 100

# Indicates if it is necessary to thin and warm-up the chain
n_warmup = 20e3
n_thin = 16

# Indicates if the summary of MCMC output will be computed (TRUE), or loaded from previously computed results (FALSE)
compute_summary = FALSE

plot_mcmc = TRUE

compute_model_select = FALSE
plot_model_eval = TRUE

# Parallel processing
parallel_comp = FALSE
if(parallel_comp){
  n_cores = 6
  options(cores=n_cores)
  doParallel::registerDoParallel()
  getDoParWorkers()
}

arc_dataset = c("NMeso")
ManureLevels = c("low","medium","high")

#####################
# SEMI-MODULAR INFERENCE #
#####################

# Levels of influence in SMI
power_eta_all <- c( 0.01,0.99,
                    seq(0.8,0.9,by=0.02),
                    seq(0.10,1.00,by=0.10) )
power_eta_all <- sort( unique(round(power_eta_all,6)) )


if( compute_mcmc ) {
  force_compute_mcmc = FALSE
  foreach( eta_PO_i = seq_along(power_eta_all), .errorhandling="pass" ) %dopar% {
    # eta_PO_i=match(0.82,power_eta_all)

    cat("power_w_PO =", power_eta_all[eta_PO_i] ,"\n")

    ##### LOAD AGRICURB DATA #####
    # Variables that will be transformed in log-scale
    vars_log_i = c("Rainfall")

    # Load data #
    Agricurb_data <- aistats2020smi::get_agricurb_data( arc_datasets="NMeso",
                                                      vars_log=vars_log_i,
                                                      vars_scale_mean_0_var_1=c("Size","Rainfall") )

    # Separates data into modern and archaeological
    data_arc = Agricurb_data %>%
      dplyr::filter(dataset!="modern") %>%
      dplyr::select(c( "Size",
                       "Site",
                       "Rainfall_min",
                       "Rainfall_max",
                       "Category",
                       "normd15N"))
    data_mod = Agricurb_data %>%
      dplyr::filter(dataset=="modern") %>%
      dplyr::select(c( "Site",
                       "ManureLevel",
                       "Rainfall",
                       "Category",
                       "normd15N"))

    # File with MCMC results #
    file_i = paste(out_dir, "/NMeso_smi_bayes_",formatC(power_eta_all[eta_PO_i],digits=3,format="f",flag="0"),".rds",sep="")

    # Compute MCMC If the file doesn't exist or is forced to run
    if( !file.exists(file_i) | force_compute_mcmc ) {
      # Set seed for reproducilibity
      set.seed(0)

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
                                                    power_w_PO = power_eta_all[eta_PO_i],

                                                    # Prior specification: FLat improper priors
                                                    prior_spec_PO="flat",
                                                    prior_spec_HM="flat",


                                                    theta_min_max=theta_min_max,

                                                    # Include Random effects?
                                                    PO_site_rnd_eff = TRUE,
                                                    HM_site_rnd_eff = TRUE,

                                                    # Number of iterations in the main chain
                                                    n_iter = n_iter_mcmc,
                                                    # Number of iterations in the cut secondary chain
                                                    n_iter_sub = n_iter_mcmc_stage2,
                                                    # Number of iterations used for warm-up and thinning
                                                    n_warmup = n_warmup,
                                                    n_thin = n_thin,


                                                    # number of cycles for adapting the proposal distribution of parameters in the PO
                                                    n_epoch_adapt=5,
                                                    n_iter_adapt=2000,

                                                    # Preserve imputed values as an output
                                                    keep_imp = TRUE,

                                                    # Preserve evaluation of the posterior log likelihood
                                                    keep_ll = TRUE,

                                                    # Output file where results are saved
                                                    out_file_rds = file_i )

    }
    NULL
  }
}


if( compute_summary ) {
  agric_smi_summary <- foreach::foreach( eta_PO_i = seq_along(power_eta_all), .errorhandling="pass", .combine="rbind" ) %dopar% {
      # eta_PO_i=match(0.82,power_eta_all)

      file_i = paste(out_dir, "/NMeso_smi_bayes_",formatC(power_eta_all[eta_PO_i],digits=3,format="f",flag="0"),".rds",sep="")
      # file.exists(file_i)

      mcmc_res = NULL
      if( file.exists(file_i) ) {
        mcmc_res <- readRDS( file=file_i )
      }

      # Summary of posteriors
      summary_theta_table <- data.frame(NULL)
      if( !is.null(mcmc_res) & !any(is.na(mcmc_res$theta_mcmc)) ){

        cat("power_w_PO =", power_eta_all[eta_PO_i] ,"\n")

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
        summary_theta_table <- summary_theta_table %>% tidyr::gather("statistic","value",-param)

        # Probability of negative values for gamma_po #
        aux <- apply( mcmc_res$theta_mcmc<=0, 2, mean )
        aux <- data.frame(param=names(aux),statistic="prob_leq_0",value=as.numeric(aux))
        summary_theta_table <- rbind( summary_theta_table,
                                      aux %>% dplyr::filter( stringr::str_detect(param,"^gamma_po_1" ) ) )
        rm(aux)

        # binding posterior summary with predictive performance measures
        if(!is.null(mcmc_res$predict_summary)){
          summary_theta_table <- rbind(summary_theta_table,mcmc_res$predict_summary)
        }

        # Meta data identifying this experiment #
        summary_theta_table$impute_spec = "smi"
        summary_theta_table$imp_playpen = 0
        summary_theta_table$eta_PO = power_eta_all[eta_PO_i]
        summary_theta_table$priorPO = "flat"
        summary_theta_table$priorHM = "flat"
        summary_theta_table$arc_dataset = "NMeso"

      }

      summary_theta_table
    }
  save( agric_smi_summary, file="./data/agric_smi_summary.rda" )

  ### GAMMA credibility intervals ###
  # Joy Division style #
  agric_smi_gamma_post_all = foreach( eta_PO_i = seq_along(power_eta_all), .errorhandling="pass", .combine="rbind" ) %dopar% {
    # dataset_i=1
    # eta_PO_i=match(0.840,power_eta_all)

    file_i = paste(out_dir, "/NMeso_smi_bayes_",formatC(power_eta_all[eta_PO_i],digits=3,format="f",flag="0"),".rds",sep="")
    # file.exists(file_i)

    mcmc_res = NULL
    gamma_post_i = NULL
    if( file.exists(file_i) ) {
      mcmc_res <- readRDS( file=file_i )
      gamma_post_i = data.frame( arc_dataset=arc_dataset[dataset_i],
                                 eta=power_eta_all[eta_PO_i],
                                 mcmc_res[["theta_mcmc"]][,"gamma_po_1",drop=F] )
    }
    gamma_post_i
  }
  agric_smi_gamma_post_all$arc_dataset = factor(agric_smi_gamma_post_all$arc_dataset,levels=c("NMeso"))
  save( agric_smi_gamma_post_all, file="./data/agric_smi_gamma_post_all.rda" )
}

# plotting MCMC results #
if(plot_mcmc) {

  # Loading summary of SMI posteriors
  # load( file="./data/agric_smi_summary.rda" )

  # Setting nice ggplot theme settings
  aistats2020smi::set_ggtheme()

  ### GAMMA credibility intervals ###
  p_gamma_interv <- aistats2020smi::agric_smi_summary %>%
    mutate( arc_dataset=factor(arc_dataset,levels=c("NMeso")) ) %>%
    filter( param=="gamma_po_1" ) %>%
    tidyr::spread(statistic,value) %>%
    ggplot( aes(x=eta_PO,fill=arc_dataset) ) +
    geom_ribbon( aes(ymax=q97.5, ymin=q2.5), alpha=0.25 ) +
    geom_ribbon( aes(ymax=q75, ymin=q25), alpha=0.50 ) +
    geom_line( aes(y=Mean, col=arc_dataset) ) +
    geom_hline(yintercept=0, lty=2) +
    geom_vline( xintercept=0.82, lty=4, col="purple" ) +
    coord_cartesian( xlim=c(0,1) ) +
    theme(legend.position="none") + # Remove legend
    labs( #title="Posterior mean and intervals of gamma",
      x="eta: PO module influence on ManureLevel", y="gamma" )
  # print(p_gamma_interv)
  ggsave( plot=p_gamma_interv,
          filename="agric_smi_gamma_post_interval.pdf",
          device="pdf", width=20, height=8, units="cm" )

  ### GAMMA credibility intervals ###
  # Joy Division style #
  load( file="./data/agric_smi_gamma_post_all.rda" )
  p_gamma_joy = agric_smi_gamma_post_all %>%
    ggplot( aes(x=gamma_po_1,y=eta, group=eta, fill=arc_dataset) ) +
    geom_vline(xintercept=0)+
    ggridges::geom_density_ridges(scale = 4, alpha=0.9) +
    scale_y_continuous(expand = c(0, 0)) +     # will generally have to set the `expand` option
    scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
    coord_cartesian(xlim=c(-4,1)) +
    # coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
    # facet_wrap( ~arc_dataset, scales = "free_x" ) +
    ggridges::theme_ridges() +
    labs( y="eta", x="gamma" ) +
    theme(legend.position="none") # Remove legend
  # print(p_gamma_joy)
  ggsave( plot=p_gamma_joy,
          filename="agric_smi_gamma_post.pdf",
          device="pdf", width=12,height=12, units="cm")

  ### GAMMA BF for negative values ###
  p_gamma_neg <- aistats2020smi::agric_smi_summary %>%
    mutate( arc_dataset=factor(arc_dataset,levels=c("NMeso")) ) %>%
    filter( param=="gamma_po_1",
            statistic=="prob_leq_0") %>%
    mutate( BF=value/(1-value) ) %>%
    select( one_of(c("arc_dataset","eta_PO","BF"))) %>%
    ggplot() +
    geom_line( aes(x=eta_PO, y=BF,col=arc_dataset) ) +
    geom_vline( xintercept=0.82, lty=4, col="purple" ) +
    geom_point( aes(x=eta_PO, y=BF,col=arc_dataset) ) +
    coord_cartesian( xlim=c(0,1) ) +
    labs( x="eta\n(degree of influence of PO module on imputation)", y="Bayes Factor for gamma < 0" )+
    theme(legend.position="none") # Remove legend
  # print(p_gamma_neg)
  ggsave( plot=p_gamma_neg,
          filename="agric_smi_gamma_leq_0_BF.pdf",
          device="pdf", width=12,height=8, units="cm")

}

### Compute model selection ###
if( compute_model_select ) {

  agric_smi_model_eval = foreach::foreach( eta_PO_i = seq_along(power_eta_all), .errorhandling="pass",.combine = rbind ) %dopar% {
      # eta_PO_i=match(0.6,power_eta_all)

      # list.files(out_dir)

      cat("power_w_PO =", power_eta_all[eta_PO_i] ,"\n")

      file_i = paste(out_dir, "/NMeso_smi_bayes_",formatC(power_eta_all[eta_PO_i],digits=3,format="f",flag="0"),".rds",sep="")
      # file.exists(file_i)

      mcmc_res = NULL
      model_eval_i = data.frame(NULL)
      if( file.exists(file_i) ) {
        mcmc_res <- readRDS( file=file_i )
      }

      if( !is.null(mcmc_res) ){
        # elpd for predicting archaelogical data in HM model
        waic_i = loo::waic( mcmc_res$loglik_mcmc[,1:mcmc_res$n_obs_arc,2] )
        loo_i = loo::loo( mcmc_res$loglik_mcmc[,1:mcmc_res$n_obs_arc,2] )

        model_eval_i = data.frame(rbind(waic_i$estimates,loo_i$estimates))
        model_eval_i$stat = rownames(model_eval_i)
        rownames(model_eval_i) = NULL

        model_eval_i$arc_dataset = "NMeso"
        model_eval_i$eta_PO = power_eta_all[eta_PO_i]
        if( (length(dim(model_eval_i))!=2) | ncol(model_eval_i)!=5 ) {
          cat("Problem with power_w_PO =", power_eta_all[eta_PO_i] ,"\n")
          cat(dim(model_eval_i))
          model_eval_i = NULL
        }
      }

      model_eval_i
    }

  agric_smi_model_eval$arc_dataset = factor(agric_smi_model_eval$arc_dataset, levels="NMeso")
  save( agric_smi_model_eval, file="./data/agric_smi_model_eval.rda" )

}

# Select best eta
if(T){
  # Loading elpd estimates for SMI posteriors
  # load( file="./data/agric_smi_model_eval.rda" )

  elpd_approx = c('elpd_waic','elpd_loo')[1]

  aux = aistats2020smi::agric_smi_model_eval %>%
    dplyr::filter(stat==elpd_approx) %>%
    dplyr::filter(arc_dataset=="NMeso")
  gp_elpd = GPfit::GP_fit(X=aux$eta_PO,Y=aux$Estimate)
  elpd_gp_estimate = data.frame( eta_PO=seq(0,1,0.02),
                                 arc_dataset="NMeso",
                                 elpd_hat=predict(gp_elpd, xnew=seq(0,1,0.02))$Y_hat )
  rm(aux,gp_elpd)

  elpd_gp_estimate$arc_dataset = factor(elpd_gp_estimate$arc_dataset,levels="NMeso")

  elpd_gp_best = elpd_gp_estimate %>%
    group_by(arc_dataset) %>%
    filter(elpd_hat == max(elpd_hat)) %>%
    as.data.frame()
}
# elpd_gp_best

# plotting model evaluation results #
if( plot_model_eval ){
  aistats2020smi::set_ggtheme()
  p_elpd <- aistats2020smi::agric_smi_model_eval %>%
    dplyr::filter(stat==elpd_approx) %>%
    ggplot() +
    geom_point( aes(x=eta_PO,y=Estimate, col=arc_dataset) ) +
    geom_line( aes(x=eta_PO,y=elpd_hat,col=arc_dataset), lty=3, data=elpd_gp_estimate ) +
    geom_vline( aes(xintercept=eta_PO), col="purple", lty=2, data=elpd_gp_best ) +
    geom_text(aes(x=elpd_gp_best[,'eta_PO'],y=-min(Estimate),label=elpd_gp_best[,'eta_PO']))+
    theme(legend.position="none") + # Remove legend
    labs(y="-ELPD",x="eta")
  # print(p_elpd)
  ggsave( plot=p_elpd,
          filename=paste("agric_smi_model_eval_",elpd_approx,".pdf",sep=""),
          device="pdf", width=12,height=8, units="cm")


  # MCMC for the best SMI posterior #

  # File with MCMC results #
  eta_star_i = elpd_gp_best[elpd_gp_best$arc_dataset=="NMeso","eta_PO"]
  # file_i = "./data/NMeso_mcmc_smi_best.rda"
  # load( file=file_i )

  # Plot convergence analysis for parameters of interest
  all_par = colnames(aistats2020smi::NMeso_mcmc_smi_best$theta_mcmc)
  interest_par = c("alpha_po_1","alpha_po_2","gamma_po_1","sigma_eta_PO",all_par[substr(all_par,1,7)=='beta_hm'],"sigma_hm","sigma_hm_eta")
  interest_par_tempered = c("alpha_po_1","alpha_po_2","gamma_po_1")

  mcmc_1 = aistats2020smi::NMeso_mcmc_smi_best$theta_mcmc[,interest_par]
  mcmc_2 = aistats2020smi::NMeso_mcmc_smi_best$theta_mcmc_tempered[,interest_par_tempered]; colnames(mcmc_2)=paste(colnames(mcmc_2),'_tilde',sep="")
  mcmc_interest <- coda::as.mcmc( cbind(mcmc_1,mcmc_2) )

  coda::effectiveSize(mcmc_interest)

  model_S <- ggmcmc::ggs( mcmc_interest )
  # Plot MCMC traces for the main parameters in the model
  p_mcmc_trace = ggmcmc::ggs_traceplot(model_S) + labs( title=paste("NMeso eta=",formatC(eta_star_i,digits=3,format="f",flag="0"),sep="") )
  ggsave( plot=p_mcmc_trace,
          filename="agric_model_mcmc_trace.pdf",
          device="pdf", width=30,height=30, units="cm")

  # Autocorrelation
  p_mcmc_autocor = ggmcmc::ggs_autocorrelation(model_S) + labs( title=paste("NMeso eta=",formatC(eta_star_i,digits=3,format="f",flag="0"),sep="") )
  ggsave( plot=p_mcmc_autocor,
          filename="agric_model_mcmc_autocorrelation.pdf",
          device="pdf", width=30,height=20, units="cm")

  # Plot MCMC Cross correlation for the main parameters in the model
  cor_mat = round(cor( mcmc_interest ),2); diag(cor_mat)=NA;
  p_mcmc_crosscor = cor_mat[colnames(mcmc_interest),rev(colnames(mcmc_interest))] %>%
    ggcorrplot::ggcorrplot( method = "circle", type='full', lab = TRUE, lab_size=2) +
    labs( title=paste("NMeso eta=",formatC(eta_star_i,digits=3,format="f",flag="0"),sep="") )
  ggsave( plot=p_mcmc_crosscor,
          filename="agric_model_mcmc_crosscorrelation.pdf",
          device="pdf", width=40,height=20, units="cm")

  # Posterior distribution of gamma
  aistats2020smi::set_ggtheme()
  p <- aistats2020smi::NMeso_mcmc_smi_best[["theta_mcmc"]][,"gamma_po_1",drop=F] %>%
    as.data.frame() %>%
    mutate(param='gamma_po_1') %>%
    ggplot() +
    geom_density( aes(x=gamma_po_1, fill=param), alpha=0.5 ) +
    geom_density( aes(x=gamma_po_1, col=param) ) +
    geom_vline(xintercept=0,lty=2) +
    theme(legend.position="none") + # Remove legend
    labs(x="gamma")
  # print(p)
  ggsave( plot=p,
          filename="gamma_po_1_smi.pdf",
          device="pdf", width=15,height=10, units="cm", dpi="print")

  # Posterior probability of negative gamma #
  p_gamma_leq_0_smi <- mean( aistats2020smi::NMeso_mcmc_smi_best[["theta_mcmc"]][,"gamma_po_1"] < 0 )
  p_gamma_leq_0_smi

  # Bayes factor for negative gamma #
  BF_gamma_leq_0_smi = p_gamma_leq_0_smi / (1-p_gamma_leq_0_smi)
  BF_gamma_leq_0_smi

}
