rm(list = ls())
options(scipen=999, stringsAsFactors=FALSE)
set.seed(0)

# pkgbuild::compile_dll()
# Rcpp::compileAttributes()
# devtools::document()
# devtools::install()

# loading required packages #
req.pck <- c( "aistats2020smi",
              "MASS","nlme","coda","ordinal","LearnBayes",
              "plyr","tidyverse", "ggplot2","GGally","lattice",
              "foreach","doParallel",
              "GPfit","rBayesianOptimization",
              "RColorBrewer","knitr","ggmcmc","cowplot","gridExtra","ggcorrplot" )
req.pck_bool <- sapply(X=req.pck,FUN=require,character.only=T)
if(!all(req.pck_bool)) {
  sapply(X=req.pck[!req.pck_bool],FUN=install.packages,character.only=T);
  sapply(X=req.pck,FUN=require,character.only=T)
}

out_dir <- './inst'

# setwd(out_dir)

# Indicates if the MCMC will be computed (TRUE), or loaded from previously computed results (FALSE)
compute_mcmc = TRUE
n_iter_mcmc = 500e3
n_iter_mcmc_stage2 = 100

# Indicates if it is necessary to thin and warm-up the chain
n_warmup = 20e3
n_thin = 16

# Indicates if the summary of MCMC output will be computed (TRUE), or loaded from previously computed results (FALSE)
compute_summary = FALSE

plot_mcmc = FALSE

compute_model_eval = FALSE
plot_model_eval = FALSE

plot_ManureLevel_imp = FALSE

# Parallel processing
parallel_comp = TRUE
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
    # eta_PO_i=match(0.84,power_eta_all)
    
    cat("power_w_PO =", power_eta_all[eta_PO_i] ,"\n")
    
    ##### LOAD AGRICURB DATA #####
    # Variables that will be transformed in log-scale
    vars_log_i = c("Rainfall")
    
    # Load data #
    Agricurb_data <- agricurbayes::get_agricurb_data( arc_datasets="NMeso",
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
  summary_theta <- foreach::foreach( dataset_i = seq_along(arc_dataset), .errorhandling="pass", .combine="rbind"  ) %:%
    foreach( eta_PO_i = seq_along(power_eta_all), .errorhandling="pass", .combine="rbind" ) %dopar% {
      # dataset_i=1
      # eta_PO_i=match(0.840,power_eta_all)
      
      file_i = paste(out_dir, arc_dataset[dataset_i],"_smi_bayes_",formatC(power_eta_all[eta_PO_i],digits=3,format="f",flag="0"),".rds",sep="")
      # file.exists(file_i)
      
      mcmc_res = NULL
      if( file.exists(file_i) ) {
        mcmc_res <- readRDS( file=file_i )
      }
      
      # Summary of posteriors
      summary_theta_table <- data.frame(NULL)
      if( !is.null(mcmc_res) & !any(is.na(mcmc_res$theta_mcmc)) ){
        
        cat("dataset =",arc_dataset[dataset_i],", power_w_PO =", power_eta_all[eta_PO_i] ,"\n")
        
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
        summary_theta_table$arc_dataset = arc_dataset[dataset_i]
        
      }
      
      summary_theta_table
    }
  
  saveRDS( summary_theta, file=paste( out_dir, "agric_data_smi_summary.rds",sep="") )
  
}

# plotting MCMC results #
if(plot_mcmc) {
  
  # Loading summary of SMI posteriors
  summary_theta <- readRDS( file=paste( out_dir, "agric_data_smi_summary.rds",sep="") )
  
  # Setting nice ggplot theme settings
  agricurbayes::set_ggtheme()
  
  ### GAMMA credibility intervals ###
  p <- summary_theta %>%
    mutate( arc_dataset=factor(arc_dataset,levels=c("NMeso","Aegean","SWGermany")) ) %>%
    filter( param=="gamma_po_1" ) %>%
    tidyr::spread(statistic,value) %>%
    ggplot( aes(x=eta_PO,fill=arc_dataset) ) +
    geom_ribbon( aes(ymax=q97.5, ymin=q2.5), alpha=0.25 ) +
    geom_ribbon( aes(ymax=q75, ymin=q25), alpha=0.50 ) +
    geom_line( aes(y=Mean, col=arc_dataset) ) +
    geom_hline(yintercept=0, lty=2) +
    geom_vline( aes(xintercept=vl), lty=4, col="purple",
                data = data.frame( arc_dataset=factor(c("NMeso","Aegean","SWGermany"),levels=c("NMeso","Aegean","SWGermany")),
                                   vl=c(0.84,0.74,1.00) ) ) +
    facet_wrap( ~arc_dataset ) +
    coord_cartesian(xlim=c(0,1),ylim=5*c(-1,1)) +
    theme(legend.position="none") + # Remove legend
    labs( #title="Posterior mean and intervals of gamma",
      x="eta: PO module influence on ManureLevel", y="gamma" )
  # print(p)
  ggsave( plot=p,
          filename=paste(out_dir,"gamma_po_1_smi_interval_alldatasets.png",sep=""),
          device="png", width=20,height=8, units="cm")
  
  ### GAMMA Probability of negative values ###
  p <- summary_theta %>%
    mutate( arc_dataset=factor(arc_dataset,levels=c("NMeso","Aegean","SWGermany")) ) %>%
    filter( param=="gamma_po_1",
            statistic=="prob_leq_0") %>%
    select(one_of(c("arc_dataset","eta_PO","statistic","value"))) %>%
    tidyr::spread(statistic,value) %>%
    ggplot() +
    geom_line( aes(x=eta_PO, y=prob_leq_0,col=arc_dataset) ) +
    coord_cartesian(ylim=c(0,1),xlim=c(0,1)) +
    labs( #title="Posterior probability of a negative gamma",
      #subtitle="effect of Size on Manure Level",
      x="eta\n(degree of influence of PO module on imputation)", y="Pr( gamma<0 )" )
  # print(p)
  ggsave( plot=p,
          filename=paste(out_dir,"gamma_po_1_smi_prob_leq_0_alldatasets.png",sep=""),
          device="png", width=12,height=8, units="cm")
  
}

### Compute model evaluation ###
# Predict delta15N responses for archaeological data #
if( compute_model_eval ) {
  
  model_eval = foreach::foreach( dataset_i = seq_along(arc_dataset), .errorhandling="pass",.combine = rbind ) %:%
    foreach::foreach( eta_PO_i = seq_along(power_eta_all), .errorhandling="pass",.combine = rbind ) %dopar% {
      # dataset_i=3
      # eta_PO_i=match(0.6,power_eta_all)
      
      # list.files(out_dir)
      
      cat("dataset =",arc_dataset[dataset_i],", power_w_PO =", power_eta_all[eta_PO_i] ,"\n")
      
      file_i = paste(out_dir, arc_dataset[dataset_i],"_smi_bayes_",formatC(power_eta_all[eta_PO_i],digits=3,format="f",flag="0"),".rds",sep="")
      # file.exists(file_i)
      
      mcmc_res = NULL
      model_eval_i = data.frame(NULL)
      if( file.exists(file_i) ) {
        mcmc_res <- readRDS( file=file_i )
      }
      
      if( !is.null(mcmc_res) ){
        
        # plot(mcmc_res$theta_mcmc[,"gamma_po_1"])
        # plot(mcmc_res$theta_mcmc[,"alpha_po_1"])
        # plot(mcmc_res$theta_mcmc[,"alpha_po_2"])
        # plot( data.frame(mcmc_res$theta_mcmc[,"alpha_po_1"],mcmc_res$theta_mcmc[,"alpha_po_2"]),
        #       main=paste("dataset=",arc_dataset[dataset_i],", eta_PO=", power_eta_all[eta_PO_i] ,sep=""),
        #       xlab="alpha 1",ylab="alpha 2",
        #       pch=20,col="#0000FF30" );abline(a=0,b=1,col=2)
        # plot( data.frame(mcmc_res$theta_mcmc[,"gamma_po_1"],mcmc_res$theta_mcmc[,"alpha_po_2"]),
        #       main=paste("dataset=",arc_dataset[dataset_i],", eta_PO=", power_eta_all[eta_PO_i] ,sep=""),
        #       xlim=c(-5,5),xlab="gamma",
        #       ylim=c(-5,5),ylab="alpha 2",
        #       pch=20,col="#0000FF30" )
        # dim(mcmc_res$loglik_mcmc)
        
        
        # PO for archaelogical data
        # mcmc_res$loglik_mcmc[,1:mcmc_res$n_obs_arc,1]
        # traceplot(mcmc_res$theta_mcmc[,"gamma_po_1"])
        # plot(apply(mcmc_res$loglik_mcmc[,1:mcmc_res$n_obs_arc,1],1,sum),type="l") # trace of log likelihood in the model
        # plot(apply(mcmc_res$loglik_mcmc[,1:mcmc_res$n_obs_arc,1],2,mean))# trace of log likelihood for individual imputed rows
        
        # HM for archaelogical data
        # mcmc_res$loglik_mcmc[,1:mcmc_res$n_obs_arc,2]
        # plot(apply(mcmc_res$loglik_mcmc[,1:mcmc_res$n_obs_arc,2],1,sum),type="l")
        # HM for modern data
        # mcmc_res$loglik_mcmc[,-(1:mcmc_res$n_obs_arc),2]
        
        # elpd for predicting archaelogical data in HM model
        waic_i = loo::waic( mcmc_res$loglik_mcmc[,1:mcmc_res$n_obs_arc,2] )
        loo_i = loo::loo( mcmc_res$loglik_mcmc[,1:mcmc_res$n_obs_arc,2] )
        
        model_eval_i = data.frame(rbind(waic_i$estimates,loo_i$estimates))
        model_eval_i$stat = rownames(model_eval_i)
        rownames(model_eval_i) = NULL
        
        model_eval_i$arc_dataset = arc_dataset[dataset_i]
        model_eval_i$eta_PO = power_eta_all[eta_PO_i]
        if( (length(dim(model_eval_i))!=2) | ncol(model_eval_i)!=5 ) {
          cat("Problem with dataset =",arc_dataset[dataset_i],", power_w_PO =", power_eta_all[eta_PO_i] ,"\n")
          cat(dim(model_eval_i))
          model_eval_i = NULL
        }
      }
      
      model_eval_i
    }
  
  model_eval$arc_dataset = factor(model_eval$arc_dataset, levels=arc_dataset)
  saveRDS( model_eval, file=paste( out_dir, "agric_data_smi_model_eval.rds",sep="") )
  
}

# Select best eta
if(T){
  # Loading elpd estimates for SMI posteriors
  model_eval <- readRDS( file=paste( out_dir, "agric_data_smi_model_eval.rds",sep="") )
  
  elpd_approx = c('elpd_waic','elpd_loo')[1]
  
  elpd_gp_estimate = foreach::foreach( dataset_i = seq_along(arc_dataset), .errorhandling="pass",.combine = rbind ) %dopar% {
    # dataset_i = 3
    
    aux = model_eval %>%
      dplyr::filter(stat==elpd_approx) %>%
      dplyr::filter(arc_dataset==c("NMeso","Aegean","SWGermany")[dataset_i])
    gp_elpd = GPfit::GP_fit(X=aux$eta_PO,Y=aux$Estimate)
    # eta_PO_optim = seq(0,1,0.02)[which.max(predict(gp_elpd, xnew=seq(0,1,0.01))$Y_hat)]
    # plot(gp_elpd); abline(v=eta_PO_optim,col="purple",lty=3); mtext(text=eta_PO_optim,side=3,at=eta_PO_optim)
    elpd_hat_i = data.frame( eta_PO=seq(0,1,0.02),
                             arc_dataset=arc_dataset[dataset_i],
                             elpd_hat=predict(gp_elpd, xnew=seq(0,1,0.02))$Y_hat )
    elpd_hat_i
  }
  elpd_gp_estimate$arc_dataset = factor(elpd_gp_estimate$arc_dataset,levels=arc_dataset)
  
  elpd_gp_best = elpd_gp_estimate %>%
    group_by(arc_dataset) %>%
    filter(elpd_hat == max(elpd_hat)) %>%
    as.data.frame()
}
# elpd_gp_best

# plotting model evaluation results #
if( plot_model_eval ){
  agricurbayes::set_ggtheme()
  p <- model_eval %>%
    dplyr::filter(stat==elpd_approx) %>%
    # dplyr::filter(arc_dataset=='Aegean') %>%
    # dplyr::filter(eta_PO<1) %>%
    ggplot() +
    # geom_line( aes(x=eta_PO,y=Estimate,col=arc_dataset) ) +
    geom_point( aes(x=eta_PO,y=Estimate, col=arc_dataset) ) +
    geom_line( aes(x=eta_PO,y=elpd_hat,col=arc_dataset), lty=3, data=elpd_gp_estimate ) +
    facet_wrap( ~arc_dataset, scales = "free" )+
    geom_vline( aes(xintercept=eta_PO), col="purple", lty=2, data=elpd_gp_best ) +
    geom_text(aes(x=eta_PO,y=elpd_hat,label=eta_PO), data=elpd_gp_best)+
    # annotate("text", x = 2:5, y = 25, label = paste("eta="))+
    theme(legend.position="none") + # Remove legend
    labs(y="elpd estimate")
  # print(p)
  ggsave( plot=p,
          filename=paste(out_dir,"agric_data_smi_model_eval_",elpd_approx,".png",sep=""),
          device="png", width=25,height=10, units="cm")
  
  
  # Gather the best posteriors for each dataset #
  p_mcmc_trace = p_mcmc_autocor = p_mcmc_crosscor = list(NULL)
  mcmc_best_smi <- foreach::foreach( dataset_i = seq(arc_dataset) ) %do% {
    # dataset_i=1
    
    # File with MCMC results #
    eta_star_i = elpd_gp_best[elpd_gp_best$arc_dataset==arc_dataset[dataset_i],"eta_PO"]
    file_i = paste(out_dir, arc_dataset[dataset_i],"_smi_bayes_",formatC(eta_star_i,digits=3,format="f",flag="0"),".rds",sep="")
    
    mcmc_res_i <- readRDS( file=file_i )
    
    # Plot convergence analysis for parameters of interest
    all_par = colnames(mcmc_res_i$theta_mcmc)
    interest_par = c("alpha_po_1","alpha_po_2","gamma_po_1",all_par[substr(all_par,1,7)=='beta_hm'],"sigma_hm")
    interest_par_tempered = c("alpha_po_1","alpha_po_2","gamma_po_1")
    
    mcmc_1 = mcmc_res_i$theta_mcmc[,interest_par]
    mcmc_2 = mcmc_res_i$theta_mcmc_tempered[,interest_par_tempered]; colnames(mcmc_2)=paste(colnames(mcmc_2),'_tilde',sep="")
    # Thinning and warm-up #
    mcmc_interest <- as.mcmc( cbind(mcmc_1,mcmc_2) )
    
    
    # effectiveSize(mcmc_interest)
    
    model_S <- ggmcmc::ggs( mcmc_interest )
    # Trace plot
    p_mcmc_trace[[dataset_i]] = ggs_traceplot(model_S) + labs( title=paste(arc_dataset[dataset_i]," eta=",formatC(eta_star_i,digits=3,format="f",flag="0"),sep="") )
    # Autocorrelation
    p_mcmc_autocor[[dataset_i]] = ggs_autocorrelation(model_S) + labs( title=paste(arc_dataset[dataset_i]," eta=",formatC(eta_star_i,digits=3,format="f",flag="0"),sep="") )
    # Cross correlation
    cor_mat = round(cor( mcmc_interest ),2); diag(cor_mat)=NA; 
    p_mcmc_crosscor[[dataset_i]] = cor_mat[colnames(mcmc_interest),rev(colnames(mcmc_interest))] %>%
      ggcorrplot::ggcorrplot( method = "circle", type='full', lab = TRUE, lab_size=2) +
      labs( title=paste(arc_dataset[dataset_i]," eta=",formatC(eta_star_i,digits=3,format="f",flag="0"),sep="") )
    rm(all_par,interest_par,model_S,cor_mat)
    
    # names(mcmc_res_i)
    mcmc_res_i
  }
  
  # Plot MCMC traces for the main parameters in the model
  p_mcmc_trace_all = cowplot::plot_grid(
    p_mcmc_trace[[1]],p_mcmc_trace[[2]],p_mcmc_trace[[3]],
    labels = NULL,
    hjust = -1,
    nrow = 1
  )
  ggsave( plot=p_mcmc_trace_all,
          filename=paste(out_dir,"agric_model_mcmc_trace.png",sep=""),
          device="png", width=30,height=30, units="cm")
  
  # Plot MCMC autocorrelation for the main parameters in the model
  p_mcmc_autocor_all = cowplot::plot_grid(
    p_mcmc_autocor[[1]],p_mcmc_autocor[[2]],p_mcmc_autocor[[3]],
    labels = NULL,
    hjust = -1,
    nrow = 1
  )
  ggsave( plot=p_mcmc_autocor_all,
          filename=paste(out_dir,"agric_model_mcmc_autocorrelation.png",sep=""),
          device="png", width=30,height=20, units="cm")
  
  # Plot MCMC Cross correlation for the main parameters in the model
  p_mcmc_crosscor_all = cowplot::plot_grid(
    p_mcmc_crosscor[[1]],p_mcmc_crosscor[[2]],p_mcmc_crosscor[[3]],
    labels = NULL,
    hjust = -1,
    nrow = 1
  )
  ggsave( plot=p_mcmc_crosscor_all,
          filename=paste(out_dir,"agric_model_mcmc_crosscorrelation.png",sep=""),
          device="png", width=40,height=20, units="cm")
  
  # Comparison of posterior distribution of gamma
  mcmc_best_smi_gamma = foreach::foreach( dataset_i = seq(arc_dataset), .combine = cbind ) %do% {
    mcmc_best_smi[[dataset_i]][["theta_mcmc"]][,"gamma_po_1"]
  }
  colnames(mcmc_best_smi_gamma) = arc_dataset
  
  agricurbayes::set_ggtheme()
  p <- mcmc_best_smi_gamma %>%
    as.data.frame() %>%
    reshape2::melt( measure.vars=arc_dataset,
                    variable.name="dataset" ) %>%
    # dplyr::filter(dataset!="NMeso") %>%
    ggplot() +
    geom_density( aes(x=value, fill=dataset), alpha=0.25 ) +
    geom_vline(xintercept=0,lty=2) +
    labs(x="gamma")
  # print(p)
  ggsave( filename=paste(out_dir,"gamma_po_1_smi_alldatasets.png",sep=""),
          plot=p, device="png", width=15,height=10, units="cm", dpi="print")
  
  # Posterior probability of negative gamma #
  p_gamma_leq_0_smi <- apply( mcmc_best_smi_gamma < 0, 2, mean )
  p_gamma_leq_0_smi
  
  # Bayes factor for negative gamma #
  BF_gamma_leq_0_smi = p_gamma_leq_0_smi / (1-p_gamma_leq_0_smi)
  BF_gamma_leq_0_smi
  
}


### Plot best bayesian imputation ###
if( plot_ManureLevel_imp ){
  foreach::foreach( dataset_i = seq_along(arc_dataset), .errorhandling="pass" ) %do% {
    # dataset_i=1
    arc_dataset_i = arc_dataset[dataset_i]
    
    ##### LOAD AGRICURB DATA #####
    # Variables that will be transformed in log-scale
    vars_log_i = c("Rainfall")
    if( is.element(arc_dataset[dataset_i],c("Aegean","SWGermany")) ) {
      vars_log_i = c( vars_log_i, "Size" )
    }
    
    # Load data #
    Agricurb_data <- agricurbayes::get_agricurb_data( arc_datasets=arc_dataset[dataset_i],
                                                      vars_log=vars_log_i )
    
    # Separates data into modern and archaeological
    data_arc = Agricurb_data %>%
      dplyr::filter(dataset!="modern") %>%
      dplyr::select(c( "Size",
                       "Site",
                       "Rainfall_min",
                       "Rainfall_max",
                       "Category",
                       "normd15N"))
    
    # File with MCMC results #
    eta_star_i = round(elpd_gp_best[elpd_gp_best$arc_dataset==arc_dataset_i,"eta_PO"],6)
    file_i = paste(out_dir, arc_dataset[dataset_i],"_smi_bayes_",formatC(eta_star_i,digits=3,format="f",flag="0"),".rds",sep="")
    
    mcmc_res <- readRDS( file=file_i )
    
    # Calculates probabilities of imputed values for all levels in ManureLevel
    ManureLevel_post_prob = foreach::foreach( ManureLevel_i=1:3, .combine=rbind) %dopar% {
      # ManureLevel_i=1
      apply(mcmc_res$ManureLevel_imp_mcmc==ManureLevel_i,2,mean)
    }
    rownames(ManureLevel_post_prob)=ManureLevels
    ManureLevel_post_prob = ManureLevel_post_prob %>%
      t() %>% as.data.frame() %>%
      dplyr::mutate( Size=data_arc$Size,
                     Site=data_arc$Site,
                     prob_leq_low=low,
                     prob_leq_medium=low+medium )
    
    # Calculates posterior prob of ManureLevels<=ManureLevels_i for different sizes
    n_out=100
    Size_seq = seq(min(ManureLevel_post_prob$Size),max(ManureLevel_post_prob$Size),length.out=n_out)
    Size_seq_norm = (Size_seq - mean(Agricurb_data$Size,na.rm=T))/sd(Agricurb_data$Size,na.rm=T)
    comb <- function(...) {
      mapply('rbind', ..., SIMPLIFY=FALSE)
    }
    PO_post_probs = foreach::foreach( iter_i = 1:nrow(mcmc_res$theta_mcmc), .combine='comb', .multicombine=TRUE ) %dopar% {
      # iter_i=1
      
      # parameters in PO module for iteration iter_i
      PO_param_i = mcmc_res[["theta_mcmc"]][iter_i,c("alpha_po_1","alpha_po_2","gamma_po_1")]
      
      # Probabilities
      prob_iter_i = data.frame(
        Size=Size_seq,
        low=loglik_PO_i_cpp( Y=matrix( 1, n_out, 1), X=matrix( Size_seq_norm, n_out, 1), alpha=PO_param_i[1:2], beta=PO_param_i[3] ),
        medium=loglik_PO_i_cpp( Y=matrix( 2, n_out, 1), X=matrix( Size_seq_norm, n_out, 1), alpha=PO_param_i[1:2], beta=PO_param_i[3] ),
        high=loglik_PO_i_cpp( Y=matrix( 3, n_out, 1), X=matrix( Size_seq_norm, n_out, 1), alpha=PO_param_i[1:2], beta=PO_param_i[3] )
      ) %>% exp() %>%
        dplyr::mutate( prob_leq_low=low,
                       prob_leq_medium=low+medium )
      
      # output
      list(prob_iter_i$prob_leq_low , prob_iter_i$prob_leq_medium)
    }
    names(PO_post_probs) = c("p_leq_low","p_leq_medium")
    # hist(PO_post_probs[[2]][,100],xlim=c(0,1))
    
    probs_aux=c(2.5,25,50,75,97.5)/100
    PO_post_probs_summary = foreach::foreach( i=1:2 ) %dopar% {
      summary_i = data.frame( Size=Size_seq,
                              Mean=apply(PO_post_probs[[i]],2,mean),
                              apply(PO_post_probs[[i]],2,quantile,probs=probs_aux)%>%t() )
      colnames(summary_i)[-(1:2)] = paste("q",probs_aux*100,sep="")
      summary_i
    }
    names(PO_post_probs_summary) = c("p_leq_low","p_leq_medium")
    
    # Setting nice ggplot theme settings
    agricurbayes::set_ggtheme()
    
    p_leq_low = ManureLevel_post_prob %>%
      ggplot() +
      geom_ribbon( aes(x=Size, ymax=q97.5, ymin=q2.5), alpha=0.125, data=PO_post_probs_summary[["p_leq_low"]] ) +
      geom_ribbon( aes(x=Size, ymax=q75, ymin=q25), alpha=0.250, data=PO_post_probs_summary[["p_leq_low"]] ) +
      geom_line( aes(x=Size, y=Mean), data=PO_post_probs_summary[["p_leq_low"]] ) +
      geom_jitter(aes(x=Size, y=prob_leq_low, col=Site, pch=Site), alpha=0.5) +
      scale_shape_manual(values=seq(0,15)) +
      theme(legend.position="right") + # Right legend
      # theme(legend.position="none") + # Remove legend
      coord_cartesian(ylim=c(0,1)) +
      labs( x=set_names(c("Size","log(1+Size)","log(1+Size)"),arc_dataset)[arc_dataset_i],
            y="Pr( ManureLevel <= low )")
    p_leq_medium = ManureLevel_post_prob %>%
      ggplot() +
      geom_ribbon( aes(x=Size, ymax=q97.5, ymin=q2.5), alpha=0.125, data=PO_post_probs_summary[["p_leq_medium"]] ) +
      geom_ribbon( aes(x=Size, ymax=q75, ymin=q25), alpha=0.250, data=PO_post_probs_summary[["p_leq_medium"]] ) +
      geom_line( aes(x=Size, y=Mean), data=PO_post_probs_summary[["p_leq_medium"]] ) +
      geom_jitter(aes(x=Size, y=prob_leq_medium, col=Site, pch=Site), alpha=0.5) +
      scale_shape_manual(values=seq(0,15)) +
      theme(legend.position="right") + # Right legend
      # theme(legend.position="none") + # Remove legend
      coord_cartesian(ylim=c(0,1)) +
      labs( x=set_names(c("Size","log(1+Size)","log(1+Size)"),arc_dataset)[arc_dataset_i],
            y="Pr( ManureLevel <= medium )")
    
    prow = cowplot::plot_grid(
      p_leq_low + theme(legend.position="none"),
      p_leq_medium + theme(legend.position="none"),
      # align = 'vh',
      labels = c('A', 'B'),
      hjust = -1,
      nrow = 1
    )
    
    # extract the legend from one of the plots
    legend <- get_legend(
      # create some space to the left of the legend
      p_leq_low + theme(legend.box.margin = margin(0, 0, 0, 12))
    )
    
    # add the legend to the row we made earlier. Give it one-third of 
    # the width of one plot (via rel_widths).
    p <- plot_grid(prow, legend, rel_widths = c(10, 1.75))
    
    # now add the title
    title <- ggdraw() +  draw_label( arc_dataset_i,
                                     fontface = 'bold', x = 0, hjust = 0 ) +
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      theme( plot.margin = margin(0, 0, 0, 7) )
    p_title = plot_grid( title, p, ncol = 1,
                         # rel_heights values control vertical title margins
                         rel_heights = c(0.1, 1) )
    
    ggsave( plot=p_title,
            filename=paste(out_dir,"arc_ManureLevel_vs_Size_smi_",arc_dataset_i,".png",sep=""),
            device="png", width=20,height=8, units="cm")
    
    NULL
  }
}
