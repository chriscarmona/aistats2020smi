rm(list = ls())
options(scipen=999, stringsAsFactors=FALSE)
set.seed(0)

# pkgbuild::compile_dll()
# Rcpp::compileAttributes()
# roxygen2::roxygenise()
# devtools::document()
# devtools::install()

# loading required packages #
req.pck <- c( "aistats2020smi",
              "MASS","nlme","coda","ordinal","LearnBayes",
              "plyr","tidyverse", "ggplot2","GGally","lattice",
              "foreach","doParallel","doRNG",
              "GPfit","rBayesianOptimization",
              "RColorBrewer","knitr","ggmcmc","cowplot","gridExtra","ggcorrplot" )
req.pck_bool <- sapply(X=req.pck,FUN=require,character.only=T)
if(!all(req.pck_bool)) {
  sapply(X=req.pck[!req.pck_bool],FUN=install.packages,character.only=T);
  sapply(X=req.pck,FUN=require,character.only=T)
}

out_dir <- './inst'

# Parallel processing
parallel_comp = TRUE
if(parallel_comp){
  n_cores = 6
  options(cores=n_cores)
  doParallel::registerDoParallel()
  getDoParWorkers()
}

# Data from Plummer (2014)
aistats2020smi::HPV %>% head()
# HPV <- data.frame(
#   ncases = c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194), # Y
#   Npop = c(26983, 250930, 829348, 157775, 150467, 352445, 553066, 26751, 75815, 150302, 354993, 3683043, 507218), # T
#   nhpv = c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4), # Z
#   Npart = c(111, 71, 162, 188, 145, 215, 166, 37, 173, 143, 229, 696, 93) ) # N
# save(HPV,file="./data/HPV.rda")

# Levels of influence in SMI
power_eta_all <- c( 0.01,0.05,0.99,
                    seq(0,1.00,by=0.10) )
power_eta_all <- sort( unique(round(power_eta_all,6)) )

compute_mcmc=FALSE
if(compute_mcmc) {
  foreach( eta_i = seq_along(power_eta_all) ) %do% {
    # eta_i <- 1
    eta_pois <- power_eta_all[eta_i]
    cat(eta_pois,", ")
    file_i <- paste(out_dir,"/HPV_partial_cut_stan_",formatC(eta_pois,digits=3,format="f",flag="0"),".rds",sep="")
    if( !file.exists(file_i) ){
      set.seed(0)
      hpv_smi_mcmc_i <- aistats2020smi::mcmc_hpv( HPV=aistats2020smi::HPV,

                                                  # Number of iterations
                                                  n_iter_mcmc = 5000, # main chain
                                                  n_iter_warmup = 1000,
                                                  n_chains_mcmc = 4, # Number of chains
                                                  n_iter_mcmc_stage_2 = 500, # Subchain

                                                  # Cut rate
                                                  eta_pois = eta_pois,
                                                  eta_binom = 1,

                                                  mcmc_file = file_i,
                                                  n_cores=n_cores )
    }
    NULL
  }
}

compute_summary=FALSE
if(compute_summary){
  hpv_smi_mcmc_all = foreach( eta_i = seq_along(power_eta_all), .combine="rbind" ) %do% {
    # eta_i <- 1
    eta_pois <- power_eta_all[eta_i]
    cat(eta_pois,", ")
    file_i <- paste(out_dir,"/HPV_partial_cut_stan_",formatC(eta_pois,digits=3,format="f",flag="0"),".rds",sep="")
    hpv_smi_mcmc_i = NULL
    if( file.exists(file_i) ){
      hpv_smi_mcmc_i <- readRDS(file = file_i)
      hpv_smi_mcmc_i$eta = eta_pois
    }
    hpv_smi_mcmc_i
  }
  save(hpv_smi_mcmc_all, file="./data/hpv_smi_mcmc_all.rda")
}

plot_mcmc=TRUE
if(plot_mcmc) {
  # load(file="./data/hpv_smi_mcmc_all.rda")

  # Yellow and black
  col_aux = colorRampPalette(c("#000000","#FFD500"))( 2 )

  # Joint posterior of theta1 and theta2 under the full and cut model
  p_theta_joint_cut_full = aistats2020smi::hpv_smi_mcmc_all %>%
    dplyr::select(one_of(c("eta",'theta1','theta2'))) %>%
    filter( eta %in% c(0,1) ) %>%
    mutate( eta=as.factor(eta)) %>%
    ggplot()+
    geom_point( aes(x=theta1,y=theta2,col=eta), alpha=0.1 ) +
    coord_cartesian(xlim=c(-3,-1),ylim=c(0,40)) +
    scale_color_manual(values=col_aux, name="eta")+
    guides(colour=guide_legend(override.aes = list(alpha=1)))+
    # labs(title="Posterior distribution of theta",subtitle="Partial Cut method") +
    theme_bw() + theme(legend.position = "bottom")
  ggsave( plot=p_theta_joint_cut_full,
          filename="hpv_smi_theta_joint_post.pdf",
          device="pdf",width=10,height=10,units="cm")

  # Marginal posterior of theta1 and theta2 under SMI for all eta
  aistats2020smi::set_ggtheme()
  p_theta_post = aistats2020smi::hpv_smi_mcmc_all %>%
    dplyr::select(one_of(c("eta",'theta1','theta2'))) %>%
    tidyr::pivot_longer(c('theta1','theta2'), names_to = "param") %>%
    ggplot( aes(x=value,y=eta, group=eta) ) +
    ggridges::geom_density_ridges(scale = 3, alpha=0.5 ) +
    scale_y_continuous(expand = c(0, 0)) +     # will generally have to set the `expand` option
    scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
    facet_wrap( ~param, scales = "free_x" ) +
    ggridges::theme_ridges()
  ggsave( plot=p_theta_post,
          filename="hpv_smi_theta_post.pdf",
          device="pdf", width=12,height=12, units="cm")
}


compute_model_eval = FALSE
if(compute_model_eval) {
  n_obs_hpv = nrow(aistats2020smi::HPV)
  hpv_smi_model_eval <- foreach( eta_i = seq_along(power_eta_all), .combine=rbind ) %dorng% {
    # eta_i <- 1
    loglik_aux <- apply( aistats2020smi::hpv_smi_mcmc_all %>% filter(eta==power_eta_all[eta_i]), 1,
                         FUN = function(x,HPV) {
                           aistats2020smi::hpv_loglik( data=HPV,
                                                       theta=x[paste("theta",1:2,sep="")],
                                                       phi=x[paste("phi_",1:n_obs_hpv,sep="")] ) },
                         HPV = aistats2020smi::HPV )
    loglik_aux <- t(loglik_aux); rownames(loglik_aux)<-NULL
    # waic_aux <- data.frame( poisson=MissBayes::waic( loglik_aux[,1:n_obs_hpv] ),
    #                         binomial=MissBayes::waic( loglik_aux[,n_obs_hpv+1:n_obs_hpv] ) )
    ll_pois <- loglik_aux[,1:n_obs_hpv]
    ll_binom <- loglik_aux[,n_obs_hpv+1:n_obs_hpv]
    waic_aux <- data.frame( poisson = c( loo::waic( ll_pois )$estimates[,1],
                                         loo::loo( ll_pois, r_eff=loo::relative_eff(exp(ll_pois),chain_id=rep(1,nrow(ll_pois))) )$estimates[,1] ),
                            binomial = c( loo::waic( ll_binom )$estimates[,1],
                                          loo::loo( ll_binom, r_eff=loo::relative_eff(exp(ll_binom),chain_id=rep(1,nrow(ll_binom))) )$estimates[,1] ) )

    waic_aux <- waic_aux %>%
      dplyr::mutate( score_id = gsub("loo","psis",rownames(waic_aux))) %>%
      tidyr::gather( key=module,value=score,-score_id) %>%
      dplyr::mutate( eta_pois = power_eta_all[eta_i] )
    waic_aux
  }
  save( hpv_smi_model_eval, file="./data/hpv_smi_model_eval.rda")
}

plot_model_eval = TRUE
if(plot_model_eval) {
  aistats2020smi::set_ggtheme()
  aistats2020smi::hpv_smi_model_eval %>%
    dplyr::mutate( module = dplyr::recode(module, 'poisson'='-elpd poisson', 'binomial'='-elpd binomial') ) %>%
    dplyr::filter(score_id=='elpd_waic') %>%
    ggplot( aes(x=eta_pois, y=-score) ) +
    geom_line( col="red" ) +
    geom_point( col="red" ) +
    facet_wrap( vars(module), ncol=2, scales = "free_y" ) +
    theme( axis.title.y=element_blank() )
}
