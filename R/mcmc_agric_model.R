#' @title
#'     MCMC procedure for hierarchical model of agricultural practices in early civilizations
#'
#' @description
#'     This function perform MCMC sampling for the hierarchical model of agricultural data.
#'     There are two modules in this model: Gaussian response (HM), and Proportional Odds (PO) prior for missing ManureLevel.
#'
#' @param data_arc Data frame. Archaeological data
#' @param data_mod Data frame. Modern data
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
#' @useDynLib aistats2020smi
#'
#' @export

mcmc_agric_model <- function( data_arc,
                              data_mod,
                              
                              impute_spec=c("full","cut","smi")[1],
                              
                              power_w_PO=1, power_w_HM=1,
                              
                              prior_spec_PO=c("flat","proper")[2],
                              prior_spec_HM=c("flat","proper")[2],
                              
                              PO_site_rnd_eff=TRUE,
                              HM_site_rnd_eff=TRUE,
                              
                              n_iter=100000,
                              n_iter_sub=10,
                              
                              theta_ini=NULL,
                              theta_min_max=NULL,
                              theta_prop_int=NULL,
                              theta_prop_kernel=c("norm","unif")[1],
                              
                              n_epoch_adapt=0,
                              n_iter_adapt=2000,
                              
                              ManureLevel_arc_imp_ini=NULL,
                              Rainfall_arc_imp_ini=NULL,
                              
                              imp_playpen=FALSE,
                              
                              gibbs_hm=TRUE,
                              
                              PO_expand=FALSE,
                              POmixda=TRUE,
                              lambda_prop_int=0.2,
                              lambda_mix_PO_prior=c(1,1),
                              
                              keep_imp=FALSE,
                              keep_ll=FALSE,
                              elpd=FALSE,
                              
                              out_file_rds=NULL,
                              log_file=NULL,
                              devel=FALSE ) {
  
  # on.exit(browser())
  
  # PO1:
  #    ManureLevel ∼ 1 + Size + (1|Site)
  #    link = "logistic"
  # HM1:
  #    normd15N ∼ 1 + Rainfall + ManureLevel + (1|Site)
  #    weights = varIdent(form=~1|Category)
  
  # This makes model.matrix avoid lossing rows when NA values are present
  prev_na.action <- options('na.action')
  options(na.action='na.pass')
  on.exit(options(na.action=prev_na.action$na.action))
  
  if( !is.null(out_file_rds) ) {
    # Checking that its possible to save the results #
    mcmc_res <- NULL
    saveRDS( mcmc_res, file=out_file_rds )
    file.remove(out_file_rds)
  }
  
  mcmc_clock <- Sys.time()
  if( !is.null(log_file) ) {
    cat("**** MCMC for agricultural model *****\n\n",
        
        "----- Model -----\n",
        "impute_spec = ",impute_spec,"\n",
        
        "power_w_PO = ",power_w_PO,"\n",
        "power_w_HM = ",power_w_HM,"\n",
        
        "prior_spec_PO = ",prior_spec_PO,"\n",
        "prior_spec_HM = ",prior_spec_HM,"\n",
        
        "PO_site_rnd_eff = ",PO_site_rnd_eff,"\n",
        "HM_site_rnd_eff = ",HM_site_rnd_eff,"\n",
        
        "\n\n",
        
        "----- MCMC -----\n\n",
        "n_iter = ",n_iter,"\n",
        
        "n_epoch_adapt = ", n_epoch_adapt,"\n",
        "n_iter_adapt = ", n_iter_adapt,"\n",
        
        "keep_imp = ",keep_imp,"\n",
        
        "\n\n",
        
        "---------------------------\n\n",
        "Starting time:\n",as.character(mcmc_clock),"\n\n",
        "---------------------------\n\n",
        file=log_file)
  }
  
  ManureLevels <- c("low","medium","high")
  k.M <- length(ManureLevels) # number of categories in "ManureLevel"
  
  Sites_arc <- sort(unique(data_arc$Site))
  n_sites_arc <- length(Sites_arc) # number of categories in "Sites"
  Sites_mod <- sort(unique(data_mod$Site))
  n_sites_mod <- length(Sites_mod) # number of categories in "Sites"
  
  # indicator of variance offset
  if( any(!is.element(data_arc$Category, c("wheat","barley"))) ) {stop("There are more cereal categories than expected.")}
  data_arc$ind_v <- ifelse(data_arc$Category=="wheat",1,0)
  if( any(!is.element(data_mod$Category, c("wheat","barley"))) ) {stop("There are more cereal categories than expected.")}
  data_mod$ind_v <- ifelse(data_mod$Category=="wheat",1,0)
  
  # Identification of parameters in the model #
  theta_names <- c( paste("alpha_po_",1:(k.M-1),sep=""), # "alpha_po" intercepts of the PO model
                    "gamma_po_1", # "gamma_po_1" fixed effect of PO model. Size on ManureLevel
                    paste("eta_po_",Sites_arc,sep=""), # random effects of Site on ManureLevel
                    "sigma_eta_PO", # "sigma_eta_PO" variance of random effects of Site on ManureLevel
                    paste("beta_hm_",1:(1+k.M),sep=""), # "beta_hm" fixed effects of the HM model.
                    paste("eta_hm_",c(Sites_arc,Sites_mod),sep=""), # random effects for Site in the HM model
                    "sigma_hm", # HM model variance
                    "v", # HM variance offset for ind_v rows
                    "sigma_hm_eta" ) # variance of random effects for Site in the HM model
  
  # Identification of parameters only in the PO module #
  theta_names_po <- c( paste("alpha_po_",1:(k.M-1),sep=""), # "alpha_po" intercepts of the PO model
                       "gamma_po_1", # "gamma_po_1" fixed effect of PO model. Size on ManureLevel
                       paste("eta_po_",Sites_arc,sep=""), # random effects of Site on ManureLevel
                       "sigma_eta_PO") # "sigma_eta_PO" variance of random effects of Site on ManureLevel
  
  # number of parameters #
  n_par <- length(theta_names)
  n_par_po <- length(theta_names_po)
  
  ### Min and max possible values of theta ###
  aux <- theta_min_max
  
  theta_min_max <- matrix(c(-Inf,Inf), nrow=n_par, ncol=2, byrow=T )
  rownames(theta_min_max) <- theta_names
  colnames(theta_min_max) <- c("theta.min","theta.max")
  # PO model #
  theta_min_max["sigma_eta_PO",1] <- 0 # positive variance of random effects
  # HM model #
  theta_min_max["sigma_hm",1] <- 0 # positive variance
  theta_min_max["v",1] <- 0 # positive variance offset
  theta_min_max["sigma_hm_eta",1] <- 0 # positive variance of random effects
  
  if( !is.null(aux) ) {
    # theta_min_max
    theta_min_max[rownames(aux),] <- aux[]
    # theta_min_max
  }
  rm(aux)
  
  # Length of the interval for the uniform proposal of each theta
  aux <- theta_prop_int
  
  theta_prop_int <- setNames(rep(1,n_par),theta_names)
  
  # the uniform interval for the "variance" terms goes
  # sigma * unif( 1/theta_prop_int , theta_prop_int )
  theta_prop_int[c("sigma_eta_PO","sigma_hm","v","sigma_hm_eta")] <- 1
  
  if( !is.null(aux) ) {
    theta_prop_int[names(aux)] <- aux[]
  }
  if( !all(is.element(theta_names,names(theta_prop_int))) ) { stop('The proposal intervals in "theta_prop_int" are not defined for all the parameters') }
  theta_prop_int <- theta_prop_int[theta_names]
  rm(aux)
  
  ### Initial values ###
  theta_mcmc_ini <- setNames(rep(NA,n_par),theta_names)
  
  # HM parameters #
  data_mod$ManureLevel <- factor( match(data_mod$ManureLevel,c("low","medium","high")), levels=1:3)
  data_mod$Site <- as.factor( data_mod$Site )
  data_mod$ind_v <- as.factor( data_mod$ind_v )
  model_hm <- nlme::lme( normd15N ~ Rainfall + ManureLevel,
                         random=~1|Site,
                         weights=varIdent(form=~1|ind_v),
                         data=data_mod, method='REML' )
  beta_hm = fixed.effects(model_hm)
  eta_hm = random.effects(model_hm)
  
  theta_mcmc_ini[paste("beta_hm_",1:4,sep="")] <- fixed.effects(model_hm)
  theta_mcmc_ini["sigma_hm_eta"] <- intervals(model_hm)$reStruct$Site[[2]]
  theta_mcmc_ini["sigma_hm"] <- intervals(model_hm)$sigma[2]
  theta_mcmc_ini["v"]<-(intervals(model_hm)$varStruct[2])^2
  theta_mcmc_ini[paste("eta_hm_",rownames(eta_hm),sep="")] <- eta_hm[,1]
  theta_mcmc_ini[paste("eta_hm_",unique(data_arc$Site),sep="")] <- 0
  
  ### Imputing Missing data ##
  data_arc_imp = data_arc
  
  # Rainfall #
  Rainfall_arc_imp_ini = (data_arc$Rainfall_min + data_arc$Rainfall_max )/2
  data_arc_imp$Rainfall = Rainfall_arc_imp_ini
  
  # ManureLevel #
  data_arc_imp$ManureLevel = factor(NA,levels=1:3)
  yhat = matrix(NA,nrow(data_arc),3)
  for(i in 1:3) {
    data_arc_imp$ManureLevel[] = i
    Xm <- model.matrix( normd15N ~ Rainfall + as.factor(ManureLevel), data=data_arc_imp )
    yhat[,i] <- Xm%*%beta_hm
  }
  ManureLevel_arc_imp_ini <- apply((yhat-data_arc_imp$normd15N)^2,1,which.min)
  data_arc_imp$ManureLevel = factor( ManureLevel_arc_imp_ini, levels=1:3 )
  
  # PO parameters #
  data_arc_imp$Site <- as.factor(data_arc_imp$Site)
  if(T) {
    model_po <- ordinal::clmm( ManureLevel ~ 1 + Size + (1|Site) , data=data_arc_imp )
    theta_mcmc_ini[c("alpha_po_1","alpha_po_2")] <- model_po$alpha
    theta_mcmc_ini["gamma_po_1"] <- model_po$beta
    eta_po <- random.effects(model_po)$Site
    theta_mcmc_ini["sigma_eta_PO"] <- unlist(model_po$ST,use.names=FALSE)
    theta_mcmc_ini[paste("eta_po_",rownames(eta_po),sep="")] <- eta_po[,1]
    
    rm(model_po)
  } else {
    model_po <- ordinal::clm( ManureLevel ~ 1 + Size, data=data_arc_imp )
    theta_mcmc_ini[c("alpha_po_1","alpha_po_2")] <- model_po$alpha
    theta_mcmc_ini["gamma_po_1"] <- model_po$beta
    theta_mcmc_ini["sigma_eta_PO"] <- 1
    theta_mcmc_ini[substr(names(theta_mcmc_ini),1,6)=="eta_po"] <- 0
    rm(model_po)
  }
  
  theta_singleimpute <- theta_mcmc_ini
  
  # User defined initial values #
  if(!is.null(theta_ini)) {
    theta_mcmc_ini[ names(theta_ini) ] <- theta_ini[]
  }
  
  
  # Keep only columns with relevant information for the C++ routine #
  data_arc <- data_arc[,c( "Size",
                           "Site",
                           "Rainfall_min",
                           "Rainfall_max",
                           "ind_v",
                           "normd15N")]
  data_mod <- data_mod[,c( "Site",
                           "ManureLevel",
                           "Rainfall",
                           "ind_v",
                           "normd15N")]
  
  # transform Site to integer #
  data_arc$Site <- match(data_arc$Site,c(Sites_arc,Sites_mod))
  data_mod$Site <- match(as.character(data_mod$Site),c(Sites_arc,Sites_mod))
  
  # transform ManureLevel to integer #
  data_mod$ManureLevel <- as.numeric(data_mod$ManureLevel)
  data_mod$ind_v <- as.numeric(data_mod$ind_v)
  
  ### Run several epochs of the MCMC to addapt proposal distibutions ###
  accept_rate_adapt = NA
  if(n_epoch_adapt>0){
    opt_accept_rate <- 0.234
    
    # Create matrices for adaptation
    accept_rate_adapt <- theta_ini_adapt <- theta_prop_int_adapt <- matrix(NA,n_epoch_adapt,length(theta_names))
    colnames(accept_rate_adapt) <- colnames(theta_ini_adapt) <- colnames(theta_prop_int_adapt) <- theta_names
    if(gibbs_hm){
      # We dont need to optimize parameters in HM module proposals when using Gibbs sampling there
      # identify names in HM module
      cond = (substr(theta_names,1,7)=="beta_hm") | (substr(theta_names,1,6)=="eta_hm")
      gibbs_hm_param = theta_names[cond]
      # Remove such parameters from adaptive MCMC 
      accept_rate_adapt = accept_rate_adapt[,setdiff(colnames(accept_rate_adapt),gibbs_hm_param)]
    }
    
    # Initializing
    theta_ini_adapt[1,] <- theta_mcmc_ini[colnames(theta_ini_adapt)]
    theta_prop_int_adapt[1,] <- theta_prop_int[colnames(theta_prop_int_adapt)]
    
    for(adapt_i in 1:n_epoch_adapt) {
      # adapt_i <- 1
      # all(names(theta_mcmc_ini)==names(theta_prop_int))
      # if(adapt_i==n_epoch_adapt){browser()}
      mcmc_res <- mcmc_PO_HM_powered( data_arc=as.matrix(data_arc),
                                      data_mod=as.matrix(data_mod),
                                      
                                      power_w_PO=power_w_PO,
                                      power_w_HM=power_w_HM,
                                      
                                      prior_spec_PO=prior_spec_PO,
                                      prior_spec_HM=prior_spec_HM,
                                      
                                      PO_site_rnd_eff=PO_site_rnd_eff,
                                      HM_site_rnd_eff=HM_site_rnd_eff,
                                      
                                      n_iter=n_iter_adapt,
                                      
                                      theta = theta_ini_adapt[1,],
                                      theta_min_max = theta_min_max,
                                      theta_prop_int = theta_prop_int_adapt[adapt_i,],
                                      theta_prop_kernel = theta_prop_kernel,
                                      
                                      ManureLevel_imp=ManureLevel_arc_imp_ini,
                                      Rainfall_imp=Rainfall_arc_imp_ini,
                                      
                                      imp_playpen=imp_playpen,
                                      
                                      gibbs_hm=gibbs_hm,
                                      
                                      keep_imp=FALSE,
                                      keep_ll=FALSE,
                                      
                                      check_mcmc=FALSE,
                                      verbose=FALSE )
      
      colnames(mcmc_res$theta_mcmc) <- theta_names
      
      # acceptance rate
      accept_rate_all <- apply(abs(mcmc_res$theta_mcmc[-1,]-mcmc_res$theta_mcmc[-nrow(mcmc_res$theta_mcmc),])>0,2,mean)
      accept_rate_adapt[adapt_i,] <- accept_rate_all[colnames(accept_rate_adapt)]
      rm(accept_rate_all)
      
      # check that the chain is moving
      not_moving_param <- accept_rate_adapt[adapt_i,]==0
      if(!PO_site_rnd_eff){
        not_moving_param <- not_moving_param[substr(names(not_moving_param),1,6)!="eta_po"]
        not_moving_param <- not_moving_param[names(not_moving_param)!="sigma_eta_PO"]
      }
      if(!HM_site_rnd_eff){
        not_moving_param <- not_moving_param[substr(names(not_moving_param),1,6)!="eta_hm"]
        not_moving_param <- not_moving_param[names(not_moving_param)!="sigma_hm_eta"]
      }
      if( any(not_moving_param,na.rm=TRUE) ){stop("There's a problem with the MCMC (code 1)")}
      rm(not_moving_param)
      
      if( adapt_i < n_epoch_adapt ){
        # theta_ini
        theta_ini_adapt[adapt_i+1,] <- mcmc_res$theta_mcmc[nrow(mcmc_res$theta_mcmc),]
        
        # theta_prop_int
        if(adapt_i==1){
          aux <- log(accept_rate_adapt[adapt_i,]/opt_accept_rate)
          aux <- aux/abs(aux)
          theta_prop_int_adapt[2,names(aux)] <- theta_prop_int_adapt[1,names(aux)] * 2^(aux)
          if(!PO_site_rnd_eff){
            theta_prop_int_adapt[,substr(colnames(theta_prop_int_adapt),1,6)=="eta_po"] <- 1
            theta_prop_int_adapt[,"sigma_eta_PO"] <- 1
          }
          if(!HM_site_rnd_eff){
            theta_prop_int_adapt[,substr(colnames(theta_prop_int_adapt),1,6)=="eta_hm"] <- 1
            theta_prop_int_adapt[,"sigma_hm_eta"] <- 1
          }
        } else {
          for(par_i in seq_along(theta_names)) {
            # par_i=1
            if( is.element( theta_names[par_i] , colnames(accept_rate_adapt) ) ) {
              data_accept_rate = data.frame( y = accept_rate_adapt[1:adapt_i,theta_names[par_i]],
                                             x = theta_prop_int_adapt[1:adapt_i,theta_names[par_i]] )
              aux <- lm( y~x, data_accept_rate)$coeff
              
              # There should be a negative relation between the length of the proposal interval and the acceptance rate
              if( (aux[2]<0)&!is.na(aux[2]) ) {
                theta_prop_int_adapt[adapt_i+1,theta_names[par_i]] <- ( opt_accept_rate - aux[1] ) / aux[2]
              } else {
                # If the relation is different, then vary the size of the interval randomly
                theta_prop_int_adapt[adapt_i+1,theta_names[par_i]] <- theta_prop_int_adapt[adapt_i,theta_names[par_i]] * runif(1,0.5,2)
              }
            }
          }
          # Calculates the percentage change for the interval
          aux <- theta_prop_int_adapt[adapt_i+1,]/theta_prop_int_adapt[adapt_i,]
          # limit the size of such change
          theta_prop_int_adapt[adapt_i+1,!is.na(aux)&(aux>2)] <- theta_prop_int_adapt[adapt_i,!is.na(aux)&(aux>2)]*2
          theta_prop_int_adapt[adapt_i+1,!is.na(aux)&(aux<0.5)] <- theta_prop_int_adapt[adapt_i,!is.na(aux)&(aux<0.5)]*0.5
          
        }
      } else {
        theta_mcmc_ini <- mcmc_res$theta_mcmc[nrow(mcmc_res$theta_mcmc),]
        
        # check that gibbs hm does not generate consequences on addaptation
        theta_prop_int <- theta_prop_int_adapt[nrow(theta_prop_int_adapt),]
        theta_prop_int[colnames(accept_rate_adapt)] <- theta_prop_int_adapt[ cbind(apply(abs(accept_rate_adapt-opt_accept_rate),2,which.min) , match(colnames(accept_rate_adapt),colnames(theta_prop_int_adapt)) ) ]
        # Fill NAs, including parameters that do not use MH (eg updated using Gibbs)
        theta_prop_int[is.na(theta_prop_int)] = 1
      }
      
    }
  }
  
  theta_mcmc_ini <- theta_singleimpute
  if(!is.null(theta_ini)) {
    theta_mcmc_ini[ names(theta_ini) ] <- theta_ini[]
  }
  
  ### RUN THE MAIN MCMC ###
  mcmc_res <- mcmc_PO_HM_powered( data_arc=as.matrix(data_arc),
                                  data_mod=as.matrix(data_mod),
                                  
                                  power_w_PO=power_w_PO,
                                  power_w_HM=power_w_HM,
                                  
                                  prior_spec_PO=prior_spec_PO,
                                  prior_spec_HM=prior_spec_HM,
                                  
                                  PO_site_rnd_eff=PO_site_rnd_eff,
                                  HM_site_rnd_eff=HM_site_rnd_eff,
                                  
                                  n_iter=n_iter,
                                  
                                  theta=theta_mcmc_ini,
                                  theta_min_max=theta_min_max,
                                  theta_prop_int=theta_prop_int,
                                  theta_prop_kernel=theta_prop_kernel,
                                  
                                  ManureLevel_imp=ManureLevel_arc_imp_ini,
                                  Rainfall_imp=Rainfall_arc_imp_ini,
                                  
                                  imp_playpen=imp_playpen,
                                  
                                  gibbs_hm=gibbs_hm,
                                  
                                  keep_imp=keep_imp,
                                  keep_ll=keep_ll,
                                  
                                  check_mcmc=FALSE,
                                  verbose=FALSE )
  colnames(mcmc_res$theta_mcmc) <- theta_names
  
  mcmc_res[["theta_mcmc"]] <- coda::mcmc( mcmc_res[["theta_mcmc"]] )
  if(keep_imp) {
    mcmc_res[["ManureLevel_imp_mcmc"]] <- coda::mcmc( mcmc_res[["ManureLevel_imp_mcmc"]] )
    mcmc_res[["Rainfall_imp_mcmc"]] <- coda::mcmc( mcmc_res[["Rainfall_imp_mcmc"]] )
  }
  # Saving results (preliminary) #
  saveRDS( mcmc_res, file=out_file_rds )
  
  if( is.element(impute_spec,c("smi","cut")) ) {
    
    # Model matrices
    X = matrix(data_arc$Size,ncol=1)
    X_eta = model.matrix(~-1+Site,data=data.frame(Site=as.factor(data_arc$Site)))
    
    comb <- function(...) {
      mapply('rbind', ..., SIMPLIFY=FALSE)
    }
    mcmc_res_smi = foreach( iter_i = 1:n_iter, .combine='comb', .multicombine=TRUE ) %do% {
      # iter_i = 1
      mcmc_res_sub = mcmc_PO( Y = mcmc_res[["ManureLevel_imp_mcmc"]][iter_i,],
                              X = X,
                              X_eta = X_eta,
                              
                              K=3,
                              
                              power_w=1,
                              
                              prior_spec=prior_spec_PO,
                              
                              rnd_eff=PO_site_rnd_eff,
                              
                              n_iter=n_iter_sub,
                              
                              theta = mcmc_res$theta_mcmc[iter_i,theta_names_po],
                              theta_min_max = theta_min_max[theta_names_po,],
                              theta_prop_int = theta_prop_int[theta_names_po],
                              theta_prop_kernel = theta_prop_kernel,
                              
                              keep_ll = keep_ll,
                              
                              check_mcmc = FALSE,
                              verbose = FALSE, quiet=TRUE )
      
      # Only keep the last iteration
      list( mcmc_res_sub$theta_mcmc[n_iter_sub,] ,
            mcmc_res_sub$loglik_mcmc[n_iter_sub,] )
    }
    names(mcmc_res_smi) = c("theta_mcmc","loglik_mcmc")
    colnames(mcmc_res_smi[["theta_mcmc"]]) = theta_names_po
    mcmc_res_smi[["theta_mcmc"]] = coda::mcmc( mcmc_res_smi[["theta_mcmc"]] )
    
    # Substitute second stage samples in the main "theta_mcmc" #
    mcmc_res[["theta_mcmc_tempered"]] = mcmc_res[["theta_mcmc"]][,theta_names_po]
    mcmc_res[["theta_mcmc"]][,theta_names_po] = mcmc_res_smi[["theta_mcmc"]]
    # also for second stage loglik_mcmc #
    mcmc_res[["loglik_mcmc_tempered"]] = mcmc_res[["loglik_mcmc"]][,(1:nrow(data_arc)),1]
    mcmc_res[["loglik_mcmc"]][,(1:nrow(data_arc)),1] = mcmc_res_smi[["loglik_mcmc"]]
  }
  
  # elpd calculation
  if(keep_ll&elpd){
    loo_res <- list()
    loo_res[[1]] <- loo::loo( mcmc_res$loglik_mcmc[,(1:nrow(data_arc)),1] )
    loo_res[[2]] <- loo::loo( mcmc_res$loglik_mcmc[,(1:nrow(data_arc)),2] )
    loo_res[[3]] <- loo::loo( mcmc_res$loglik_mcmc[,-(1:nrow(data_arc)),2] )
    names(loo_res) <- c("PO","HM_arc","HM_mod")
    mcmc_res$loo = loo_res
    
    waic_res <- list()
    waic_res[[1]] <- loo::waic( mcmc_res$loglik_mcmc[,(1:nrow(data_arc)),1] )
    waic_res[[2]] <- loo::waic( mcmc_res$loglik_mcmc[,(1:nrow(data_arc)),2] )
    waic_res[[3]] <- loo::waic( mcmc_res$loglik_mcmc[,-(1:nrow(data_arc)),2] )
    names(waic_res) <- c("PO","HM_arc","HM_mod")
    mcmc_res$waic = waic_res
    
    # summarise predictive measures
    mcmc_res$predict_summary = foreach::foreach(i=1:3,.combine = rbind) %do%{
      loo_summary <- mcmc_res$loo[[i]]$estimates %>% as.data.frame() %>%
        tibble::rownames_to_column("param") %>%
        tidyr::pivot_longer(-param,names_to='statistic',values_to="value") %>%
        dplyr::mutate(param=paste(param,names(mcmc_res$loo)[i],sep="_"))
      waic_summary <- mcmc_res$waic[[i]]$estimates %>% as.data.frame() %>%
        tibble::rownames_to_column("param") %>%
        tidyr::pivot_longer(-param,names_to='statistic',values_to="value") %>%
        dplyr::mutate(param=paste(param,names(mcmc_res$waic)[i],sep="_"))
      rbind(loo_summary,waic_summary) %>%
        as.data.frame()
    }
    
  }
  
  ### Adding information about the run ###
  mcmc_res$n_obs_arc = nrow(data_arc)
  mcmc_res$n_obs_mod = nrow(data_mod)
  mcmc_res$n_obs = mcmc_res$n_obs_arc + mcmc_res$n_obs_mod
  mcmc_res$theta_ini = theta_mcmc_ini
  mcmc_res$theta_singleimpute =theta_singleimpute
  
  mcmc_res$theta_min_max = theta_min_max
  mcmc_res$theta_prop_int = theta_prop_int
  mcmc_res$theta_prop_kernel = theta_prop_kernel
  mcmc_res$impute_spec = impute_spec
  mcmc_res$accept_rate_adapt = accept_rate_adapt
  
  if( !is.null(log_file) ) {
    cat("\n\n---------------------------\n\n",
        "Finish time:\n",as.character(Sys.time()),"\n\n",
        file=log_file, append=TRUE)
  }
  
  if( !is.null(out_file_rds) ) {
    # Saving the results #
    saveRDS( mcmc_res, file=out_file_rds )
  }
  
  return(mcmc_res)
  
}
