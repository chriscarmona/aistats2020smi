#' @title
#'     Function to generate synthetic agricultural data according to the model PO1_HM1_HM1
#'
#' @description
#'     \code{get_synth_agric_data} generates synthetic data using the specified model parameters and data observed in covariates from archaelogical and modern data.
#'
#' @param theta data.frame/matrix with one row. Defining the parameters in the model. colnames(theta) should specify the parameter.
#' @param arc_datasets indicated the datasets used for sampling the archaeological covariates
#' @param covariates Character. One of c("true","nice"), specyfing the nature of covariates in the synthetic data.
#'
#' @details
#' The hierarchical model consists of two parts: The Proportional Odds (PO1) component, and the Gaussian Linear model (HM1).
#'   PO1:
#'      ManureLevel ∼ 1 + log(Size) + (1|Site)
#'      link = "logistic"
#'   HM1:
#'      normd15N ∼ 1 + log(Rainfall) + ManureLevel + (1|Site)
#'      weights = varIdent(form=~1|Category)
#'
#' @return
#' \code{get_synth_agric_data} returns a list with the following objects:
#' \describe{
#'     \item{\code{synth_agric_data}}{ A dataframe similar to AgricurbAll, but the response obbeys the PO1_HM1_HM1 model..}
#'     \item{\code{synth_agric_data_complete}}{ \code{synth_agric_data} without missing data..}
#'     \item{\code{ManureLevel_probs}}{Probabilities of observing each ManureLevel}
#' }
#'
#' @examples
#'
#' rnorm(1)
#'
#' @import magrittr
#'
#' @importFrom reshape melt
#' @importFrom stats rnorm runif
#'
#' @export
#'

get_synth_agric_model <-function( theta=NULL,
                                  arc_datasets = c("NMeso"),
                                  covariates = c("true","nice")[1],
                                  vars_log = NULL,
                                  vars_sqrt = NULL,
                                  vars_scale_int_0_1 = NULL,
                                  vars_scale_mean_0_var_1 = NULL ) {

  # on.exit(browser())
  # Loading real data

  synth_agric_data <- get_agricurb_data( arc_datasets=arc_datasets,
                                         vars_log = vars_log,
                                         vars_sqrt = vars_sqrt,
                                         vars_scale_int_0_1=vars_scale_int_0_1,
                                         vars_scale_mean_0_var_1=vars_scale_mean_0_var_1 )

  # number of observations
  n_obs <- nrow(synth_agric_data)

  # modern data rows
  id_arch <- synth_agric_data$dataset!="modern"

  if(covariates=="nice") {
    # Relation between Site and Size in true data
    synth_agric_data %>%
      dplyr::filter(dataset!="modern") %>%
      ggplot() +
      geom_point(aes(x=Site,y=Size))

    synth_agric_data %>%
      dplyr::filter(dataset=="modern") %>%
      ggplot() +
      geom_point(aes(x=Site,y=Rainfall))

    synth_agric_data[id_arch,"Size"] <- stats::rnorm( sum(id_arch) )
    synth_agric_data[,"Rainfall"] <- stats::rnorm( nrow(synth_agric_data) )
  }

  # indicator of variance offset
  synth_agric_data$ind_v <- NA
  synth_agric_data[synth_agric_data$Category=="wheat","ind_v"] <- 1
  synth_agric_data[synth_agric_data$Category=="barley","ind_v"] <- 0

  ##### Parameters #####

  ManureLevels <- c("low","medium","high")
  k.M <- length(ManureLevels) # number of categories in "ManureLevels"
  Sites <- sort(unique(synth_agric_data$Site))
  Sites_arc <- sort(unique(synth_agric_data$Site[synth_agric_data$dataset!="modern"]))
  n_Site <- length(Sites) # number of categories in "Sites"

  theta_names <- c( paste("alpha_po_",1:(k.M-1),sep=""), # "alpha_po" intercepts of the PO model
                    "gamma_po_1", # "gamma_po_1" fixed effect of PO model. Size on ManureLevel
                    "sigma_eta_PO", # "sigma_eta_PO" variance of random effects of Site on ManureLevel
                    paste("eta_po_",Sites_arc,sep=""), # random effects of Site on ManureLevel
                    paste("beta_hm_",1:(1+1+(k.M-1)),sep=""), # "beta_hm" fixed effects of the HM model.
                    "sigma_hm", # HM model variance
                    "v", # HM variance offset for ind_v rows
                    "sigma_hm_eta", # variance of random effects for Site in the HM model
                    paste("eta_hm_",Sites,sep="") ) # random effects for Site in the HM model

  # Load default parameter values
  theta_default <- aistats2020smi::theta_agric_model_example

  if( !is.null(theta) ) {
    theta_default[names(theta)] <- theta
    theta = theta_default
  }

  if( !all( is.element(theta_names,names(theta)) ) ){
    stop("theta is not correctly specified")
  }

  if(covariates=="true") {
    # Imputing missing rainfall
    synth_agric_data[is.na(synth_agric_data$Rainfall),"Rainfall"] <- apply( synth_agric_data[is.na(synth_agric_data$Rainfall),c("Rainfall_min","Rainfall_max")],1,function(x){ stats::runif(1,as.numeric(x)[1],as.numeric(x)[2]) })
  }

  ##### Response #####

  ### PO1 ###
  X_aux <- as.matrix(synth_agric_data[id_arch,"Size",drop=F])
  Z_aux <- categ_to_dummie(match(synth_agric_data[id_arch,"Site"],Sites_arc),1:length(Sites_arc))

  # synth_agric_data$Site <- factor(synth_agric_data$Site,levels=Sites)
  # X_aux <- model.matrix(~-1+Size,data=synth_agric_data[id_arch,])
  # Z_aux <- model.matrix(~-1+Site,data=synth_agric_data[id_arch,])

  ManureLevel_probs <- matrix(NA,sum(id_arch),k.M)

  for( i in 1:sum(id_arch) ) {
    # i<-61
    for( j in 1:k.M ) {
      # j<-1
      ManureLevel_probs[i,j] <- exp( loglik_PO_cpp( Y = as.matrix(match(ManureLevels[j],ManureLevels)),
                                                    X = cbind(X_aux,Z_aux)[i,,drop=F],
                                                    alpha = theta[paste("alpha_po_",1:(k.M-1),sep="")],
                                                    beta = theta[c( "gamma_po_1",paste("eta_po_",Sites_arc,sep="") )] ) )
    }
  }
  rm(X_aux,Z_aux)

  ManureLevel_probs_allobs <- matrix(NA,n_obs,k.M)
  ManureLevel_probs_allobs[!id_arch,] <- categ_to_dummie(match(synth_agric_data[!id_arch,"ManureLevel"],ManureLevels),1:k.M)
  ManureLevel_probs_allobs[id_arch,] <- ManureLevel_probs

  if(!all(apply(ManureLevel_probs_allobs,1,sum)-1<1e-7)) {stop("There is a problem simulating from the PO model")}

  synth_agric_data[id_arch,"ManureLevel"] <- apply( X=ManureLevel_probs,
                                                     MARGIN = 1,
                                                     FUN = function(x){sample(ManureLevels,size=1,prob=x)} )
  rm(i,j)

  #rm(ManureLevel_probs)
  #table(synth_agric_data$ManureLevel)
  #table(as.numeric(synth_agric_data$ManureLevel))

  ### HM1 ###

  X_aux <- as.matrix( cbind( 1,
                             synth_agric_data[,"Rainfall",drop=F],
                             categ_to_dummie( match(synth_agric_data[,"ManureLevel"],ManureLevels),1:k.M)[,-1] ))
  Z_aux <- categ_to_dummie(match(synth_agric_data[,"Site"],Sites),1:n_Site)
  #X_aux <- model.matrix(~Rainfall+ManureLevel,data=synth_agric_data)
  #Z_aux <- model.matrix(~-1+Site,data=synth_agric_data)

  synth_agric_data$normd15N <- NA
  #synth_agric_data$normd15N_exp <- X_aux%*%t(theta[1,paste("beta_hm_",1:(1+1+(k.M-1)),sep=""),drop=F]) + Z_aux%*%t(theta[1,paste("eta_hm_",Sites,sep=""),drop=F])
  synth_agric_data$normd15N_exp <- cbind(X_aux,Z_aux) %*% matrix( theta[c(paste("beta_hm_",1:(1+1+(k.M-1)),sep=""),paste("eta_hm_",Sites,sep=""))], ncol=1 )

  rm(X_aux,Z_aux)

  # sum(is.na(synth_agric_data$normd15N_exp))

  synth_agric_data[synth_agric_data$ind_v==0,"normd15N"] <- stats::rnorm( n=sum(synth_agric_data$ind_v==0,na.rm=T),
                                                                          mean=synth_agric_data[synth_agric_data$ind_v==0,"normd15N_exp"],
                                                                          sd=theta["sigma_hm"] )
  synth_agric_data[synth_agric_data$ind_v==1,"normd15N"] <- stats::rnorm( n=sum(synth_agric_data$ind_v==1,na.rm=T),
                                                                          mean=synth_agric_data[synth_agric_data$ind_v==1,"normd15N_exp"],
                                                                          sd=sqrt(theta["v"]) * theta["sigma_hm"] )


  synth_agric_data <- synth_agric_data[,c("dataset","Site","Size","ManureLevel","Rainfall","Rainfall_min","Rainfall_max","Category","normd15N")]
  ### Synthetic Missing data ###

  synth_agric_data_complete <- synth_agric_data

  # Archaelogical
  synth_agric_data[id_arch,"Rainfall"] <- NA
  synth_agric_data[id_arch,"ManureLevel"] <- NA
  # Modern
  synth_agric_data[!id_arch,"Size"] <- NA

  return( list( synth_agric_data=synth_agric_data,
                synth_agric_data_complete=synth_agric_data_complete,
                ManureLevel_probs=ManureLevel_probs_allobs,
                theta=theta ) )
}
