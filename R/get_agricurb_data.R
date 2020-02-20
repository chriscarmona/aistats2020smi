#' @title
#'     Function to obtain the real agricultural data in a convenient format
#'
#' @description
#'     The function \code{get_agricurb_data} pre-process the agricultural data into a convenient format ready to analisis.
#'
#'     The formatting process is as follows:
#' \enumerate{
#'     \item Keep modern data and as many archaeological datasets as specified in \code{arc_datasets}.
#'     \item Eliminate rows with unknown \code{Size} in the archaeological datasets.
#'     \item Eliminate rows with \code{Category} different to \code{wheat} and \code{barley}.
#'     \item Modify \code{Rainfall} in the archaeological data sets. Creates \code{Rainfall_min=Rainfall+minRainfall} and \code{Rainfall_max=Rainfall+maxRainfall}, and after that \code{Rainfall=NA}.
#'     \item Transform the variables listed in \code{vars_log} to a logaritmic scale.
#'     \item Rescale the variables listed in \code{vars_scale_int_0_1} to the interval \code{[0,1]}.
#'     \item Rescale the variables listed in \code{vars_scale_mean_0_var_1} to have empirical mean 0 and empirical variance 1.
#' }
#'
#' @param arc_datasets indicated the datasets used for sampling the archaeological covariates.
#' @param vars_log indicated the variables to be rescale to have mean \code{0} and variance \code{1}.
#' @param vars_scale_int_0_1 indicated the variables to be rescale to the interval \code{[0,1]}.
#' @param vars_scale_mean_0_var_1 indicated the variables to be rescale to have mean \code{0} and variance \code{1}.
#'
#' @return
#'    Returns a data frame similar to AgricurbAll ready to analyse with the current methods, such as \code{mcmc_PO1_HM1_HM1}.
#'
#' @examples
#'
#' require(ggplot2)
#'
#' Agricurb_data <- get_agricurb_data( arc_datasets = "NMeso",
#'                                          vars_log = c("Rainfall"),
#'                                          vars_scale_int_0_1 = c("Size", "Rainfall"),
#'                                          vars_scale_mean_0_var_1 = NULL )
#' Agricurb_data$ManureLevel <- factor(Agricurb_data$ManureLevel,levels=c("low","medium","high"))
#'
#' # archaeological data identifier
#' id_arch <- !is.element(Agricurb_data$dataset,"modern")
#'
#' #####
#' # Checking data characteristics
#' #####
#'
#' # HM model #
#' # Relation between Rainfall and normd15N
#' # (assuming missing data)
#' p <- ggplot( aes(x=Rainfall,y=normd15N), data=Agricurb_data)
#' p <- p + geom_point(aes(colour=ManureLevel), pch=20, size=3, alpha=0.7)
#' p
#'
#' @import ggplot2
#' @import stats
#'
#' @export
#'

get_agricurb_data <- function( arc_datasets="NMeso",
                                    vars_log = NULL,
                                    vars_sqrt = NULL,
                                    vars_scale_int_0_1 = NULL,
                                    vars_scale_mean_0_var_1 = NULL ) {

  # Starts with the Modern dataset
  Agricurb_data <- aistats2020smi::Modern

  # Merges archaeological datasets
  for (i in 1:length(arc_datasets)) {
    Agricurb_data <- merge( Agricurb_data, eval( parse(text = paste("aistats2020smi::",arc_datasets[i],sep="")) ), all=T )
  }

  # Eliminates rows with unknown "Size" in the archaeological datasets.
  aux_vec <- !(!is.element(Agricurb_data$dataset,c("modern")) & is.na(Agricurb_data$Size))
  Agricurb_data <- Agricurb_data[aux_vec,]

  # Eliminates rows with "Category" different to "wheat" or "barley."
  aux_vec <- is.element(Agricurb_data$Category,c("wheat","barley"))
  Agricurb_data <- Agricurb_data[aux_vec,]

  # Indicator variable for rows with archaelogical info
  id_arch <- !is.element(Agricurb_data$dataset,c("modern"))


  # Calculates min and max values of "Rainfall" in the archaeological data sets.
  # Creates:
  #    Rainfall_min = Rainfall + minRainfall
  #    Rainfall_max = Rainfall + maxRainfall
  # and after that Rainfall=NA
  Agricurb_data[id_arch,"Rainfall_guess"] <- Agricurb_data[id_arch,"Rainfall"]
  Agricurb_data[id_arch,"Rainfall"] <- NA
  Agricurb_data$Rainfall_max <- Agricurb_data$Rainfall_guess + Agricurb_data$maxRainfall
  Agricurb_data$Rainfall_min <- Agricurb_data$Rainfall_guess + Agricurb_data$minRainfall


  # Transform the variables listed in vars_log to a logaritmic scale: log(x+1)
  # vars_log = c("Rainfall")
  if(!is.null(vars_log)){
    if(is.element("Rainfall",vars_log) ) {
      vars_log <- union(vars_log,c("Rainfall_min","Rainfall_max"))
    }
    Agricurb_data[,vars_log] <- log(1+Agricurb_data[,vars_log])
  }

  # Transform the variables listed in vars_sqrt, applying a square root.
  # vars_sqrt = c("Size")
  if(!is.null(vars_sqrt)){
    if(is.element("Rainfall",vars_sqrt) ) {
      vars_sqrt <- union(vars_sqrt,c("Rainfall_min","Rainfall_max"))
    }
    Agricurb_data[,vars_sqrt] <- sqrt(Agricurb_data[,vars_sqrt])
  }
  
  # Rescale the variables listed in "vars_scale_int_0_1" to the interval [0,1]
  # vars_scale_int_0_1 = c("Size","Rainfall")
  if( !is.null(vars_scale_int_0_1) ) {
    if( is.element("Rainfall",vars_scale_int_0_1) ) {
      # different treatment for rainfall, it has many variables associated

      # Rescaling Rainfall
      aux_vec1 <- min(Agricurb_data[,c("Rainfall","Rainfall_min","Rainfall_max")],na.rm=T)
      aux_vec2 <- max(Agricurb_data[,c("Rainfall","Rainfall_min","Rainfall_max")],na.rm=T)
      Agricurb_data[,c("Rainfall","Rainfall_min","Rainfall_max")] <- (Agricurb_data[,c("Rainfall","Rainfall_min","Rainfall_max")]-aux_vec1)/(aux_vec2-aux_vec1)
      rm(aux_vec1,aux_vec2)
      vars_scale_int_0_1 <- setdiff(vars_scale_int_0_1,"Rainfall")
    }

    aux_vec1 <- apply(Agricurb_data[,vars_scale_int_0_1,drop=F], 2, min,na.rm=T)
    aux_vec2 <- apply(Agricurb_data[,vars_scale_int_0_1,drop=F], 2, max,na.rm=T)
    aux_mat1 <- matrix(aux_vec1,nrow=nrow(Agricurb_data),ncol=length(aux_vec1),byrow=T)
    aux_mat2 <- matrix(aux_vec2-aux_vec1,nrow=nrow(Agricurb_data),ncol=length(aux_vec1),byrow=T)

    Agricurb_data[,vars_scale_int_0_1] <- (Agricurb_data[,vars_scale_int_0_1,drop=F] - aux_mat1)/aux_mat2
  }

  # Rescale the variables listed in "vars_scale_mean_0_var_1" to have empirical mean 0 and empirical variance 1.
  # vars_scale_mean_0_var_1 = c("normd15N")
  if( !is.null(vars_scale_mean_0_var_1) ) {
    aux_vec1 <- apply(Agricurb_data[,vars_scale_mean_0_var_1,drop=F], 2, mean,na.rm=T)
    aux_vec2 <- apply(Agricurb_data[,vars_scale_mean_0_var_1,drop=F], 2, stats::sd,na.rm=T)
    if(is.element("Rainfall",vars_scale_mean_0_var_1) ) {
      vars_scale_mean_0_var_1 <- union(vars_scale_mean_0_var_1,c("Rainfall_min","Rainfall_max"))
      aux_vec1[match(c("Rainfall_min","Rainfall_max"),vars_scale_mean_0_var_1)] <- aux_vec1[match("Rainfall",vars_scale_mean_0_var_1)]
      aux_vec2[match(c("Rainfall_min","Rainfall_max"),vars_scale_mean_0_var_1)] <- aux_vec2[match("Rainfall",vars_scale_mean_0_var_1)]
    }
    aux_mat1 <- matrix(aux_vec1,nrow=nrow(Agricurb_data),ncol=length(aux_vec1),byrow=T)
    aux_mat2 <- matrix(aux_vec2,nrow=nrow(Agricurb_data),ncol=length(aux_vec2),byrow=T)

    Agricurb_data[,vars_scale_mean_0_var_1] <- (Agricurb_data[,vars_scale_mean_0_var_1,drop=F] - aux_mat1)/aux_mat2
  }

  # Keep only relevant columns
  Agricurb_data <- Agricurb_data[,c("dataset","Site","Phase","Date","Size","ManureLevel","Rainfall","Rainfall_min","Rainfall_max",'Rainfall_guess',"Category","normd15N")]

  return( Agricurb_data )
}
