
#' @title Joint log-likelihood for the HPV model
#'
#' @description Computes the joint log-likelihood for the HPV model
#'
#' @param data data frame with n_obs rows and four columns: "nhpv", "Npart", "ncases", "Npop"
#' @param theta matrix of dimension c(2, 1). Two parameters defining the linear relation between HPV prevalence and cancer incidence in the HPV model.
#' @param phi matrix of dimension c(n_obs, 1). The HPV prevalence parameter of each population. values in (0,1)
#' @param eta_modules Named tuple specifying the degree of influence for the "pois" and "binom" modules in the HPV model.
#' @importFrom stats dbinom dpois setNames
#' @export
hpv_loglik <- function( data,
                        theta, phi,
                        eta_modules=stats::setNames(c(1,1), c("pois", "binom")) ) {

  n_obs = nrow(data)

  check_dim <- (length(phi)==n_obs) && (length(theta)==2)
  if(!check_dim){ stop("There is a problem with the dimensions of theta and/or phi") }

  mu <- data[,"Npop"] * 1.0E-3 * exp( theta[1] + theta[2] * phi  )

  loglik <- matrix(0,nrow=n_obs,ncol=2)
  loglik[,1] <- eta_modules["pois"] * stats::dpois( x=data[,"ncases"], lambda=c(mu), log=TRUE )
  loglik[,2] <- eta_modules["binom"] * stats::dbinom( x=data[,"nhpv"], size=data[,"Npart"], prob=phi , log=TRUE )
  colnames(loglik) <- c("pois","binom")

  return( loglik )
}
