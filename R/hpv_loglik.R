
#'
#' @import nlme
#' @import ordinal
#' @import doRNG
#'
#' @export
hpv_loglik <- function( data,
                        theta, phi,
                        eta_modules=setNames(c(1,1), c("pois", "binom")) ) {

  # Calculates loglikelihood of HPV model

  # browser()
  # ncases[j] ~ dpois(mu[j])
  # mu[j] <- Npop[j] * 1.0E-3 * exp( theta[1] + phi[j] * theta[2])
  # nhpv[j] ~ dbin( phi[j], Npart[j])


  # data: data frame with n_obs rows and four columns: "nhpv", "Npart", "ncases", "Npop"
  # theta: matrix of dimension c(2, 1)
  # phi: matrix of dimension c(n_obs, 1)
  # eta_modules: vector of dimension 2

  n_obs = nrow(data)

  check_dim <- (length(phi)==n_obs) && (length(theta)==2)
  if(!check_dim){ stop("There is a problem with the dimensions of theta and/or phi") }

  mu <- data[,"Npop"] * 1.0E-3 * exp( theta[1] + theta[2] * phi  )

  loglik <- matrix(0,nrow=n_obs,ncol=2)
  loglik[,1] <- eta_modules["pois"] * dpois( x=data[,"ncases"], lambda=c(mu), log=TRUE )
  loglik[,2] <- eta_modules["binom"] * dbinom( x=data[,"nhpv"], size=data[,"Npart"], prob=phi , log=TRUE )
  colnames(loglik) <- c("pois","binom")

  return( loglik )
}
