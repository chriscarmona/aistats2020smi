#'
#' @docType data
#'
#' @name
#'     NMeso
#'
#' @title
#'     Archaeological Data from excavations in the north Mesopotamia region
#'
#' @description
#'     Archaeological Data from excavations in the north Mesopotamia region.
#'     Data provided by Amy Bogaard <amy.bogaard@arch.ox.ac.uk>
#'
#' @format
#'     A data frame with 269 rows and 15 variables.
#'     \describe{
#'         \item{ID}{cereal grain-sample ID.}
#'         \item{normd13C}{unknown.}
#'         \item{d13Csd}{unknown.}
#'         \item{d13Cair}{unknown.}
#'         \item{D13C}{unknown.}
#'         \item{normd15N}{Measured cereal grain-sample nitrogen isotope ratio.}
#'         \item{d15Nsd}{unknown.}
#'         \item{Site}{The context site name.}
#'         \item{Phase}{The context phase name.}
#'         \item{Date}{The context date.}
#'         \item{Rainfall}{The annual rainfall base estimate for the grain-sample context.}
#'         \item{Size}{The size of the site in hectares at the time of the context date.}
#'         \item{Species}{The cereal species.}
#'         \item{Category}{The species category: barley, wheat, etc.}
#'         \item{minDate}{lower bound on the unknown true date for the grain-sample context.}
#'         \item{maxDate}{upper bound on the unknown true date for the grain-sample context.}
#'         \item{maxRainfall}{upper bound on the unknown true Rainfall for the grain-sample context. Need to be added to \code{Rainfall}}
#'         \item{minRainfall}{upper bound on the unknown true Rainfall for the grain-sample context. Need to be added to \code{Rainfall}}
#'         \item{dataset}{Name of the dataset. Used to distinguish among data from other sources.}
#'     }
#'
#' @details
#'     A data frame with 269 rows and 19 variables.
#'
#' @examples
#'
#'    head(NMeso)
#'
#' @seealso
#'    \code{\link{mcmc_agric_model}}
#'

NULL
