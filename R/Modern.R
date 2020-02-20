#'
#' @docType data
#'
#' @name
#'     Modern
#'
#' @title
#'     Agricultural Data gathered in an experimental setting
#'
#' @description
#'     The Modern data were gathered in an experimental setting in which manure level is an independent variable and the isotope ratio is an experiment outcome.
#'     Data provided by Amy Bogaard <amy.bogaard@arch.ox.ac.uk>
#'
#' @usage
#'     Modern
#'
#' @format
#'     A data frame with 268 rows and 13 variables.
#'     \describe{
#'         \item{ID}{cereal grain-sample ID.}
#'         \item{normd13C}{unknown.}
#'         \item{d13Cair}{unknown.}
#'         \item{D13C}{unknown.}
#'         \item{normd15N}{Measured cereal grain-sample nitrogen isotope ratio.}
#'         \item{Site}{The grain-sample site name.}
#'         \item{Rainfall}{Annual rainfall derived from interpolation of average monthly climate data for 1960-1990, available from the WorldClim database, for the location in which the cereal grain sample was grown.}
#'         \item{Species}{The cereal species.}
#'         \item{Category}{The species category: barley, wheat, etc.}
#'         \item{Year}{The year in which the grain-sample was gathered.}
#'         \item{Irrigation}{A string describing the irrigation type.}
#'         \item{ManureLevel}{The level of the manuring treatment the cereal grain-sample received. A 3 level ordinal variable with levels “low”, “medium” and “high”}
#'         \item{dataset}{Name of the dataset. Used to distinguish among data from other sources.}
#'     }
#'
#' @details
#'     A data frame with 268 rows and 13 variables.
#'
#' @examples
#'
#'    head(Modern)
#'
#' @seealso
#'    \code{\link{mcmc_agric_model}}
#'

NULL
