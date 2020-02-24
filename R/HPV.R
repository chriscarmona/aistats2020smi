#'
#' @docType data
#'
#' @name
#'     HPV
#'
#' @title
#'     Data for ecological study of HPV and cervical cancer.
#'
#' @description
#'     Data used for investigation of the international correlation between human papillomavirus (HPV) prevalence and cervical cancer incidence (Maucort-Boulch et al. 2008).
#'     The data came from a series of surveys in 13 different countries. Incidence data come from cancer registries in the same populations. We consider the oldest age group (55–64 years) in the survey, which demonstrated the strongest correlation between HPV and cervical cancer.
#'     In population i, the outcome ncases[i] is the number of cancer cases arising from Npop[i] woman-years of follow-up; nhpv[i] is the number of women infected with high-risk HPV in a sample of size Npart[i] from the same population.
#'
#' @usage
#'     HPV
#'
#' @format
#'     A data frame with 13 rows and 4 variables.
#'     \describe{
#'         \item{ncases}{ Number of cases of cervical cancer. }
#'         \item{Npop}{ Population size covered by the cancer registry corresponding to the area of the survey }
#'         \item{nhpv}{ Number of high-risk HPV-positive participants. }
#'         \item{Npart}{ Number of participants in the HPV prevalence survey. }
#'     }
#'
#' @details
#'     A data frame with 13 rows and 4 variables.
#'
#' @examples
#'    head(HPV)
#'    HPV <- data.frame( ncases = c( 16, 215, 362, 97, 76,
#'                                   62, 710, 56, 133,28,
#'                                   62, 413, 194 ),
#'                       Npop = c( 26983, 250930, 829348, 157775, 150467,
#'                                 352445, 553066, 26751, 75815, 150302,
#'                                 354993, 3683043, 507218 ),
#'                       nhpv = c( 7, 6, 10, 10, 1,
#'                                 1, 10, 4, 35, 0,
#'                                 10, 8, 4 ),
#'                       Npart = c( 111, 71, 162, 188, 145,
#'                                  215, 166, 37, 173, 143,
#'                                  229, 696, 93 ) )
#'
#' @seealso
#'    \code{\link{mcmc_hpv}}
#'
#' @references
#'    Martyn Plummer (2015) Cuts in Bayesian graphical models.
#'    Chris U. Carmona, Geoff K. Nicholls (2020). Semi-Modular Inference: enhanced learning in multi-modular models by tempering the influence of components
#'
NULL
