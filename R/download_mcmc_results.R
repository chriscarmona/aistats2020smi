#' @title
#'     Download MCMC results from Carmona and Nicholls (2020)
#'
#' @description
#'     Downloads the MCMC results associated with the article Carmona and Nicholls (2020).
#'     The function downloads a .zip file from the authors website and store it in the specified output folder.
#'
#' @param mcmc_dir Output directory where the files with the MCMC will be stored.
#' @param url a character string naming the URL of a resource to be downloaded
#'
#' @return
#'    The function will generate a pdf corresponding to the statistical supplement of Styring et al. (2020)
#'
#' @importFrom utils download.file unzip
#'
#' @export

download_mcmc_results <- function( mcmc_dir,
                                   url = "https://chriscarmona.me/files/aistats2020smi_mcmc_results.zip"
  ) {

  mcmc_dir = normalizePath(mcmc_dir)
  cat("Downloading MCMC files to:\n",mcmc_dir,"\n")
  temp <- tempfile()
  utils::download.file( url=url, destfile=temp )
  utils::unzip(temp, exdir=mcmc_dir, junkpaths=TRUE )

  return(TRUE)
}
