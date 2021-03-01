#' @title
#'     Generate Article Carmona & Nicholls (2020)
#'
#' @description
#'     This function reproduce the article ``Semi-Modular Inference: enhanced learning in multi-modular models by tempering the influence of components''.
#'     The function copies an internal .Rnw file to the specified output folder, and then compiles it using the knitr package.
#'
#' @param out_dir Output directory where the pdf and auxiliary latex files will be stored
#' @param compute_mcmc Indicates if the MCMC must be recomputed. If TRUE, the computation time can be in the order of several hours.
#' @param mcmc_dir Directory containing the required files with MCMC results.
#' @param compile_path Directory where the article is compiled and all auxiliary files are saved into. If NULL, a temporary directory is used.
#' @param pkg_dir Internal directory where the package is installed. Used to retrieve the master Rnw file.
#'
#' @return
#'    The function will generate a pdf corresponding to the statistical supplement of Styring et al. (2020)
#'
#' @importFrom utils Sweave
#'
#' @export

generate_article <- function( out_dir,
                              compute_mcmc = FALSE,
                              mcmc_dir = NULL,
                              compile_path = NULL,
                              pkg_dir=system.file(package = "aistats2020smi")
) {

  # If mcmc_dir is not specified, the MCMC results will be downloaded ta a temporary directory and discarded afterwards
  if( !compute_mcmc & is.null(mcmc_dir) ){
    mcmc_dir <- tempdir()
    download_mcmc_results( mcmc_dir = mcmc_dir )
  }

  # If compile_path is not specified, a temporary directory is created and consequently all auxiliary files are discarded
  if( is.null(compile_path) ){
    compile_path = tempdir()
  }

  # Remove trailing slash from directories
  out_dir = normalizePath( out_dir )
  mcmc_dir = normalizePath( mcmc_dir )
  pkg_dir = normalizePath( pkg_dir )
  compile_path = normalizePath( compile_path )

  if( !file.exists(out_dir) ) {
    stop(paste(out_dir,'not found'))
  }

  # Copy Rnw files to compile_path
  tex_path = paste( pkg_dir,'/article',sep='')
  file.copy( paste( tex_path,'/',list.files(tex_path),sep=''), compile_path, overwrite=TRUE )

  # Set variables inside the Rnw file
  rnw_file = "smi_aistats_2020.Rnw"

  # Compile article within the directory compile_path
  cwd = getwd()
  setwd(compile_path)

  cat("Compiling ",rnw_file,"...\n")
  pdf_path = paste(out_dir,'/Carmona and Nicholls - 2020 - Semi-Modular Inference.pdf',sep='')

  tx  <- readLines( rnw_file )
  tx  <- gsub( pattern = "mcmc_dir = 'set_mcmc_dir'",
                replacement = paste("mcmc_dir = '",mcmc_dir,"'",sep=''),
                x = tx)
  tx  <- gsub( pattern = "compute_mcmc = FALSE",
               replacement = paste("compute_mcmc = ",compute_mcmc,sep=''),
               x = tx)
  writeLines(tx, con=rnw_file)

  utils::Sweave( file=rnw_file, encoding='utf-8', quiet=TRUE )
  tools::texi2pdf( gsub('.Rnw','.tex',rnw_file), quiet=TRUE )

  file.copy( from=gsub('.Rnw','.pdf',rnw_file), to=pdf_path )

  # Return to original directory
  setwd(cwd)

  cat('Article and supplement generated succesfully\n',pdf_path,'\n')

  return(TRUE)

}
