rm(list = ls())
options(scipen=999, stringsAsFactors=FALSE)
set.seed(0)

# pkgbuild::compile_dll()
# Rcpp::compileAttributes()
# devtools::document()
# devtools::install()

# loading required packages #
req.pck <- c( "aistats2020smi",
              "MASS","nlme","coda","ordinal","LearnBayes",
              "plyr","tidyverse", "ggplot2","GGally","lattice",
              "foreach","doParallel","doRNG",
              "GPfit","rBayesianOptimization",
              "RColorBrewer","knitr","ggmcmc","cowplot","gridExtra","ggcorrplot" )
req.pck_bool <- sapply(X=req.pck,FUN=require,character.only=T)
if(!all(req.pck_bool)) {
  sapply(X=req.pck[!req.pck_bool],FUN=install.packages,character.only=T);
  sapply(X=req.pck,FUN=require,character.only=T)
}

out_dir <- './inst'

# Parallel processing
parallel_comp = TRUE
if(parallel_comp){
  n_cores = 6
  options(cores=n_cores)
  doParallel::registerDoParallel()
  getDoParWorkers()
}

# Data from Plummer (2014)
aistats2020smi::HPV %>% head()
# HPV <- data.frame(
#   ncases = c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194), # Y
#   Npop = c(26983, 250930, 829348, 157775, 150467, 352445, 553066, 26751, 75815, 150302, 354993, 3683043, 507218), # T
#   nhpv = c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4), # Z
#   Npart = c(111, 71, 162, 188, 145, 215, 166, 37, 173, 143, 229, 696, 93) ) # N
# save(HPV,file="./data/HPV.rda")

# Levels of influence in SMI
power_eta_all <- c( 0.01,0.99,
                    # seq(0.8,0.9,by=0.02),
                    seq(0.10,1.00,by=0.10) )
power_eta_all <- sort( unique(round(power_eta_all,6)) )

hpv_pc_mcmc <- foreach( eta_i = seq_along(power_eta_all) ) %do% {
  # eta_i <- 1
  eta_pois <- power_eta_all[eta_i]
  cat(eta_pois,", ")
  file_i <- paste(out_dir,"/HPV_partial_cut_stan_",formatC(eta_pois,digits=3,format="f",flag="0"),".rds",sep="")
  if(!file.exists(file_i)){
    set.seed(0)
    hpv_mcmc_smi <- aistats2020smi::mcmc_hpv( HPV=aistats2020smi::HPV,

                                              # Number of iterations
                                              n_iter_mcmc = 500, # main chain
                                              n_iter_warmup = 50,
                                              n_chains_mcmc = 4, # Number of chains
                                              n_iter_mcmc_stage_2 = 500, # Subchain

                                              # Cut rate
                                              eta_pois = eta_pois,
                                              eta_binom = 1,

                                              mcmc_file = file_i,
                                              n_cores=n_cores )
  } else {
    hpv_pc_mcmc_i <- readRDS(file=file_i)
  }
  hpv_pc_mcmc_i
}
