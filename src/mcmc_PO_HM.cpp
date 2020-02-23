#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "aistats2020smi_shared.h"

// [[Rcpp::interfaces(r, cpp)]]

arma::colvec loglik_prior_PO_HM( const std::string prior_spec_PO,
                                 const std::string prior_spec_HM,
                                 
                                 const arma::colvec Y_PO,
                                 const arma::mat X_PO,
                                 const arma::mat X_eta_PO,
                                 
                                 const arma::colvec alpha_PO,
                                 const arma::colvec gamma_PO,
                                 const arma::colvec eta_PO,
                                 const double sigma_eta_PO,
                                 
                                 const arma::colvec Y_HM,
                                 const arma::mat X_HM,
                                 const arma::mat X_eta_HM,
                                 const arma::colvec ind_v_HM,
                                 
                                 const arma::colvec beta_HM,
                                 const arma::colvec eta_HM,
                                 const double sigma_HM,
                                 const double v_HM,
                                 const double sigma_eta_HM,
                                 
                                 const double power_w_PO=1,
                                 const double power_w_HM=1 ) {
  
  // output: probs
  arma::vec probs = arma::zeros<arma::vec>(4);
  
  // "probs" contents:
  // 0: loglik_PO
  // 1: logprior_PO
  // 2: loglik_HM
  // 3: logprior_HM
  
  
  probs[0] = power_w_PO * loglik_PO_cpp( Y_PO,
                                         join_horiz(X_PO,X_eta_PO),
                                         alpha_PO,
                                         join_vert(gamma_PO,eta_PO) );
  
  probs[1] = logprior_PO_cpp( alpha_PO,
                              gamma_PO,
                              sigma_eta_PO,
                              eta_PO,
                              prior_spec_PO );
  
  
  
  probs[2] = power_w_HM * loglik_HM_cpp( Y_HM,
                                         join_horiz(X_HM,X_eta_HM),
                                         join_vert(beta_HM,eta_HM),
                                         sigma_HM,
                                         v_HM,
                                         ind_v_HM );
  
  probs[3] = logprior_HM_cpp( beta_HM,
                              sigma_HM,
                              v_HM,
                              sigma_eta_HM,
                              eta_HM,
                              prior_spec_HM );
  
  
  return probs;
  
}



arma::mat loglik_i_PO_HM( const arma::colvec Y_PO,
                          const arma::mat X_PO,
                          const arma::mat X_eta_PO,
                          
                          const arma::colvec alpha_PO,
                          const arma::colvec gamma_PO,
                          const arma::colvec eta_PO,
                          const double sigma_eta_PO,
                          
                          const arma::colvec Y_HM,
                          const arma::mat X_HM,
                          const arma::mat X_eta_HM,
                          const arma::colvec ind_v_HM,
                          
                          const arma::colvec beta_HM,
                          const arma::colvec eta_HM,
                          const double sigma_HM,
                          const double v_HM,
                          const double sigma_eta_HM ) {
  
  const unsigned int n_obs = X_HM.n_rows;
  const unsigned int n_obs_arc = X_PO.n_rows;
  const unsigned int n_obs_mod = n_obs-n_obs_arc;
  
  // output: loglik_i
  arma::mat loglik_i = arma::zeros<arma::mat>(n_obs,2);
  
  loglik_i.submat(0,0,n_obs_arc-1,0) = loglik_PO_i_cpp( Y_PO,
                  join_horiz(X_PO,X_eta_PO),
                  alpha_PO,
                  join_vert(gamma_PO,eta_PO) );
  
  loglik_i.col(1) = loglik_HM_i_cpp( Y_HM,
               join_horiz(X_HM,X_eta_HM),
               join_vert(beta_HM,eta_HM),
               sigma_HM,
               v_HM,
               ind_v_HM );
  
  return loglik_i;
  
}



bool check_po_separability( const arma::vec Y,
                            const arma::vec X,
                            const unsigned int K=3 ){
  
  // Returns true in there is (quasi-)complete separability
  
  if( Y.n_rows != X.n_rows ) {
    throw std::range_error("X and Y do not have the same size!");
  }
  
  bool bool_separable=false;
  unsigned int i=0;
  
  arma::vec aux_vec;
  bool bool_overlap=true;
  
  arma::vec len_range_X_Y = arma::zeros<arma::vec>(K);
  
  arma::mat range_X_Y = arma::zeros<arma::mat>(K,2);
  
  // Ensuring values for all ManureLevels
  aux_vec = arma::unique( Y );
  if( aux_vec.n_rows<K ) {
    bool_separable = true;
  } else {
    // Avoiding Complete Separation //
    // Range of Size for each value of ManureLevels
    for( i=0; i<K; i++ ) {
      range_X_Y(i,0) = X.rows( find( Y ==(i+1)) ).min();
      range_X_Y(i,1) = X.rows( find( Y ==(i+1)) ).max();
    }
    len_range_X_Y = range_X_Y.col(1)-range_X_Y.col(0);
    
    // bool_overlap = any( len_range_X_Y==0 ) | (( range_X_Y.max()-range_X_Y.min() ) > arma::sum( len_range_X_Y ));
    
    // overlap for consecutive classes
    // bool_overlap = true;
    for( i=1; i<K; i++ ) {
      bool_overlap = bool_overlap & ( ( range_X_Y.rows(arma::span(i-1,i)).max()-range_X_Y.rows(arma::span(i-1,i)).min() ) < arma::sum( len_range_X_Y.rows(arma::span(i-1,i)) ) );
    }
    // If there is overlap = No (quasi-)complete separation
    bool_separable = !bool_overlap;
    
  }
  
  return bool_separable;
  
}



// [[Rcpp::export]]
Rcpp::List mcmc_PO_HM_powered( const arma::mat data_arc,
                               const arma::mat data_mod,
                               
                               double power_w_PO,
                               double power_w_HM,
                               
                               const std::string prior_spec_PO,
                               const std::string prior_spec_HM,
                               
                               const bool PO_site_rnd_eff,
                               const bool HM_site_rnd_eff,
                               
                               const unsigned int n_iter,
                               const unsigned int n_warmup,
                               const unsigned int n_thin,
                               
                               arma::colvec theta,
                               const arma::mat theta_min_max,
                               arma::colvec theta_prop_int,
                               const std::string theta_prop_kernel,
                               
                               arma::colvec ManureLevel_imp,
                               arma::colvec Rainfall_imp,
                               
                               const bool imp_playpen=false,
                               
                               const bool gibbs_hm=true,
                               
                               const bool keep_imp=true,
                               const bool keep_ll=true,
                               
                               const bool check_mcmc=false,
                               const bool verbose=false ) {
  /*
   The matrix "data_arc" must contain the columns:
   0: Size
   1: Site
   2: Rainfall_min
   3: Rainfall_max
   4: ind_v
   5: normd15N
   
   The matrix "data_mod" must contain the columns:
   0: Site
   1: ManureLevel
   2: Rainfall
   3: ind_v
   4: normd15N
   */
  
  /*
   OUTPUT:
   theta_mcmc : matrix with the chain of parameters. dim(theta_mcmc)=(n_iter,n_par)
   ManureLevel_imp_mcmc : if keep_imp=true, matrix with the chain of imputed ManureLevel for archaeological data. dim(theta_mcmc)=(n_iter,n_obs_arc)
   Rainfall_imp_mcmc : if keep_imp=true, matrix with the chain of imputed Rainfall for archaeological data. dim(theta_mcmc)=(n_iter,n_obs_arc)
   loglik_mcmc: if keep_ll=true, matrix with log likelihood of PO and HM modules
   */
  
  
  // Defining Auxiliar variables
  unsigned int iter_i=0;
  
  // unsigned int i=0;
  unsigned int j=0;
  unsigned int theta_j=0;
  // int aux_int=0;
  // double aux_double=0;
  arma::vec aux_vec;
  // arma::colvec aux_colvec;
  arma::uvec aux_uvec;
  arma::mat aux_mat;
  bool aux_bool=false;
  
  // number of rows in data
  const unsigned int n_obs_arc = data_arc.n_rows;
  const unsigned int n_obs_mod = data_mod.n_rows;
  
  // number of parameters
  unsigned int n_par = theta.n_rows;
  
  // MCMC warm-up and thinning //
  // vector with the iterations that will be kept and returned by the function
  arma::uvec iter_out = arma::regspace<arma::uvec>(n_warmup+1,n_thin,n_iter);
  arma::uvec iter_i_save;
  
  // OUTPUT
  arma::mat theta_mcmc;
  arma::mat ManureLevel_imp_mcmc;
  arma::mat Rainfall_imp_mcmc;
  arma::cube loglik_mcmc;
  
  theta_mcmc = arma::zeros<arma::mat>(iter_out.n_rows,n_par);
  if(keep_imp){
    ManureLevel_imp_mcmc = arma::zeros<arma::mat>(iter_out.n_rows,n_obs_arc);
    Rainfall_imp_mcmc = arma::zeros<arma::mat>(iter_out.n_rows,n_obs_arc);
  } else {
    ManureLevel_imp_mcmc = arma::zeros<arma::mat>(1,1);
    Rainfall_imp_mcmc = arma::zeros<arma::mat>(1,1);
  }
  if(keep_ll){
    loglik_mcmc = arma::zeros<arma::cube>(iter_out.n_rows,n_obs_arc+n_obs_mod,2);
  } else {
    loglik_mcmc = arma::zeros<arma::cube>(1,1,1);
  }
  
  // Sites
  arma::vec sites_arc = arma::unique(data_arc.col(1));
  arma::vec sites_mod = arma::unique(data_mod.col(0));
  unsigned int n_sites_arc = sites_arc.n_rows;
  unsigned int n_sites_mod = sites_mod.n_rows;
  
  
  if( n_par != (11+2*n_sites_arc+n_sites_mod) ) {
    throw std::range_error("theta is not consistent with the data!");
  }
  if( theta_min_max.n_rows != n_par ) {
    throw std::range_error("theta_min_max is not consistent with the data!");
  }
  if( theta_prop_int.n_rows != n_par ) {
    throw std::range_error("theta_prop_int is not consistent with the data!");
  }
  
  // Sites design matrix
  arma::mat Site_PO_dummie = categ_to_dummie( data_arc.col(1), sites_arc );
  
  arma::mat Site_HM_dummie = join_vert( categ_to_dummie( data_arc.col(1), join_vert( sites_arc, sites_mod ) ),
                                        categ_to_dummie( data_mod.col(0), join_vert( sites_arc, sites_mod ) ) );
  
  // ManureLevel
  arma::vec ManureLevels;
  ManureLevels << 1 << 2 << 3 << arma::endr;
  
  // ManureLevel dummy matrix
  arma::mat ManureLevel_HM_dummie;
  ManureLevel_HM_dummie = join_vert( categ_to_dummie( ManureLevel_imp, ManureLevels ),
                                     categ_to_dummie( data_mod.col(1), ManureLevels ) );
  
  ManureLevel_HM_dummie.shed_col(0);
  
  // behaviour of power_w_PO depending on the type of imputation
  if( gibbs_hm & (power_w_HM!=1) ) {
    throw std::range_error("Gibbs sampling only implemented for power_w_HM=1");
  }
  
  if( (theta_prop_kernel!="norm") & (theta_prop_kernel!="unif") ){
    throw std::range_error("theta_prop_kernel not supported!");
  }
  // vectors with probabilities used for Metropolis-Hastings steps
  arma::colvec cur_probs = arma::zeros(4);
  arma::colvec prop_probs = arma::zeros(4);
  arma::colvec check_probs = arma::zeros(4);
  // auxiliar value for the proposal random walk
  double rho = 0;
  // Metropolis Hastings ratio
  double MH_ratio = 0;
  // Sampling order
  arma::ucolvec sampling_order;
  
  
  // Covariate and response matrices
  arma::colvec Y_PO;
  arma::colvec Y_PO_prop;
  arma::mat X_PO;
  arma::colvec Y_HM;
  arma::mat X_HM;
  arma::mat X_HM_prop;
  
  arma::colvec ind_v_HM;
  
  // Matrices for gibbs sampling
  arma::mat Y_HM_gibbs;
  arma::mat Sigma_Y_HM;
  arma::mat Sigma_Y_HM_inv;
  arma::colvec mu_prior_gibbs;
  arma::mat Sigma_prior_gibbs;
  
  arma::colvec mu_post_gibbs;
  arma::mat Sigma_post_gibbs;
  
  arma::colvec sample_post_gibbs;
  
  Y_PO = ManureLevel_imp;
  
  Y_HM = join_vert( data_arc.col(5),
                    data_mod.col(4) );
  
  X_HM = arma::ones<arma::mat>(n_obs_arc+n_obs_mod,1);
  X_HM.insert_cols( X_HM.n_cols, join_vert( Rainfall_imp,
                                            data_mod.col(2) ) ); // columns: Rainfall_imp
  X_HM.insert_cols( X_HM.n_cols, ManureLevel_HM_dummie );
  
  ind_v_HM = join_vert( data_arc.col(4),
                        data_mod.col(3) );
  
  Sigma_Y_HM = arma::eye<arma::mat>(Y_HM.n_rows,Y_HM.n_rows);
  Sigma_Y_HM.diag().fill( pow(theta( 8+2*n_sites_arc+n_sites_mod ),2) );
  aux_vec = arma::ones<arma::vec>(Y_HM.n_rows);
  aux_vec.rows(find(ind_v_HM)).fill( theta( 9+2*n_sites_arc+n_sites_mod ) );
  Sigma_Y_HM.diag() = aux_vec % Sigma_Y_HM.diag();
  Sigma_Y_HM_inv = arma::inv_sympd(Sigma_Y_HM);
  
  // proposed values of the chain
  arma::colvec theta_prop;
  arma::colvec ManureLevel_imp_prop;
  arma::colvec Rainfall_imp_prop;
  arma::mat ManureLevel_HM_dummie_prop;
  
  theta_prop = theta;
  ManureLevel_imp_prop = ManureLevel_imp;
  Rainfall_imp_prop = Rainfall_imp;
  ManureLevel_HM_dummie_prop = ManureLevel_HM_dummie;
  
  // MCMC process monitoring //
  // vector with the number of iteration when the progress percentage will be displayed
  arma::vec iter_progress = arma::linspace<arma::vec>(0, n_iter, 21);
  iter_progress = arma::round(iter_progress);
  iter_progress.shed_row(0);
  iter_progress(19) = n_iter;
  
  if(!PO_site_rnd_eff){
    theta.rows( arma::span(3,2+n_sites_arc) ).fill(0); // eta_PO
    theta( 3+n_sites_arc ) = 1.0 ; // sigma_eta_PO
  }
  
  if(!HM_site_rnd_eff){
    theta.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ).fill(0); // eta_HM
    theta( 10+2*n_sites_arc+n_sites_mod ) = 1.0 ; // sigma_eta_HM
  }
  
  // Set current probabilities //
  cur_probs = loglik_prior_PO_HM( prior_spec_PO,
                                  prior_spec_HM,
                                  
                                  Y_PO,
                                  data_arc.col(0),
                                  Site_PO_dummie,
                                  
                                  theta.rows( arma::span(0,1) ), // alpha_PO
                                  theta.row( 2 ), // gamma_PO
                                  theta.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                  theta( 3+n_sites_arc ), // sigma_eta_PO
                                  
                                  Y_HM,
                                  X_HM,
                                  Site_HM_dummie,
                                  ind_v_HM,
                                  
                                  theta.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                  theta.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                  theta( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                  theta( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                  theta( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                  
                                  power_w_PO,
                                  power_w_HM );
  if( arma::sum( cur_probs.rows(0,1) ) == -INFINITY ){ throw std::range_error("invalid model parameters, loglik=-INFINITY! (code 1)"); }
  
  //// Recording initial values in the chain ////
  iter_i = 0;
  if( arma::any( iter_out==iter_i ) ) {
    iter_i_save = arma::find( iter_out==iter_i );
    // Save theta
    theta_mcmc.row( iter_i_save(0) ) = theta.t();
    
    // Save imputed data
    if(keep_imp){
      ManureLevel_imp_mcmc.row( iter_i_save(0) ) = ManureLevel_imp.t();
      Rainfall_imp_mcmc.row( iter_i_save(0) ) = Rainfall_imp.t();
    }
    
    // Save log-likelihood
    if(keep_ll){
      loglik_mcmc.subcube(iter_i_save(0),0,0  , iter_i_save(0),n_obs_arc+n_obs_mod-1,1) = loglik_i_PO_HM( Y_PO,
                          data_arc.col(0),
                          Site_PO_dummie,
                          
                          theta.rows( arma::span(0,1) ), // alpha_PO
                          theta.row( 2 ), // gamma_PO
                          theta.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                          theta( 3+n_sites_arc ), // sigma_eta_PO
                          
                          Y_HM,
                          X_HM,
                          Site_HM_dummie,
                          ind_v_HM,
                          
                          theta.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                          theta.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                          theta( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                          theta( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                          theta( 10+2*n_sites_arc+n_sites_mod ) ); // sigma_eta_HM 
    }
  }
  
  Rcpp::Rcout << "Starting MCMC..." << std::endl << "progress:" << std::endl;
  for( iter_i=1; iter_i<n_iter; iter_i++ ) {
    if( arma::sum( cur_probs.rows(0,1) ) == -INFINITY ){ throw std::range_error("invalid model parameters, loglik=-INFINITY! (code 1)"); }
    
    if(verbose){
      Rcpp::Rcout << "iter_i=" << iter_i << std::endl;
    }
    
    // check that "cur_probs" is consistent
    if( check_mcmc ) {
      check_probs = loglik_prior_PO_HM( prior_spec_PO,
                                        prior_spec_HM,
                                        
                                        Y_PO,
                                        data_arc.col(0),
                                        Site_PO_dummie,
                                        
                                        theta.rows( arma::span(0,1) ), // alpha_PO
                                        theta.row( 2 ), // gamma_PO
                                        theta.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                        theta( 3+n_sites_arc ), // sigma_eta_PO
                                        
                                        Y_HM,
                                        X_HM,
                                        Site_HM_dummie,
                                        ind_v_HM,
                                        
                                        theta.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                        theta.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                        theta( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                        theta( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                        theta( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                        
                                        power_w_PO,
                                        power_w_HM );
      
      if( !approx_equal(cur_probs, check_probs, "absdiff", 0.00001 ) ) {
        throw std::range_error("There's a problem with the MCMC (code 1)");
      }
    }
    
    //////////
    // SAMPLING HM MODEL //
    //////////
    
    if( gibbs_hm ) {
      // Gibbs sampling for sampling beta_hm and eta_hm //
      
      // Sampling beta_hm
      Y_HM_gibbs = Y_HM - Site_HM_dummie * theta.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) );
      
      mu_prior_gibbs = arma::zeros<arma::colvec>(4); // Prior mean = zero
      Sigma_prior_gibbs = arma::eye<arma::mat>(4,4); // Prior covariance matrix = diagonal ones
      
      aux_mat = Sigma_prior_gibbs + X_HM.t() * Sigma_Y_HM_inv * X_HM;
      Sigma_post_gibbs = arma::inv_sympd(aux_mat);  // posterior covariance
      
      mu_post_gibbs = arma::inv_sympd(Sigma_prior_gibbs) * mu_prior_gibbs + X_HM.t() * ( Sigma_Y_HM_inv.diag() % Y_HM_gibbs);
      mu_post_gibbs = Sigma_post_gibbs * mu_post_gibbs;  // posterior mean
      
      sample_post_gibbs = arma::randn<arma::colvec>( 4 );
      sample_post_gibbs = mu_post_gibbs + arma::chol(Sigma_post_gibbs) * sample_post_gibbs ;
      
      theta.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ) = sample_post_gibbs;
      
      // Sampling eta_hm
      if( HM_site_rnd_eff ) {
        Y_HM_gibbs = Y_HM - X_HM * theta.rows( arma::span(4+n_sites_arc,7+n_sites_arc) );
        
        mu_prior_gibbs = arma::zeros<arma::colvec>(n_sites_arc+n_sites_mod); // Prior mean = zero
        Sigma_prior_gibbs = arma::eye<arma::mat>(n_sites_arc+n_sites_mod,n_sites_arc+n_sites_mod); // Prior covariance matrix = diagonal sigma_hm_eta^2
        Sigma_prior_gibbs.diag().fill( pow(theta( 10+2*n_sites_arc+n_sites_mod ),2) );
        
        aux_mat = Sigma_prior_gibbs + Site_HM_dummie.t() * (Sigma_Y_HM_inv * Site_HM_dummie);
        Sigma_post_gibbs = arma::inv_sympd(aux_mat);  // posterior covariance
        
        mu_post_gibbs = arma::inv_sympd(Sigma_prior_gibbs) * mu_prior_gibbs + Site_HM_dummie.t() * ( Sigma_Y_HM_inv.diag() % Y_HM_gibbs);
        mu_post_gibbs = Sigma_post_gibbs * mu_post_gibbs;  // posterior mean
        
        sample_post_gibbs = arma::randn<arma::colvec>( n_sites_arc+n_sites_mod );
        sample_post_gibbs = mu_post_gibbs + arma::chol(Sigma_post_gibbs) * sample_post_gibbs ;
        
        theta.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ) = sample_post_gibbs;
        
      }
      
      // Compute cur_probs with the updated values of theta
      cur_probs = loglik_prior_PO_HM( prior_spec_PO,
                                      prior_spec_HM,
                                      
                                      Y_PO,
                                      data_arc.col(0),
                                      Site_PO_dummie,
                                      
                                      theta.rows( arma::span(0,1) ), // alpha_PO
                                      theta.row( 2 ), // gamma_PO
                                      theta.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                      theta( 3+n_sites_arc ), // sigma_eta_PO
                                      
                                      Y_HM,
                                      X_HM,
                                      Site_HM_dummie,
                                      ind_v_HM,
                                      
                                      theta.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                      theta.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                      theta( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                      theta( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                      theta( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                      
                                      power_w_PO,
                                      power_w_HM );
      
      
    } else {
      // Metropolis Hastings for sampling beta_hm and eta_hm
      
      // Sampling beta_hm
      if(verbose){
        Rcpp::Rcout << "beta_hm, " << std::endl;
      }
      sampling_order = arma::regspace<arma::ucolvec>(4+n_sites_arc,7+n_sites_arc);
      sampling_order = arma::shuffle(sampling_order);
      for( j=0; j<sampling_order.n_rows; j++ ) {
        
        theta_j = sampling_order(j);
        
        // Reseting proposal values of theta
        theta_prop = theta;
        
        // Sampling a new value from the proposal
        if(theta_prop_kernel=="unif"){
          rho = R::runif(-theta_prop_int(theta_j)/2,theta_prop_int(theta_j)/2);
        } else if(theta_prop_kernel=="norm"){
          rho = R::rnorm(0,theta_prop_int(theta_j));
        }
        theta_prop(theta_j) = theta(theta_j) + rho;
        // bounce in min and max posible values
        theta_prop(theta_j) = bounce_limit( theta_prop(theta_j), theta_min_max(theta_j,0), theta_min_max(theta_j,1) );
        
        // MH step //
        prop_probs = loglik_prior_PO_HM( prior_spec_PO,
                                         prior_spec_HM,
                                         
                                         Y_PO,
                                         data_arc.col(0),
                                         Site_PO_dummie,
                                         
                                         theta_prop.rows( arma::span(0,1) ), // alpha_PO
                                         theta_prop.row( 2 ), // gamma_PO
                                         theta_prop.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                         theta_prop( 3+n_sites_arc ), // sigma_eta_PO
                                         
                                         Y_HM,
                                         X_HM,
                                         Site_HM_dummie,
                                         ind_v_HM,
                                         
                                         theta_prop.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                         theta_prop.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                         theta_prop( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                         theta_prop( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                         theta_prop( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                         
                                         power_w_PO,
                                         power_w_HM );
        
        // Metropolis-Hastings ratio
        MH_ratio = arma::sum( prop_probs.rows(2,3) ) - arma::sum( cur_probs.rows(2,3) );
        
        if( log(R::runif(0,1)) < MH_ratio ) {
          // Updating current values of the MCMC
          theta = theta_prop;
          cur_probs = prop_probs;
        }
      }
      
      // Sampling eta_hm
      if( HM_site_rnd_eff ) {
        if(verbose){
          Rcpp::Rcout << "eta_hm, " << std::endl;
        }
        sampling_order = arma::regspace<arma::ucolvec>(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod);
        sampling_order = arma::shuffle(sampling_order);
        for( j=0; j<sampling_order.n_rows;j++ ) {
          
          theta_j = sampling_order(j);
          
          // Reseting proposal values of theta
          theta_prop = theta;
          
          // Sampling a new value from the proposal
          if(theta_prop_kernel=="unif"){
            rho = R::runif(-theta_prop_int(theta_j)/2,theta_prop_int(theta_j)/2);
          } else if(theta_prop_kernel=="norm"){
            rho = R::rnorm(0,theta_prop_int(theta_j));
          }
          theta_prop(theta_j) = theta(theta_j) + rho;
          // bounce in min and max posible values
          theta_prop(theta_j) = bounce_limit( theta_prop(theta_j), theta_min_max(theta_j,0), theta_min_max(theta_j,1) );
          
          // MH step //
          prop_probs = loglik_prior_PO_HM( prior_spec_PO,
                                           prior_spec_HM,
                                           
                                           Y_PO,
                                           data_arc.col(0),
                                           Site_PO_dummie,
                                           
                                           theta_prop.rows( arma::span(0,1) ), // alpha_PO
                                           theta_prop.row( 2 ), // gamma_PO
                                           theta_prop.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                           theta_prop( 3+n_sites_arc ), // sigma_eta_PO
                                           
                                           Y_HM,
                                           X_HM,
                                           Site_HM_dummie,
                                           ind_v_HM,
                                           
                                           theta_prop.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                           theta_prop.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                           theta_prop( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                           theta_prop( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                           theta_prop( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                           
                                           power_w_PO,
                                           power_w_HM );
          
          // Metropolis-Hastings ratio
          MH_ratio = arma::sum( prop_probs.rows(2,3) ) - arma::sum( cur_probs.rows(2,3) );
          
          if( log(R::runif(0,1)) < MH_ratio ) {
            // Updating current values of the MCMC
            theta = theta_prop;
            cur_probs = prop_probs;
          }
        }
      }
    }
    
    // Sampling sigma_hm, v and sigma_eta_hm
    if(verbose){
      Rcpp::Rcout << "sigma_hm, v_HM, sigma_eta_HM " << std::endl;
    }
    sampling_order = arma::zeros<arma::ucolvec>(3);
    sampling_order(0) = 8+2*n_sites_arc+n_sites_mod; // sigma_HM
    sampling_order(1) = 9+2*n_sites_arc+n_sites_mod; // v_HM
    sampling_order(2) = 10+2*n_sites_arc+n_sites_mod; // sigma_eta_HM
    sampling_order = arma::shuffle(sampling_order);
    for( j=0; j<3; j++ ) {
      
      theta_j = sampling_order(j);
      
      // Reseting proposal values of theta
      theta_prop = theta;
      
      // Sampling a new value from the proposal
      if(theta_prop_kernel=="unif"){
        rho = R::runif(-theta_prop_int(theta_j)/2,theta_prop_int(theta_j)/2);
      } else if(theta_prop_kernel=="norm"){
        rho = R::rnorm(0,theta_prop_int(theta_j));
      }
      theta_prop(theta_j) = theta(theta_j) + rho;
      // bounce in min and max posible values
      theta_prop(theta_j) = bounce_limit( theta_prop(theta_j), theta_min_max(theta_j,0), theta_min_max(theta_j,1) );
      
      // MH step //
      prop_probs = loglik_prior_PO_HM( prior_spec_PO,
                                       prior_spec_HM,
                                       
                                       Y_PO,
                                       data_arc.col(0),
                                       Site_PO_dummie,
                                       
                                       theta_prop.rows( arma::span(0,1) ), // alpha_PO
                                       theta_prop.row( 2 ), // gamma_PO
                                       theta_prop.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                       theta_prop( 3+n_sites_arc ), // sigma_eta_PO
                                       
                                       Y_HM,
                                       X_HM,
                                       Site_HM_dummie,
                                       ind_v_HM,
                                       
                                       theta_prop.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                       theta_prop.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                       theta_prop( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                       theta_prop( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                       theta_prop( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                       
                                       power_w_PO,
                                       power_w_HM );
      
      // Metropolis-Hastings ratio
      MH_ratio = arma::sum( prop_probs.rows(2,3) ) - arma::sum( cur_probs.rows(2,3) );
      
      if( log(R::runif(0,1)) < MH_ratio ) {
        // Updating current values of the MCMC
        theta = theta_prop;
        cur_probs = prop_probs;
      }
    }
    
    
    //////////
    // SAMPLING PO MODEL //
    //////////
    
      
      // Sampling alpha_po
      if(verbose){
        Rcpp::Rcout << "alpha_po, " << std::endl;
      }
      sampling_order = arma::regspace<arma::ucolvec>(0,1);
      sampling_order = arma::shuffle(sampling_order);
      for( j=0; j<sampling_order.n_rows;j++ ) {
        
        theta_j = sampling_order(j);
        
        // Reseting proposal values of theta
        theta_prop = theta;
        
        // Sampling a new value from the proposal
        if(theta_prop_kernel=="unif"){
          rho = R::runif(-theta_prop_int(theta_j)/2,theta_prop_int(theta_j)/2);
        } else if(theta_prop_kernel=="norm"){
          rho = R::rnorm(0,theta_prop_int(theta_j));
        }
        theta_prop(theta_j) = theta(theta_j) + rho;
        // bounce in min and max posible values
        theta_prop(theta_j) = bounce_limit( theta_prop(theta_j), theta_min_max(theta_j,0), theta_min_max(theta_j,1) );
        
        // MH step //
        prop_probs = loglik_prior_PO_HM( prior_spec_PO,
                                         prior_spec_HM,
                                         
                                         Y_PO,
                                         data_arc.col(0),
                                         Site_PO_dummie,
                                         
                                         theta_prop.rows( arma::span(0,1) ), // alpha_PO
                                         theta_prop.row( 2 ), // gamma_PO
                                         theta_prop.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                         theta_prop( 3+n_sites_arc ), // sigma_eta_PO
                                         
                                         Y_HM,
                                         X_HM,
                                         Site_HM_dummie,
                                         ind_v_HM,
                                         
                                         theta_prop.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                         theta_prop.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                         theta_prop( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                         theta_prop( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                         theta_prop( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                         
                                         power_w_PO,
                                         power_w_HM );
        
        // Metropolis-Hastings ratio
        MH_ratio = arma::sum( prop_probs.rows(0,1) ) - arma::sum( cur_probs.rows(0,1) );
        
        if( log(R::runif(0,1)) < MH_ratio ) {
          // Updating current values of the MCMC
          theta = theta_prop;
          cur_probs = prop_probs;
        }
      }
      
      // Sampling gamma_po
      if(verbose){
        Rcpp::Rcout << "gamma_po, " << std::endl;
      }
      for( j=0; j<1;j++ ) {
        
        theta_j = 2;
        
        // Reseting proposal values of theta
        theta_prop = theta;
        
        // Sampling a new value from the proposal
        if(theta_prop_kernel=="unif"){
          rho = R::runif(-theta_prop_int(theta_j)/2,theta_prop_int(theta_j)/2);
        } else if(theta_prop_kernel=="norm"){
          rho = R::rnorm(0,theta_prop_int(theta_j));
        }
        theta_prop(theta_j) = theta(theta_j) + rho;
        // bounce in min and max posible values
        theta_prop(theta_j) = bounce_limit( theta_prop(theta_j), theta_min_max(theta_j,0), theta_min_max(theta_j,1) );
        
        // MH step //
        
        prop_probs = loglik_prior_PO_HM( prior_spec_PO,
                                         prior_spec_HM,
                                         
                                         Y_PO,
                                         data_arc.col(0),
                                         Site_PO_dummie,
                                         
                                         theta_prop.rows( arma::span(0,1) ), // alpha_PO
                                         theta_prop.row( 2 ), // gamma_PO
                                         theta_prop.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                         theta_prop( 3+n_sites_arc ), // sigma_eta_PO
                                         
                                         Y_HM,
                                         X_HM,
                                         Site_HM_dummie,
                                         ind_v_HM,
                                         
                                         theta_prop.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                         theta_prop.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                         theta_prop( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                         theta_prop( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                         theta_prop( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                         
                                         power_w_PO,
                                         power_w_HM );
        
        // Metropolis-Hastings ratio
        MH_ratio = arma::sum( prop_probs.rows(0,1) ) - arma::sum( cur_probs.rows(0,1) );
        // if(verbose){
        //   Rcpp::Rcout << "\t prop=" << arma::sum( prop_probs.rows(0,1) )  << std::endl;
        //   Rcpp::Rcout << "\t cur=" << arma::sum( cur_probs.rows(0,1) ) << std::endl;
        //   Rcpp::Rcout << "\t MH_ratio=" << MH_ratio << std::endl;
        // }
        
        if( log(R::runif(0,1)) < MH_ratio ) {
          // Updating current values of the MCMC
          theta = theta_prop;
          cur_probs = prop_probs;
        }
      }
      
      if( PO_site_rnd_eff ) {
        // Sampling eta_po
        if(verbose){
          Rcpp::Rcout << "eta_po, " << std::endl;
        }
        sampling_order = arma::regspace<arma::ucolvec>(3,2+n_sites_arc);
        sampling_order = arma::shuffle(sampling_order);
        for( j=0; j<sampling_order.n_rows;j++ ) {
          
          theta_j = sampling_order(j);
          
          // Reseting proposal values of theta
          theta_prop = theta;
          
          // Sampling a new value from the proposal
          if(theta_prop_kernel=="unif"){
            rho = R::runif(-theta_prop_int(theta_j)/2,theta_prop_int(theta_j)/2);
          } else if(theta_prop_kernel=="norm"){
            rho = R::rnorm(0,theta_prop_int(theta_j));
          }
          theta_prop(theta_j) = theta(theta_j) + rho;
          // bounce in min and max posible values
          theta_prop(theta_j) = bounce_limit( theta_prop(theta_j), theta_min_max(theta_j,0), theta_min_max(theta_j,1) );
          
          // MH step //
          
          prop_probs = loglik_prior_PO_HM( prior_spec_PO,
                                           prior_spec_HM,
                                           
                                           Y_PO,
                                           data_arc.col(0),
                                           Site_PO_dummie,
                                           
                                           theta_prop.rows( arma::span(0,1) ), // alpha_PO
                                           theta_prop.row( 2 ), // gamma_PO
                                           theta_prop.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                           theta_prop( 3+n_sites_arc ), // sigma_eta_PO
                                           
                                           Y_HM,
                                           X_HM,
                                           Site_HM_dummie,
                                           ind_v_HM,
                                           
                                           theta_prop.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                           theta_prop.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                           theta_prop( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                           theta_prop( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                           theta_prop( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                           
                                           power_w_PO,
                                           power_w_HM );
          
          // Metropolis-Hastings ratio
          MH_ratio = arma::sum( prop_probs.rows(0,1) ) - arma::sum( cur_probs.rows(0,1) );
          
          if( log(R::runif(0,1)) < MH_ratio ) {
            // Updating current values of the MCMC
            theta = theta_prop;
            cur_probs = prop_probs;
          }
        }
        
        // Sampling sigma_eta_PO
        if(verbose){
          Rcpp::Rcout << "sigma_eta_PO, " << std::endl;
        }
        for( j=0; j<1; j++ ) {
          
          theta_j = 3+n_sites_arc;
          
          // Reseting proposal values of theta
          theta_prop = theta;
          
          // Sampling a new value from the proposal
          if(theta_prop_kernel=="unif"){
            rho = R::runif(-theta_prop_int(theta_j)/2,theta_prop_int(theta_j)/2);
          } else if(theta_prop_kernel=="norm"){
            rho = R::rnorm(0,theta_prop_int(theta_j));
          }
          theta_prop(theta_j) = theta(theta_j) + rho;
          // bounce in min and max posible values
          theta_prop(theta_j) = bounce_limit( theta_prop(theta_j), theta_min_max(theta_j,0), theta_min_max(theta_j,1) );
          
          
          // MH step //
          prop_probs = loglik_prior_PO_HM( prior_spec_PO,
                                           prior_spec_HM,
                                           
                                           Y_PO,
                                           data_arc.col(0),
                                           Site_PO_dummie,
                                           
                                           theta_prop.rows( arma::span(0,1) ), // alpha_PO
                                           theta_prop.row( 2 ), // gamma_PO
                                           theta_prop.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                           theta_prop( 3+n_sites_arc ), // sigma_eta_PO
                                           
                                           Y_HM,
                                           X_HM,
                                           Site_HM_dummie,
                                           ind_v_HM,
                                           
                                           theta_prop.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                           theta_prop.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                           theta_prop( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                           theta_prop( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                           theta_prop( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                           
                                           power_w_PO,
                                           power_w_HM );
          
          // Metropolis-Hastings ratio
          MH_ratio = arma::sum( prop_probs.rows(0,1) ) - arma::sum( cur_probs.rows(0,1) );
          
          if( log(R::runif(0,1)) < MH_ratio ) {
            // Updating current values of the MCMC
            theta = theta_prop;
            cur_probs = prop_probs;
          }
        }
      }
    
    
    // check that "cur_probs" is consistent
    if( check_mcmc ) {
      check_probs = loglik_prior_PO_HM( prior_spec_PO,
                                        prior_spec_HM,
                                        
                                        Y_PO,
                                        data_arc.col(0),
                                        Site_PO_dummie,
                                        
                                        theta.rows( arma::span(0,1) ), // alpha_PO
                                        theta.row( 2 ), // gamma_PO
                                        theta.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                        theta( 3+n_sites_arc ), // sigma_eta_PO
                                        
                                        Y_HM,
                                        X_HM,
                                        Site_HM_dummie,
                                        ind_v_HM,
                                        
                                        theta.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                        theta.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                        theta( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                        theta( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                        theta( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                        
                                        power_w_PO,
                                        power_w_HM );
      
      if( !approx_equal(cur_probs, check_probs, "absdiff", 0.00001 ) ) {
        throw std::range_error("There's a problem with the MCMC (code 2)");
      }
    }
    
    //////////
    // IMPUTING ManureLevels //
    //////////
    if(verbose){
      Rcpp::Rcout << "ManureLevels, " << std::endl;
    }
    
    // Sampling one by one //
    sampling_order = arma::regspace<arma::ucolvec>(1,n_obs_arc);
    sampling_order = arma::shuffle(sampling_order);
    
    for( j=0; j<n_obs_arc; j++ ){
      
      // Sampling new random value //
      ManureLevel_imp_prop = ManureLevel_imp;
      aux_vec = arma::shuffle(ManureLevels);
      ManureLevel_imp_prop( sampling_order(j)-1 ) = aux_vec(0);
      
      // Check that the proposal is different than the current value
      aux_bool = ManureLevel_imp_prop( sampling_order(j)-1 ) != ManureLevel_imp( sampling_order(j)-1 );
      // Check that the proposal does not make PO module fall into (quasi-complete) separability
      if(imp_playpen){
        aux_bool = aux_bool & !check_po_separability( ManureLevel_imp_prop, data_arc.col(0), 3 );
      }
      
      if( aux_bool ) {
        
        // Imputing sampled value //
        Y_PO_prop = ManureLevel_imp_prop;
        ManureLevel_HM_dummie_prop = join_vert( categ_to_dummie( ManureLevel_imp_prop, ManureLevels ),
                                                categ_to_dummie( data_mod.col(1), ManureLevels ) );
        ManureLevel_HM_dummie_prop.shed_col(0);
        
        X_HM_prop = arma::ones<arma::mat>(n_obs_arc+n_obs_mod,1);
        X_HM_prop.insert_cols( X_HM_prop.n_cols, join_vert( Rainfall_imp,
                                                            data_mod.col(2) ) ); // columns: Rainfall_imp
        X_HM_prop.insert_cols( X_HM_prop.n_cols, ManureLevel_HM_dummie_prop );
        
        // MH step //
        prop_probs = loglik_prior_PO_HM( prior_spec_PO,
                                         prior_spec_HM,
                                         
                                         Y_PO_prop,
                                         data_arc.col(0),
                                         Site_PO_dummie,
                                         
                                         theta.rows( arma::span(0,1) ), // alpha_PO
                                         theta.row( 2 ), // gamma_PO
                                         theta.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                         theta( 3+n_sites_arc ), // sigma_eta_PO
                                         
                                         Y_HM,
                                         X_HM_prop,
                                         Site_HM_dummie,
                                         ind_v_HM,
                                         
                                         theta.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                         theta.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                         theta( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                         theta( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                         theta( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                         
                                         power_w_PO,
                                         power_w_HM );
        
        // Metropolis-Hastings ratio
        MH_ratio = arma::sum( prop_probs ) - arma::sum( cur_probs );
        
        if( log(R::runif(0,1)) < MH_ratio ) {
          // Updating current values of the MCMC
          ManureLevel_imp = ManureLevel_imp_prop;
          Y_PO = Y_PO_prop;
          ManureLevel_HM_dummie = ManureLevel_HM_dummie_prop;
          X_HM = X_HM_prop;
          cur_probs = prop_probs;
        }
        
      }
      
      // check that "cur_probs" is consistent
      if( check_mcmc ) {
        check_probs = loglik_prior_PO_HM( prior_spec_PO,
                                          prior_spec_HM,
                                          
                                          Y_PO,
                                          data_arc.col(0),
                                          Site_PO_dummie,
                                          
                                          theta.rows( arma::span(0,1) ), // alpha_PO
                                          theta.row( 2 ), // gamma_PO
                                          theta.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                          theta( 3+n_sites_arc ), // sigma_eta_PO
                                          
                                          Y_HM,
                                          X_HM,
                                          Site_HM_dummie,
                                          ind_v_HM,
                                          
                                          theta.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                          theta.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                          theta( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                          theta( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                          theta( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                          
                                          power_w_PO,
                                          power_w_HM );
        
        if( !approx_equal(cur_probs, check_probs, "absdiff", 0.00001 ) ) {
          throw std::range_error("There's a problem with the MCMC (code 3)");
        }
      }
      
    }
    
    
    //////////
    // IMPUTING Rainfall //
    //////////
    
    if(verbose){
      Rcpp::Rcout << "Rainfall, " << std::endl;
    }
    aux_vec = data_arc.col(3) - data_arc.col(2);
    aux_vec = aux_vec % arma::randu(n_obs_arc,1);
    aux_vec += data_arc.col(2);
    Rainfall_imp_prop = aux_vec;
    
    X_HM_prop = X_HM;
    X_HM_prop.col(1) = join_vert( Rainfall_imp_prop, data_mod.col(2) );
    
    prop_probs = loglik_prior_PO_HM( prior_spec_PO,
                                     prior_spec_HM,
                                     
                                     Y_PO,
                                     data_arc.col(0),
                                     Site_PO_dummie,
                                     
                                     theta.rows( arma::span(0,1) ), // alpha_PO
                                     theta.row( 2 ), // gamma_PO
                                     theta.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                                     theta( 3+n_sites_arc ), // sigma_eta_PO
                                     
                                     Y_HM,
                                     X_HM_prop,
                                     Site_HM_dummie,
                                     ind_v_HM,
                                     
                                     theta.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                                     theta.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                                     theta( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                                     theta( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                                     theta( 10+2*n_sites_arc+n_sites_mod ), // sigma_eta_HM
                                     
                                     power_w_PO,
                                     power_w_HM );
    
    // Metropolis-Hastings ratio
    MH_ratio = arma::sum( prop_probs ) - arma::sum( cur_probs );
    
    if( log(R::runif(0,1)) < MH_ratio ) {
      // Updating current values of the MCMC
      Rainfall_imp = Rainfall_imp_prop;
      X_HM = X_HM_prop;
      cur_probs = prop_probs;
    }
    
    //// Recording updated values in the chain ////
    if( arma::any( iter_out==iter_i ) ) {
      iter_i_save = arma::find( iter_out==iter_i );
      // Save theta
      theta_mcmc.row( iter_i_save(0) ) = theta.t();
      
      // Save imputed data
      if(keep_imp){
        ManureLevel_imp_mcmc.row( iter_i_save(0) ) = ManureLevel_imp.t();
        Rainfall_imp_mcmc.row( iter_i_save(0) ) = Rainfall_imp.t();
      }
      
      //// Compute and save log-likelihoods with current parameter values ////
      if(keep_ll){
        loglik_mcmc.subcube(iter_i_save(0),0,0  , iter_i_save(0),n_obs_arc+n_obs_mod-1,1) = loglik_i_PO_HM( Y_PO,
                            data_arc.col(0),
                            Site_PO_dummie,
                            
                            theta.rows( arma::span(0,1) ), // alpha_PO
                            theta.row( 2 ), // gamma_PO
                            theta.rows( arma::span(3,2+n_sites_arc) ), // eta_PO
                            theta( 3+n_sites_arc ), // sigma_eta_PO
                            
                            Y_HM,
                            X_HM,
                            Site_HM_dummie,
                            ind_v_HM,
                            
                            theta.rows( arma::span(4+n_sites_arc,7+n_sites_arc) ), // beta_HM
                            theta.rows( arma::span(8+n_sites_arc,7+2*n_sites_arc+n_sites_mod) ), // eta_HM
                            theta( 8+2*n_sites_arc+n_sites_mod ), // sigma_HM
                            theta( 9+2*n_sites_arc+n_sites_mod ), // v_HM
                            theta( 10+2*n_sites_arc+n_sites_mod ) ); // sigma_eta_HM 
      }
    }
    
    // MCMC progress monitoring //
    if ( any(iter_progress==(iter_i+1) ) ) {
      Rcpp::Rcout << round(100*(iter_i+1)/n_iter) << "% ";
    }
  }
  
  Rcpp::Rcout << std::endl << "...MCMC Completed!" << std::endl ;
  
  return Rcpp::List::create( Rcpp::Named("theta_mcmc") = theta_mcmc,
                             Rcpp::Named("ManureLevel_imp_mcmc") = ManureLevel_imp_mcmc,
                             Rcpp::Named("Rainfall_imp_mcmc") = Rainfall_imp_mcmc,
                             Rcpp::Named("loglik_mcmc") = loglik_mcmc );
  
}
