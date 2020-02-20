#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "aistats2020smi_shared.h"

// [[Rcpp::interfaces(r, cpp)]]

arma::colvec loglik_prior_PO( const std::string prior_spec,
                              
                              const arma::colvec Y,
                              const arma::mat X,
                              const arma::mat X_eta,
                              
                              const arma::colvec alpha,
                              const arma::colvec beta,
                              const arma::colvec eta,
                              const double sigma_eta,
                              
                              const double power_w=1 ) {
  
  // output: probs
  arma::vec probs = arma::zeros<arma::vec>(2);
  
  // "probs" contents:
  // 0: loglik
  // 1: logprior
  
  probs[0] = power_w * loglik_PO_cpp( Y,
                                      join_horiz(X,X_eta),
                                      alpha,
                                      join_vert(beta,eta) );
  
  probs[1] = logprior_PO_cpp( alpha,
                              beta,
                              sigma_eta,
                              eta,
                              prior_spec );
  
  return probs;
  
}



arma::mat loglik_i_PO( const arma::colvec Y,
                       const arma::mat X,
                       const arma::mat X_eta,
                       
                       const arma::colvec alpha,
                       const arma::colvec beta,
                       const arma::colvec eta,
                       const double sigma_eta ) {
  
  const unsigned int n_obs = X.n_rows;
  
  // output: loglik_i
  arma::mat loglik_i = arma::zeros<arma::mat>(n_obs,1);
  
  loglik_i.col(0) = loglik_PO_i_cpp( Y,
               join_horiz(X,X_eta),
               alpha,
               join_vert(beta,eta) );
  
  return loglik_i;
  
}



// [[Rcpp::export]]
Rcpp::List mcmc_PO( const arma::colvec Y,
                    const arma::mat X,
                    arma::mat X_eta,
                    
                    const unsigned int K,
                    
                    const double power_w,
                    
                    const std::string prior_spec,
                    
                    const bool rnd_eff,
                    
                    const unsigned int n_iter,
                    
                    arma::colvec theta,
                    const arma::mat theta_min_max,
                    arma::colvec theta_prop_int,
                    const std::string theta_prop_kernel,
                    
                    const bool keep_ll=true,
                    
                    const bool check_mcmc=false,
                    const bool verbose=false,
                    const bool quiet=false ) {
  /*
   INPUT:
   Y: vector with ordinal responses
   X: matrix with data for "fixed" effects (i.e. effect prior determined by prior_spec)
   X_eta: matrix with data for "random" effects (i.e. Gaussian prior with variance sigma_eta^2)
   K: integer>0 with the number of ordinal values to consider in Y
   */
  
  /*
   OUTPUT:
   theta_mcmc : matrix with the chain of parameters. dim(theta_mcmc)=(n_iter,n_par)
   */
  
  
  // Defining Auxiliar variables
  unsigned int iter_i=0;
  
  // unsigned int i=0;
  unsigned int j=0;
  unsigned int theta_j=0;
  // int aux_int=0;
  // double aux_double=0;
  arma::vec aux_vec;
  
  // number of rows in data
  const unsigned int n_obs = Y.n_rows;
  
  if(!rnd_eff){
    X_eta = arma::zeros<arma::mat>(n_obs,1);
  }
  
  if( (n_obs != X.n_rows) | (n_obs != X_eta.n_rows) ) {
    throw std::range_error("number of rows in response Y and covariates X,X_eta are not consistent!");
  }
  
  // number of parameters
  unsigned int n_par = theta.n_rows;
  unsigned int n_par_fix = X.n_cols;
  unsigned int n_par_rnd = X_eta.n_cols;
  
  if( n_par != ((K-1)+n_par_fix+n_par_rnd+1) ) {
    throw std::range_error("theta is not consistent with the data!");
  }
  
  
  // OUTPUT
  arma::mat theta_mcmc;
  arma::mat loglik_mcmc;
  
  theta_mcmc = arma::zeros<arma::mat>(n_iter,n_par);
  if(keep_ll){
    loglik_mcmc = arma::zeros<arma::mat>(n_iter,n_obs);
  } else {
    loglik_mcmc = arma::zeros<arma::mat>(1,1);
  }
  
  // METROPOLIS-HASTINGS //
  if( (theta_prop_kernel!="norm") & (theta_prop_kernel!="unif") ){
    throw std::range_error("theta_prop_kernel not supported!");
  }
  // vectors with probabilities used for Metropolis-Hastings steps
  arma::colvec cur_probs = arma::zeros(2);
  arma::colvec prop_probs = arma::zeros(2);
  arma::colvec check_probs = arma::zeros(2);
  // auxiliar value for the proposal random walk
  double rho = 0;
  // Metropolis Hastings ratio
  double MH_ratio = 0;
  // Sampling order
  arma::ucolvec sampling_order;
  
  // proposed values of the chain
  arma::colvec theta_prop;
  theta_prop = theta;
  
  // MCMC process monitoring //
  // vector with the number of iteration when the progress percentage will be displayed
  arma::vec iter_progress = arma::linspace<arma::vec>(0, n_iter, 21);
  iter_progress = arma::round(iter_progress);
  iter_progress.shed_row(0);
  iter_progress(19) = n_iter;
  
  // initialize random effects at zero
  // if(rnd_eff){
  //   theta.rows( arma::span(K-1+n_par_fix,K-1+n_par_fix+n_par_rnd-1) ).fill(0); // eta_PO
  // }
  
  // Set current probabilities //
  cur_probs = loglik_prior_PO( prior_spec,
                               
                               Y,
                               X,
                               X_eta,
                               
                               theta.rows( arma::span( 0, (K-1)-1 ) ), // alpha
                               theta.rows( arma::span( K-1, K-1+n_par_fix-1 ) ), // beta
                               theta.rows( arma::span( K-1+n_par_fix, K-1+n_par_fix+n_par_rnd-1 ) ), // eta
                               theta( K-1+n_par_fix+n_par_rnd ), // sigma_eta
                               
                               power_w );
  
  if( arma::sum( cur_probs ) == -INFINITY ){ throw std::range_error("invalid model parameters, loglik=-INFINITY! (code 1)"); }
  
  //// Recording initial values in the chain ////
  theta_mcmc.row(0) = theta.t();
  if(keep_ll){
    loglik_mcmc.row(0) = loglik_i_PO( Y,
                    X,
                    X_eta,
                    
                    theta.rows( arma::span( 0, (K-1)-1 ) ), // alpha
                    theta.rows( arma::span( K-1, K-1+n_par_fix-1 ) ), // beta
                    theta.rows( arma::span( K-1+n_par_fix, K-1+n_par_fix+n_par_rnd-1 ) ), // eta
                    theta( K-1+n_par_fix+n_par_rnd ) ).t(); // sigma_eta
  }
  
  if(!quiet){
    Rcpp::Rcout << "Starting MCMC..." << std::endl << "progress:" << std::endl;
  }
  
  for( iter_i=1; iter_i<n_iter; iter_i++ ) {
    if( arma::sum( cur_probs ) == -INFINITY ){ throw std::range_error("invalid model parameters, loglik=-INFINITY! (code 1)"); }
    
    if( verbose & (!quiet) ){
      Rcpp::Rcout << "iter_i=" << iter_i << std::endl;
    }
    
    // check that "cur_probs" is consistent
    if( check_mcmc ) {
      check_probs = loglik_prior_PO( prior_spec,
                                     
                                     Y,
                                     X,
                                     X_eta,
                                     
                                     theta.rows( arma::span( 0, (K-1)-1 ) ), // alpha
                                     theta.rows( arma::span( K-1, K-1+n_par_fix-1 ) ), // beta
                                     theta.rows( arma::span( K-1+n_par_fix, K-1+n_par_fix+n_par_rnd-1 ) ), // eta
                                     theta( K-1+n_par_fix+n_par_rnd ), // sigma_eta
                                     
                                     power_w );
      
      if( !approx_equal(cur_probs, check_probs, "absdiff", 0.00001 ) ) {
        throw std::range_error("There's a problem with the MCMC (code 1)");
      }
    }
    
    //////////
    // SAMPLING PO MODEL //
    //////////
    
    // Sampling alpha
    if( verbose & (!quiet) ){
      Rcpp::Rcout << "alpha, " << std::endl;
    }
    sampling_order = arma::regspace<arma::ucolvec>(0,K-1);
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
      prop_probs = loglik_prior_PO( prior_spec,
                                    
                                    Y,
                                    X,
                                    X_eta,
                                    
                                    theta_prop.rows( arma::span( 0, (K-1)-1 ) ), // alpha
                                    theta_prop.rows( arma::span( K-1, K-1+n_par_fix-1 ) ), // beta
                                    theta_prop.rows( arma::span( K-1+n_par_fix, K-1+n_par_fix+n_par_rnd-1 ) ), // eta
                                    theta_prop( K-1+n_par_fix+n_par_rnd ), // sigma_eta
                                    
                                    power_w );
      
      // Metropolis-Hastings ratio
      MH_ratio = arma::sum( prop_probs ) - arma::sum( cur_probs );
      
      if( log(R::runif(0,1)) < MH_ratio ) {
        // Updating current values of the MCMC
        theta = theta_prop;
        cur_probs = prop_probs;
      }
    }
    
    
    // Sampling beta ("fixed" effects)
    if( verbose & (!quiet)){
      Rcpp::Rcout << "alpha, " << std::endl;
    }
    sampling_order = arma::regspace<arma::ucolvec>( K, K+n_par_fix-1 );
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
      prop_probs = loglik_prior_PO( prior_spec,
                                    
                                    Y,
                                    X,
                                    X_eta,
                                    
                                    theta_prop.rows( arma::span( 0, (K-1)-1 ) ), // alpha
                                    theta_prop.rows( arma::span( K-1, K-1+n_par_fix-1 ) ), // beta
                                    theta_prop.rows( arma::span( K-1+n_par_fix, K-1+n_par_fix+n_par_rnd-1 ) ), // eta
                                    theta_prop( K-1+n_par_fix+n_par_rnd ), // sigma_eta
                                    
                                    power_w );
      
      // Metropolis-Hastings ratio
      MH_ratio = arma::sum( prop_probs ) - arma::sum( cur_probs );
      
      if( log(R::runif(0,1)) < MH_ratio ) {
        // Updating current values of the MCMC
        theta = theta_prop;
        cur_probs = prop_probs;
      }
    }
    
    if( rnd_eff ) {
      // Sampling eta
      if(verbose & (!quiet)){
        Rcpp::Rcout << "eta, " << std::endl;
      }
      sampling_order = arma::regspace<arma::ucolvec>( K+n_par_fix, K+n_par_fix+n_par_rnd-1 );
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
        prop_probs = loglik_prior_PO( prior_spec,
                                      
                                      Y,
                                      X,
                                      X_eta,
                                      
                                      theta_prop.rows( arma::span( 0, (K-1)-1 ) ), // alpha
                                      theta_prop.rows( arma::span( K-1, K-1+n_par_fix-1 ) ), // beta
                                      theta_prop.rows( arma::span( K-1+n_par_fix, K-1+n_par_fix+n_par_rnd-1 ) ), // eta
                                      theta_prop( K-1+n_par_fix+n_par_rnd ), // sigma_eta
                                      
                                      power_w );
        
        // Metropolis-Hastings ratio
        MH_ratio = arma::sum( prop_probs ) - arma::sum( cur_probs );
        
        if( log(R::runif(0,1)) < MH_ratio ) {
          // Updating current values of the MCMC
          theta = theta_prop;
          cur_probs = prop_probs;
        }
      }
      
      // Sampling sigma_eta_PO
      if(verbose & (!quiet)){
        Rcpp::Rcout << "sigma_eta_PO, " << std::endl;
      }
      for( j=0; j<1; j++ ) {
        
        theta_j = K-1+n_par_fix+n_par_rnd;
        
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
        prop_probs = loglik_prior_PO( prior_spec,
                                      
                                      Y,
                                      X,
                                      X_eta,
                                      
                                      theta_prop.rows( arma::span( 0, (K-1)-1 ) ), // alpha
                                      theta_prop.rows( arma::span( K-1, K-1+n_par_fix-1 ) ), // beta
                                      theta_prop.rows( arma::span( K-1+n_par_fix, K-1+n_par_fix+n_par_rnd-1 ) ), // eta
                                      theta_prop( K-1+n_par_fix+n_par_rnd ), // sigma_eta
                                      
                                      power_w );
        
        // Metropolis-Hastings ratio
        MH_ratio = arma::sum( prop_probs ) - arma::sum( cur_probs );
        
        if( log(R::runif(0,1)) < MH_ratio ) {
          // Updating current values of the MCMC
          theta = theta_prop;
          cur_probs = prop_probs;
        }
      }
    }
    
    
    // check that "cur_probs" is consistent
    if( check_mcmc ) {
      check_probs = loglik_prior_PO( prior_spec,
                                     
                                     Y,
                                     X,
                                     X_eta,
                                     
                                     theta.rows( arma::span( 0, (K-1)-1 ) ), // alpha
                                     theta.rows( arma::span( K-1, K-1+n_par_fix-1 ) ), // beta
                                     theta.rows( arma::span( K-1+n_par_fix, K-1+n_par_fix+n_par_rnd-1 ) ), // eta
                                     theta( K-1+n_par_fix+n_par_rnd ), // sigma_eta
                                     
                                     power_w );
      
      if( !approx_equal(cur_probs, check_probs, "absdiff", 0.00001 ) ) {
        throw std::range_error("There's a problem with the MCMC (code 2)");
      }
    }
    
    
    //// Compute log likelihood with current parameter values ////
    if(keep_ll){
      loglik_mcmc.row(iter_i) = loglik_i_PO( Y,
                      X,
                      X_eta,
                      
                      theta.rows( arma::span( 0, (K-1)-1 ) ), // alpha
                      theta.rows( arma::span( K-1, K-1+n_par_fix-1 ) ), // beta
                      theta.rows( arma::span( K-1+n_par_fix, K-1+n_par_fix+n_par_rnd-1 ) ), // eta
                      theta( K-1+n_par_fix+n_par_rnd ) ).t(); // sigma_eta
    }
    
    //// Recording updated values in the chain ////
    theta_mcmc.row(iter_i) = theta.t();
    
    // MCMC progress monitoring //
    if ( any(iter_progress==(iter_i+1) ) & (!quiet)) {
      Rcpp::Rcout << round(100*(iter_i+1)/n_iter) << "% ";
    }
    
  }
  
  if(!quiet){
    Rcpp::Rcout << std::endl << "...MCMC Completed!" << std::endl ;
  }
  
  return Rcpp::List::create( Rcpp::Named("theta_mcmc") = theta_mcmc,
                             Rcpp::Named("loglik_mcmc") = loglik_mcmc );
  
}
