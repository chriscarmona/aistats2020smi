#include <math.h>
#include <Rmath.h>
#include <algorithm>    // std::max
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "aistats2020smi_shared.h"

// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
arma::colvec loglik_PO_i_cpp( const arma::colvec Y,
                              const arma::mat X,
                              const arma::colvec alpha,
                              const arma::colvec beta ) {
  // const bool pointwise=false


  // Proportional Odds model
  // Y must be integer valued, with K categories: 1,2,...,k_Y
  // P(Y<=y_i)=plogis(alpha_i-X%*%beta); for y_i < max(levels(Y))
  // alpha: vector with the k_Y-1 intercepts
  // beta: (n.q)x(1) vector with the slopes

  unsigned int n_obs = X.n_rows;
  unsigned int p_vars = X.n_cols;

  arma::colvec loglik = arma::zeros<arma::colvec>(n_obs);

  // check that alpha's are ordered
  int n_alpha;
  n_alpha = alpha.size();
  for( int i=1;i<n_alpha;i++ ) {
    if( alpha[i-1]>alpha[i] ) {
      loglik.fill(-INFINITY);
      return loglik;
    }
  }

  // Number of categories in Y
  unsigned int k_Y;
  k_Y = alpha.n_rows+1;

  if( Y.n_rows!=n_obs ) {
    throw std::range_error("There's a problem with the data passed to loglik_PO (code 1)");
  }
  if( beta.n_rows!=p_vars ) {
    throw std::range_error("There's a problem with the data passed to loglik_PO (code 2)");
  }

  // linear predictor
  arma::colvec lin_pred = X * beta;

  for (unsigned int i=0;i<n_obs;i++){
    if( Y(i)==1 ) {
      loglik(i) = R::plogis( alpha(Y(i)-1)-lin_pred(i),0,1,1,1 );
    } else if( Y(i)<k_Y ) {
      loglik(i) = std::log( R::plogis( alpha(Y(i)-1)-lin_pred(i),0,1,1,0 ) - R::plogis( alpha(Y(i)-2)-lin_pred(i),0,1,1,0 ) );
    } else if( Y(i)==k_Y ) {
      loglik(i) = R::plogis( alpha(Y(i)-2)-lin_pred(i),0,1,0,1 );
    } else {
      Rcpp::Rcout << "Y(i)=" << Y(i) << ", k_Y=" << k_Y << ", i=" << i << std::endl;
      throw std::range_error("There's a problem with the data passed to loglik_PO (code 3)");
    }
  }

  // if(pointwise){
  //   return loglik;
  // }

  return loglik;

}

// [[Rcpp::export]]
double loglik_PO_cpp( const arma::colvec Y,
                      const arma::mat X,
                      const arma::colvec alpha,
                      const arma::colvec beta ) {

  unsigned int n_obs = X.n_rows;
  arma::colvec loglik = arma::zeros<arma::colvec>(n_obs);

  loglik = loglik_PO_i_cpp( Y,
                            X,
                            alpha,
                            beta );
  return sum( loglik );

}

// [[Rcpp::export]]
double logprior_PO_cpp( const arma::colvec alpha,
                        const arma::colvec beta,
                        const double sigma_eta,
                        const arma::colvec eta,
                        const std::string prior_spec ) {

  int n_alpha;
  n_alpha = alpha.size();

  int n_beta;
  n_beta = beta.size();

  int n_eta;
  n_eta = eta.size();

  double logprior;
  logprior = 0;

  // check that alpha's are ordered
  for( int i=1;i<n_alpha;i++ ) {
    if( alpha[i-1]>alpha[i] ) {
      return -INFINITY;
    }
  }

  if (prior_spec=="flat") {
    // alpha's
    // constant=1
    // Johnson & Albert (1999) section 4.2.1
    logprior += 0;
    // beta
    logprior += 0;
  } else if (prior_spec=="proper") {

    // Effects: normal prior
    // beta's ~ normal(0,sigma_beta)
    for(int i=0; i<n_beta; i++ ){
      logprior += R::dnorm( beta[i], 0, 5, 1 );
    }

    if( true ) {
      // intercept: normal prior
      for(int i=0;i<n_alpha;i++) {
        logprior += R::dnorm( alpha[i], 0, 5, 1 );
      }
    } else {
      // mean intercept: normal
      // intercept distances: exponential

      double mean_alpha;
      mean_alpha=0;
      for(int i=0;i<n_alpha;i++) {
        mean_alpha += alpha[i];
      }

      double sigma_alpha;
      sigma_alpha=5;

      logprior += R::dnorm(mean_alpha,0,sigma_alpha,1);

      // Distances between alpha's tied by an exponential distribution with mean "lambda_alpha"
      double lambda_alpha;
      lambda_alpha = 1.0/5.0;

      for(int i=1; i<n_alpha; i++ ){
        logprior += R::dexp( alpha[i]-alpha[i-1], lambda_alpha, 1 );
      }
      // curve(dexp(x,lambda_alpha),0,30)
    }

  } else {
    throw std::range_error("There's a problem with prior_spec in PO module");
  }

  // Random effects
  //if(!is.null(sigma_eta)&!is.null(eta)) {
  // sigma_eta: variance of random effects
  if( sigma_eta<0 ){
    return -INFINITY;
  } else {
    // Jeffrey's prior
    // logprior += -std::log(sigma_eta);

    // Inverse Gamma
    logprior += dinvgamma(sigma_eta,2,1,1);
  }
  // eta
  for(int i=0; i<n_eta; i++ ){
    logprior += R::dnorm( eta[i], 0, sigma_eta, 1 );
  }
  //}

  return logprior;
}

// [[Rcpp::export]]
arma::colvec loglik_HM_i_cpp( const arma::colvec Y,
                              const arma::mat X,
                              const arma::colvec beta,
                              const double sigma,
                              const double v,
                              const arma::colvec ind_v ) {

  // Y_i's indep

  // HM0:
  //    normd15N ∼ 1 + std::log(Rainfall) + ManureLevel
  //    weights = varIdent(form=~1|Category)

  unsigned int n_obs = X.n_rows;
  unsigned int p_vars = X.n_cols;

  arma::colvec loglik = arma::zeros<arma::colvec>(n_obs);

  if( Y.n_rows!=n_obs ) {
    throw std::range_error("There's a problem with the data passed to loglik_HM (code 1)");
  }
  if( beta.n_rows!=p_vars ) {
    throw std::range_error("There's a problem with the data passed to loglik_HM (code 2)");
  }
  if( ind_v.n_rows!=n_obs ) {
    throw std::range_error("There's a problem with the data passed to loglik_HM (code 3)");
  }

  // check that all sigma's are positive
  if( (sigma<=0) || (v<=0) ) {
    loglik.fill( -INFINITY );
    return loglik;
  }

  arma::vec lin_pred = X*beta;

  for ( unsigned int i=0; i<n_obs; i++ ) {
    if( ind_v(i)==0 ) {
      loglik(i) = R::dnorm( Y(i), lin_pred(i), sigma , 1);
    } else {
      loglik(i) = R::dnorm( Y(i), lin_pred(i), sqrt(v)*sigma, 1 );
    }
  }
  return loglik;

}

// [[Rcpp::export]]
double loglik_HM_cpp( const arma::colvec Y,
                      const arma::mat X,
                      const arma::colvec beta,
                      const double sigma,
                      const double v,
                      const arma::colvec ind_v ) {

  unsigned int n_obs = X.n_rows;
  arma::colvec loglik = arma::zeros<arma::colvec>(n_obs);

  loglik = loglik_HM_i_cpp( Y, X,
                            beta, sigma,
                            v, ind_v );
  return sum( loglik );

}

// [[Rcpp::export]]
double logprior_HM_cpp( const arma::colvec beta,
                        const double sigma,
                        const double v,
                        const double sigma_eta,
                        const arma::colvec eta,
                        const std::string prior_spec ) {

  unsigned int n_beta;
  n_beta = beta.size();

  unsigned int n_eta;
  n_eta = eta.size();

  double logprior;
  logprior = 0;

  // check that all sigma's are positive
  if( (sigma<=0) || (sigma_eta <=0) || (v<=0) ) {
    return -INFINITY;
  }

  if (prior_spec=="flat") {
    // beta's
    logprior += 0;
    // sigma~1/sigma
    // [Hoff, 2009], section 9.2.2
    logprior += -std::log(sigma);

    // var_offset
    logprior += -std::log(v);
  } else if (prior_spec=="proper") {

    // beta's~norm(0,10) i.i.d.
    double var_beta =10;
    // [Hoff, 2009], section 9.2.2
    for( unsigned int i=0; i<n_beta; i++ ) {
      logprior += R::dnorm( beta[i], 0, sqrt(var_beta), 1 );
    }

    // sigma~1/sigma
    // [Hoff, 2009], section 9.2.2
    logprior += -std::log(sigma);

    // var_offset
    logprior += -std::log(v);
  } else {
    throw std::range_error("There's a problem with prior_spec in HM module");
  }


  // Random effects
  //if( !is.null(sigma_eta) ) {

  // sigma_eta: variance of random effects

  // Jeffrey's prior
  // logprior += -std::log(sigma_eta);

  // Inverse Gamma
  logprior += dinvgamma(sigma_eta,2,1,1);

  // eta
  for( unsigned int i=0; i<n_eta; i++ ) {
    logprior += R::dnorm( eta[i], 0, sigma_eta, 1 );
  }
  //}

  return logprior;
}
