#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(r, cpp)]]

//' @export
// [[Rcpp::export]]
Rcpp::List  SMI_post_biased_data( const arma::vec Z,
                                  const arma::vec Y,
                                  const double sigma_z,
                                  const double sigma_y,
                                  const double sigma_phi,
                                  const double sigma_theta,
                                  const double sigma_theta_tilde,
                                  const double eta=1 ) {
  /*
  Calculates the posterior mean and variance of the Semi-Modular Posterior
  P_{smi,\eta}( \varphi, \theta, \tilde\theta | Z, Y )
  
  Z, Y : data from the unbiased and biased modules, respectively
  eta : degree of influence of the unbiased module into phi
  sigma_z, sigma_y : likelihood std for Z and Y, respectively
  sigma_phi, sigma_theta : std dev for the priors of phi, theta, theta_tilde respectively
  */
  unsigned int n=Z.n_rows;
  unsigned int m=Y.n_rows;
  
  arma::mat A = arma::zeros<arma::mat>(3,3);
  
  A(0,0) = n/pow(sigma_z,2) + (1+eta) * m/pow(sigma_y,2) - ( m / ( pow(sigma_y,2) + m*pow(sigma_theta,2) ) ) + 1/pow(sigma_phi,2);
  A(1,1) = m/pow(sigma_y,2) + 1/pow(sigma_theta,2);
  A(2,2) = eta * m/pow(sigma_y,2) + 1/pow(sigma_theta_tilde,2);
  A(0,1) = A(1,0) = m/pow(sigma_y,2);
  A(0,2) = A(2,0) = eta*m/pow(sigma_y,2);
  
  arma::mat b = arma::zeros<arma::mat>(3,1);
  b(0,0) = sum(Z)/pow(sigma_z,2) + (1+eta) * sum(Y)/pow(sigma_y,2) - sum(Y) / (pow(sigma_y,2)+m*pow(sigma_theta,2));
  b(1,0) = sum(Y)/pow(sigma_y,2);
  b(2,0) = eta*sum(Y)/pow(sigma_y,2);
  
  arma::mat cov = A.i();
  arma::mat mean = cov * b;
  
  return Rcpp::List::create( Rcpp::Named("mean") = mean,
                             Rcpp::Named("cov") = cov );
}

//' @export
// [[Rcpp::export]]
Rcpp::List  SMI_pred_biased_data( const arma::vec Z,
                                  const arma::vec Y,
                                  const double sigma_z,
                                  const double sigma_y,
                                  const double sigma_phi,
                                  const double sigma_theta,
                                  const double sigma_theta_tilde,
                                  const double eta=1 ) {
  /*
   Calculates the posterior mean and variance of the Semi-Modular Posterior
   P_{smi,\eta}( \varphi, \theta, \tilde\theta | Z, Y )
   
   Z, Y : data from the unbiased and biased modules, respectively
   eta : degree of influence of the unbiased module into phi
   sigma_z, sigma_y : likelihood std for Z and Y, respectively
   sigma_phi, sigma_theta : std dev for the priors of phi, theta, theta_tilde respectively
   */
  unsigned int n=Z.n_rows;
  unsigned int m=Y.n_rows;
  
  arma::mat A = arma::zeros<arma::mat>(5,5);
  
  A(0,0) = (n+1)/pow(sigma_z,2) + (m+eta*m+1)/pow(sigma_y,2) - ( m / ( pow(sigma_y,2) + m*pow(sigma_theta,2) ) ) + 1/pow(sigma_phi,2);
  A(1,1) = (m+1)/pow(sigma_y,2) + 1/pow(sigma_theta,2);
  A(2,2) = eta * m/pow(sigma_y,2) + 1/pow(sigma_theta_tilde,2);
  A(3,3) = 1/pow(sigma_z,2);
  A(4,4) = 1/pow(sigma_y,2);
  
  A(0,1) = A(1,0) = (m+1)/pow(sigma_y,2);
  A(0,2) = A(2,0) = eta*m/pow(sigma_y,2);
  A(0,3) = A(3,0) = -1/pow(sigma_z,2);
  A(0,4) = A(4,0) = -1/pow(sigma_y,2);
  A(1,4) = A(4,1) = -1/pow(sigma_y,2);
  
  arma::mat b = arma::zeros<arma::mat>(5,1);
  b(0,0) = sum(Z)/pow(sigma_z,2) + (1+eta) * sum(Y)/pow(sigma_y,2) - sum(Y) / (pow(sigma_y,2)+m*pow(sigma_theta,2));
  b(1,0) = sum(Y)/pow(sigma_y,2);
  b(2,0) = eta*sum(Y)/pow(sigma_y,2);
  
  arma::mat cov = A.i();
  arma::mat mean = cov * b;
  
  return Rcpp::List::create( Rcpp::Named("mean") = mean.rows(3,4),
                             Rcpp::Named("cov") = cov.submat(3,3,4,4) );
}
