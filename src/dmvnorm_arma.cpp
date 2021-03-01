#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::interfaces(r, cpp)]]

const double log2pi = std::log(2.0 * M_PI);

//' @title Multivariate Normal Density Function
//' @description This function returns the log-density for a multivariate Gaussian distribution. Originally from: https://gallery.rcpp.org/articles/dmvnorm_arma/
//' @param x matrix of observations
//' @param mean mean vector
//' @param sigma positive definite covariance matrix
//' @param logd logical; whether log-density should be returned (default to FALSE)
//' @export
// [[Rcpp::export]]
arma::vec dmvnorm_arma( arma::mat x,
                        arma::rowvec mean,
                        arma::mat sigma,
                        bool logd = false) {
  unsigned int n = x.n_rows;
  unsigned int xdim = x.n_cols;
  unsigned int i;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}
