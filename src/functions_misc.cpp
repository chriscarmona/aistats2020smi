#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
arma::mat categ_to_dummie( const arma::vec x,
                           const arma::vec categ_x ) {
  // this function transforms a categorical variable into a matrix with dummie columns, one for each category in categ_x
  arma::mat dummie_mat = arma::zeros(x.n_rows,categ_x.n_rows);
  arma::vec dummie_i;
  for( unsigned int i=0; i<categ_x.n_rows; i++ ) {
    dummie_i = arma::zeros(x.n_rows,1);
    dummie_i.elem(find(x==categ_x(i))).ones();
    dummie_mat.col(i) = dummie_i;
  }
  return dummie_mat;
}

// [[Rcpp::export]]
double sign( const double val ) {
  return (0 < val) - (val < 0);
}

// [[Rcpp::export]]
double bounce_limit( double x,
                     const double a,
                     const double b) {
  while( (x<a) || (x>b) ) {
    if(x < a) {
      x = a + (a-x);
    }
    if(x > b) {
      x = b - (x-b);
    }
  }
  return x;
}

// [[Rcpp::export]]
double dinvgamma( const double x,
                  const double alpha, //shape
                  const double beta, //scale
                  const unsigned int lg=0 ) {
  double logpdf;

  if ( (alpha <= 0.0) | (beta <= 0.0) | (x<= 0.0) ) {
    if( lg==1 ) {
      return -INFINITY;
    } else if ( lg == 0 ) {
      return 0.0;
    }
  }

  logpdf = alpha * std::log(beta) - R::lgammafn(alpha) - (alpha + 1.0) * std::log(x) - (beta/x);

  if( lg==1 ) {
    return logpdf;
  } else if ( lg == 0 ) {
    return std::exp(logpdf);
  }

  return NAN;
}
