
arma::mat categ_to_dummie( const arma::vec x,
                           const arma::vec categ_x );

double sign( const double val );

double bounce_limit( double x,
                     const double a,
                     const double b);

double dinvgamma( const double x,
                  const double alpha,
                  const double beta,
                  const unsigned int lg );

arma::colvec loglik_PO_i_cpp( const arma::colvec Y,
                              const arma::mat X,
                              const arma::colvec alpha,
                              const arma::colvec beta );

double loglik_PO_cpp( const arma::colvec Y,
                      const arma::mat X,
                      const arma::colvec alpha,
                      const arma::colvec beta );

double logprior_PO_cpp( const arma::colvec alpha,
                        const arma::colvec beta,
                        const double sigma_eta,
                        const arma::colvec eta,
                        const std::string prior_spec ) ;

arma::colvec loglik_HM_i_cpp( const arma::colvec Y,
                              const arma::mat X,
                              const arma::colvec beta,
                              const double sigma,
                              const double v,
                              const arma::colvec ind_v );

double loglik_HM_cpp( const arma::colvec Y,
                      const arma::mat X,
                      const arma::colvec beta,
                      const double sigma,
                      const double v,
                      const arma::colvec ind_v );

double logprior_HM_cpp( const arma::colvec beta,
                        const double sigma,
                        const double v,
                        const double sigma_eta,
                        const arma::colvec eta,
                        const std::string prior_spec );
