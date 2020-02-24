
// saved as plummer_hpv.stan

data {
  int<lower=0> n_obs; // number of observations

  int<lower=0> nhpv[n_obs];
  int<lower=0> Npart[n_obs];
  int<lower=0> ncases[n_obs];
  real<lower=0> Npop[n_obs];
  
  real<lower=0,upper=1> eta_pois;
  real<lower=0,upper=1> eta_binom;
}

parameters {
  real theta1;
  real<lower=0> theta2;
  real<lower=0,upper=1> phi[n_obs];
}

transformed parameters {
  real<lower=0> mu[n_obs];
  for (i in 1:n_obs) {
    mu[i] = (Npop[i]/1000)* exp( theta1 + phi[i] * theta2);
  }
}

model {
  // The likelihood
  for (i in 1:n_obs) {
    if(eta_pois!=0){
      target += eta_pois * poisson_lpmf( ncases[i] | mu[i] );
    }
    if(eta_binom!=0){
      target += eta_binom * binomial_lpmf( nhpv[i] | Npart[i], phi[i] );
    }
  }
}

// generated quantities {
//   real loglik_i[n_obs];
//   for (i in 1:n_obs) {
//     if(eta_pois!=0){
//       loglik_i[i] += poisson_lpmf( ncases[i] | mu[i] );
//     }
//     if(eta_binom!=0){
//       loglik_i[i] += binomial_lpmf( nhpv[i] | Npart[i], phi[i] );
//     }
//   }
// }
