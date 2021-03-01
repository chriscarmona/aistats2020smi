
// saved as plummer_hpv.stan

data {
  int<lower=0> n_obs; // number of observations

  int<lower=0> ncases[n_obs];
  real<lower=0> Npop[n_obs];
  
  real<lower=0,upper=1> phi[n_obs];
}

parameters {
  real theta1;
  real<lower=0> theta2;
}

transformed parameters {
  real<lower=0> mu[n_obs];
  for (i in 1:n_obs) {
    mu[i] = (Npop[i]/1000)* exp(theta1 + phi[i] * theta2);
  }
}

model {
  // The likelihood
  for (i in 1:n_obs) {
    ncases[i] ~ poisson( mu[i] );
  }
}
