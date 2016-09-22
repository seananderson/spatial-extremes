# model for presence-absence data
data {
  int N;
  real y[N]; // data
}
parameters {
  real mu;
  real<lower=1> scaledf;
  real<lower=0> sigma;
  real<lower=0> W[N];
}
model {
  # priors on parameters for covariances, etc
  mu ~ normal(0,1);
  sigma ~ normal(0,1);
  scaledf ~ exponential(0.01);

  # commented out line gives same result
  #W ~ scaled_inv_chi_square(scaledf,1); # see discussion https://groups.google.com/forum/#!topic/stan-users/0F0O4hfHA8g
  #W ~ inv_gamma(scaledf/2, scaledf/2);
  W ~ gamma(scaledf/2, scaledf/2);
  for(t in 1:N) {
  y[t] ~ normal(0, sqrt((1/W[t])*sigma*sigma));
  }
}

