# model for growth
data {
  int<lower=1> N; # individuals
  real y[N,2]; #
  matrix[N,2] ages;
}
parameters {
  real sigmaA;
  real sigmaB;
  real A[2]; #asymptotic adult length, by sex
  real<lower=0> b[2]; # slope of curve, by sex
  real<lower=0> c[2]; # inflection point of curve, by sex
  real<lower=0> m[2]; # rel position, by sex
  real theta_a[ind]; # random effect in max length
  real theta_b[ind]; # random effect in growth rate
}
model {
  # priors
  A ~ normal(0, 2);
  b ~ normal(0, 2);
  c ~ normal(0, 2);
  m ~ normal(0, 2);
  sigmaA ~ cauchy(0,5);
  sigmaB ~ cauchy(0,5);
  theta_a ~ normal(0, sigmaA); # random effect in max length
  theta_b ~ normal(0, sigmaB); # random effect in max length

  for(n in 1:N) {
	  y[n] ~ gamma(gammaA, gammaA/exp(fmin(spatialEffects[1, location[n]], 200)));
  }
}

