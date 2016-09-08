# model for presence-absence data
data {
  int<lower=1> nKnots; // dimension of knots
  int<lower=1> nT; // number of observations
  matrix[nKnots,nKnots] distKnotsSq;
  matrix[nT,nKnots] y; // data
}
parameters {
  real<lower=0> gp_scale;
  real<lower=0> gp_sigmaSq;
  real<lower=0> scaledf;
}
transformed parameters {
  vector[nKnots] muZeros;
  matrix[nKnots, nKnots] SigmaKnots;
  SigmaKnots = gp_sigmaSq * exp(-gp_scale * distKnotsSq);# cov matrix between knots
  for(i in 1:nKnots) {
	muZeros[i] = 0;
  }
}
model {
  # priors on parameters for covariances, etc
  gp_scale ~ cauchy(0,5);
  gp_sigmaSq ~ cauchy(0,5);
  scaledf ~ exponential(0.1);
  for(t in 1:nT) {
  y[t] ~ multi_student_t(scaledf, muZeros, SigmaKnots);
  }
}

