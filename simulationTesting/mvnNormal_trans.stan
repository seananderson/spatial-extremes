# model for presence-absence data
data {
  int<lower=1> nKnots;
  int<lower=1> nLocs;
  int<lower=1> nT;
  matrix[nT,nLocs] y;
  matrix[nKnots,nKnots] distKnotsSq;
  matrix[nLocs,nKnots] distKnots21Sq;
}
parameters {
  real<lower=0> gp_scale;
  real<lower=0> gp_sigmaSq;
  real<lower=0> sigma;
  real<lower=0> scaledf;
  real<lower=0> W[nT];
  vector[nKnots] spatialEffectsKnots[nT];
}
transformed parameters {
	vector[nKnots] muZeros;
	vector[nLocs] spatialEffects[nT];
  matrix[nKnots, nKnots] SigmaKnots;
  matrix[nKnots, nKnots] SigmaKnots_chol;
  matrix[nLocs,nKnots] SigmaOffDiag;

  SigmaKnots = gp_sigmaSq * exp(-gp_scale * distKnotsSq);# cov matrix between knots
  SigmaOffDiag = gp_sigmaSq * exp(-gp_scale * distKnots21Sq);# cov matrix between knots and projected locs
	for(i in 1:nKnots) {
		muZeros[i] = 0;
	}
	SigmaOffDiag = SigmaOffDiag * inverse(SigmaKnots); # multiply and invert once, used below
	for(i in 1:nT) {
  spatialEffects[i] = SigmaOffDiag * spatialEffectsKnots[i];
	}
}
model {
  # priors on parameters for covariances, etc
  gp_scale ~ cauchy(0,5);
  gp_sigmaSq ~ cauchy(0,5);
  sigma ~ cauchy(0,5);
  scaledf ~ exponential(0.01);
  W ~ inv_gamma(scaledf/2, gp_sigmaSq*scaledf/2);

  for(t in 1:nT) {
  spatialEffectsKnots[t] ~ multi_normal(muZeros, W[t]*SigmaKnots);
  }

  for(t in 1:nT) {
    y[t] ~ normal(spatialEffects[t], 0.00001);
    #y[t] ~ normal(spatialEffects[t], sigma);
  }

}

