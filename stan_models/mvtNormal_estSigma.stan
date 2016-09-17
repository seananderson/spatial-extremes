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
  real<lower=1> scaledf;
  real<lower=0> sigma;
  vector[nKnots] spatialEffectsKnots[nT];
}
transformed parameters {
	vector[nKnots] muZeros;
	vector[nLocs] spatialEffects[nT];
  matrix[nKnots, nKnots] SigmaKnots;
  matrix[nLocs,nKnots] SigmaOffDiag;
  matrix[nLocs,nKnots] invSigmaKnots;
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
  gp_scale ~ normal(0,1);
  gp_sigmaSq ~ normal(0,1);
  sigma ~ normal(0,1);
  scaledf ~ gamma(2,0.1); # prior from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  for(t in 1:nT) {
  spatialEffectsKnots[t] ~ multi_student_t(scaledf, muZeros, SigmaKnots);
  }

  for(t in 1:nT) {
    y[t] ~ normal(spatialEffects[t], sigma);
  }

}

