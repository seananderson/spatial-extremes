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
  real<lower=2> scaledf;
  real<lower=0> CV;
  vector[nKnots] spatialEffectsKnots[nT];
}
transformed parameters {
	vector[nKnots] muZeros;
	vector[nLocs] spatialEffects[nT];
  matrix[nKnots, nKnots] SigmaKnots;
  matrix[nLocs,nKnots] SigmaOffDiag;
  matrix[nLocs,nKnots] invSigmaKnots;
  real<lower=0> gammaA;
  SigmaKnots = gp_sigmaSq * exp(-gp_scale * distKnotsSq);// cov matrix between knots
  SigmaOffDiag = gp_sigmaSq * exp(-gp_scale * distKnots21Sq);// cov matrix between knots and projected locs
	for(i in 1:nKnots) {
		muZeros[i] = 0;
	}
	SigmaOffDiag = SigmaOffDiag * inverse_spd(SigmaKnots); // multiply and invert once, used below
	for(i in 1:nT) {
  spatialEffects[i] = SigmaOffDiag * spatialEffectsKnots[i];
	}
	gammaA = 1/(CV*CV);
}
model {
  // priors on parameters for covariances, etc
  gp_scale ~ student_t(5, 0, 20);
  gp_sigmaSq ~ student_t(5, 0, 2);
  CV ~ normal(0, 1); // TODO why was N(0, 1)?
  scaledf ~ gamma(2, 0.1); // prior from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  for(t in 1:nT) {
    spatialEffectsKnots[t] ~ multi_student_t(scaledf, muZeros, SigmaKnots);
  }

  for(t in 1:nT) {
    for(l in 1:nLocs) {
      y[t,l] ~ gamma(gammaA, gammaA/exp(spatialEffects[t,l]));
    }
  }

}

