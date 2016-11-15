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
  real<lower=0> W[nT];
  vector[nKnots] spatialEffectsKnots_z[nT]; // has a standard normal prior below
}
transformed parameters {
	vector[nKnots] muZeros;
	vector[nLocs] spatialEffects[nT];
  matrix[nKnots, nKnots] SigmaKnots;
  matrix[nLocs,nKnots] SigmaOffDiag;
  matrix[nLocs,nKnots] invSigmaKnots;
  real<lower=0> gammaA;
  vector[nKnots] spatialEffectsKnots[nT];

  SigmaKnots = gp_sigmaSq * exp(-gp_scale * distKnotsSq);// cov matrix between knots
  SigmaOffDiag = gp_sigmaSq * exp(-gp_scale * distKnots21Sq);// cov matrix between knots and projected locs
	for(i in 1:nKnots) {
		muZeros[i] = 0;
	}
	SigmaOffDiag = SigmaOffDiag * inverse_spd(SigmaKnots); // multiply and invert once, used below

	gammaA = 1/(CV*CV);

	for(t in 1:nT) {
    spatialEffectsKnots[t] = cholesky_decompose(W[t]*SigmaKnots) * spatialEffectsKnots_z[t]; // 2nd part of multivariate Matt trick
    // spatialEffectsKnots[t] = cholesky_decompose(SigmaKnots) * spatialEffectsKnots_z[t]; // 2nd part of multivariate Matt trick
  }

  for(i in 1:nT) {
    spatialEffects[i] = SigmaOffDiag * spatialEffectsKnots[i];
	}

}
model {

  // priors on parameters for covariances, etc
  gp_scale ~ student_t(5, 0, 20);
  gp_sigmaSq ~ student_t(5, 0, 2);
  CV ~ normal(0, 1);
  scaledf ~ gamma(2, 0.1);

  for(t in 1:nT) {
    spatialEffectsKnots_z[t] ~ normal(0, 1); // 1st part of the multivariate Matt trick
  }

  #W ~ inv_gamma(scaledf/2, gp_sigmaSq*scaledf/2);
  W ~ scaled_inv_chi_square(scaledf, 1);


  for(t in 1:nT) {
    for(l in 1:nLocs) {
      y[t,l] ~ gamma(gammaA, gammaA/exp(spatialEffects[t,l]));
    }
  }

}

