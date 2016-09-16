# model for presence-absence data
data {
  int<lower=1> nKnots;
  int<lower=1> nLocs;
  int<lower=1> N;
  int<lower=1> nT;
  real y[N];
  int location[N];
  matrix[nKnots,nKnots] distKnotsSq;
  matrix[nLocs,nKnots] distKnots21Sq;
}
parameters {
  real<lower=0> gp_scale;
  real<lower=0> gp_sigmaSq;
  real<lower=0> scaledf;
  real<lower=0> gammaA;
  real<lower=0,upper=40> dft;
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
  spatialEffects[1] = SigmaOffDiag * inverse(SigmaKnots) * (spatialEffectsKnots[1]);
}
model {
  spatialEffectsKnots[1] ~ multi_student_t(dft, muZeros, SigmaKnots);

  # priors on parameters for covariances, etc
  dft ~ cauchy(0,5);
  gp_scale ~ cauchy(0,5);
  gp_sigmaSq ~ cauchy(0,5);
  gammaA ~ cauchy(0,5);
  for(n in 1:N) {
	  y[n] ~ gamma(gammaA, gammaA/exp(fmin(spatialEffects[1, location[n]], 200)));
  }
}

