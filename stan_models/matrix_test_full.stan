# model for presence-absence data
data {
  int<lower=1> nKnots;
  int<lower=1> nLocs;
  real y[10];
  matrix[nKnots,nKnots] distKnotsSq;
  matrix[nLocs*10,nKnots] distKnots21Sq;
}
parameters {
  real<lower=0> sigma;
}
transformed parameters {
  matrix[nLocs*10,nKnots] SigmaOffDiag;
	SigmaOffDiag = distKnots21Sq * distKnotsSq;
}
model {
  sigma ~ normal(0,1);
  for(i in 1:10) {
    y[i] ~ normal(0,sigma);
  }
}

