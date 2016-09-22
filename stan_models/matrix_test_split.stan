# model for presence-absence data
data {
  int<lower=1> nKnots;
  int<lower=1> nLocs;
  real y[10];
  matrix[nKnots,nKnots] distKnotsSq;
  matrix[nLocs,nKnots] distKnots21Sq1;
  matrix[nLocs,nKnots] distKnots21Sq2;
  matrix[nLocs,nKnots] distKnots21Sq3;
  matrix[nLocs,nKnots] distKnots21Sq4;
  matrix[nLocs,nKnots] distKnots21Sq5;
  matrix[nLocs,nKnots] distKnots21Sq6;
  matrix[nLocs,nKnots] distKnots21Sq7;
  matrix[nLocs,nKnots] distKnots21Sq8;
  matrix[nLocs,nKnots] distKnots21Sq9;
  matrix[nLocs,nKnots] distKnots21Sq10;
}
parameters {
  real<lower=0> sigma;
}
transformed parameters {
  matrix[nLocs,nKnots] SigmaOffDiag1;
  matrix[nLocs,nKnots] SigmaOffDiag2;
  matrix[nLocs,nKnots] SigmaOffDiag3;
  matrix[nLocs,nKnots] SigmaOffDiag4;
  matrix[nLocs,nKnots] SigmaOffDiag5;
  matrix[nLocs,nKnots] SigmaOffDiag6;
  matrix[nLocs,nKnots] SigmaOffDiag7;
  matrix[nLocs,nKnots] SigmaOffDiag8;
  matrix[nLocs,nKnots] SigmaOffDiag9;
  matrix[nLocs,nKnots] SigmaOffDiag10;
	SigmaOffDiag1 = distKnots21Sq1 * distKnotsSq;
	SigmaOffDiag2 = distKnots21Sq2 * distKnotsSq;
	SigmaOffDiag3 = distKnots21Sq3 * distKnotsSq;
	SigmaOffDiag4 = distKnots21Sq4 * distKnotsSq;
	SigmaOffDiag5 = distKnots21Sq5 * distKnotsSq;
	SigmaOffDiag6 = distKnots21Sq6 * distKnotsSq;
	SigmaOffDiag7 = distKnots21Sq7 * distKnotsSq;
	SigmaOffDiag8 = distKnots21Sq8 * distKnotsSq;
	SigmaOffDiag9 = distKnots21Sq9 * distKnotsSq;
	SigmaOffDiag10 = distKnots21Sq10 * distKnotsSq;
}
model {
  sigma ~ normal(0,1);
  for(i in 1:10) {
    y[i] ~ normal(0,sigma);
  }
}

