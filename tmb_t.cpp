#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
DATA_VECTOR(y);

PARAMETER(log_df);
//PARAMETER(location);
//PARAMETER(log_scale);

PARAMETER_VECTOR(w);

int n = y.size();

Type nll = 0.0;
Type location = 0.0;

for(int i = 0; i < n; i++) {
  Type df = exp(log_df);
  nll -= dgamma(exp(w[i]), df/2.0, df/2.0, true);
  //nll -= dnorm(y[i], location, sqrt(1/exp(w[i]) * pow(exp(log_scale),2)), true);
  nll -= dnorm(y[i], location, sqrt(1/exp(w[i])), true);
  //nll -= dt(y[i], df, true);
}

REPORT(w);
ADREPORT(w);

return nll;
}
