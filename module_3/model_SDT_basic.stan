data{
  int<lower=1> nTrials; // number of trials total
  int<lower=1, upper=2> stimidx[nTrials]; // stimulus index
  int<lower=0, upper=1> resp[nTrials]; // response
  real dprime_range[2]; // lower and bound, dprime estimate
  real k_range[2]; // lower and bound, criterion estimate
}

parameters{
  real<lower=dprime_range[1], upper=dprime_range[2]> dprime; // d' estimate for subject
  real<lower=k_range[1], upper=k_range[2]> k; // criterion estimate for subject
}

model{
  vector[2] pCW; // probability of a clockwise judgement for each stimulus
  dprime ~ uniform(dprime_range[1],dprime_range[2]); // flat prior on d'
  k ~ uniform(k_range[1],k_range[2]); // flat prior on criterion 
  pCW[1] = 1 - normal_cdf(k,-0.5*dprime,1); // probability clockwise for CCW
  pCW[2] = 1 - normal_cdf(k,0.5*dprime,1); // probability clockwise for CW
  for (n in 1:nTrials) {
    resp[n] ~ bernoulli(pCW[stimidx[n]]); // likelihood function
  }
}
