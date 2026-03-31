data{
  int<lower=1> nObs; // number of observers
  int<lower=1> nTrials; // number of trials total
  int<lower=1> stimidx[nTrials]; // stimulus index
  int<lower=0, upper=1> resp[nTrials,nObs]; // response
  real dprime_range[2]; // lower and bound, dprime estimate
  real k_range[2]; // lower and bound, criterion estimate
}

parameters{
  real<lower=0> dprime_pop_sd; // d' sd at the population level (hyperparameter)
  real<lower=0> k_pop_sd; // k sd at the population level (hyperparameter)
  real<lower=dprime_range[1], upper=dprime_range[2]> dprime[nObs]; // d' estimate for subject
  real<lower=k_range[1], upper=k_range[2]> k[nObs]; // criterion estimate for subject
}

model{
  vector[2] pCW; // probability of a clockwise judgement for each stimulus
  dprime_pop_sd ~ lognormal(0,1); // lognormal prior on sd (hyperprior)
  k_pop_sd ~ lognormal(0,1); // lognormal prior on sd (hyperprior)
  for (s in 1:nObs){
    dprime[s] ~ normal(1,dprime_pop_sd); // hierarchical prior on d'
    k[s] ~ normal(0,k_pop_sd); // hierarchical prior on criterion 
    pCW[1] = 1 - normal_cdf(k[s],-0.5*dprime[s],1); // probability clockwise for CCW
    pCW[2] = 1 - normal_cdf(k[s],0.5*dprime[s],1); // probability clockwise for CW
    for (n in 1:nTrials) {
      resp[n,s] ~ bernoulli(pCW[stimidx[n]]); // likelihood function
    }
  }
  
}
