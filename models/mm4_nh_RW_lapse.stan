// A rescrola wagner model with lapse term
// Author: Ondrej Zika, 4th Sept, 2023


data {
  int<lower=1> nSub;     // Number of participants
  int<lower=1> nVis;     // Number of visits
  int<lower=1> maxTr;
  int<lower=0, upper=1> O[nVis, nSub, maxTr]; // Delivered outcomes
  real V[nVis, nSub, maxTr];      // Values reported by participants
  int<lower=1> nObserved[nVis, nSub]; // Variable containing number of trials
  int<lower=0> nOuts; // number of outcome types (should be 2 for sh/nosh)
}

parameters {
  real<lower=0, upper=1> alpha[nVis, nSub];  // Learning rate for each participant, visit
  real<lower=0, upper=1> lapse[nVis, nSub];  // Lapse term
  real<lower=0, upper=1> Vpred[nVis, nSub, maxTr]; // model prediction for each trial
  real<lower=0> sigma_V[nSub];       // reporting noise, separate for each participant 
}

transformed parameters {
  real pe[nVis, nSub, maxTr];
  for (v in 1:nVis) {
    for (s in 1:nSub) {
      for (t in 1:nObserved[v,s]-1) {
        pe[v,s,t] = O[v,s,t] - V[v,s,t];
      }
    }
  }
}

model {
  for (v in 1:nVis) {
    for (s in 1:nSub) {
      // Prior on subject-level alpha
  
      alpha[v,s] ~ beta(1,1);
      lapse[v,s] ~ beta(1,1);

      sigma_V[s] ~ cauchy(0, 2.5); 
      Vpred[v, s, 1] ~ beta(1, 1);
      
      for (t in 1:nObserved[v,s]-1) {
        V[v,s,t+1] ~ normal((1-lapse[v,s])*(V[v,s,t] + alpha[v,s]*pe[v,s,t]) + (0.5*lapse[v,s]), sigma_V[s]);
        // Also gather predictions
        Vpred[v,s,t+1] ~ normal((1-lapse[v,s])*(V[v,s,t] + alpha[v,s]*pe[v,s,t]) + (0.5*lapse[v,s]), sigma_V[s]);
      }
    }
  }
}

generated quantities {
  real log_lik[nVis, nSub, maxTr];  // Log likelihood for each trial
  for (v in 1:nVis) {
    for (s in 1:nSub) {
      for (t in 1:nObserved[v,s]) {
        // Probability of the actual data given the model (P(D|M))
        log_lik[v,s,t] = normal_lpdf(V[v,s,t] | Vpred[v,s,t] , sigma_V[s]);
      }
    }
  };
  
}