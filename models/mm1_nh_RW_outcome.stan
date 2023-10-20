// A rescrola wagner model with separate learning rates for shock and no-shock

// Author: Ondrej Zika, 4th Sept, 2023

data {
  int<lower=1> nSub;     // Number of participants
  int<lower=1> nVis;     // Number of visits
  int<lower=1> maxTr;     // maxTrials in a given visit
  int<lower=0, upper=1> O[nVis, nSub, maxTr]; // Delivered outcomes
  real V[nVis, nSub, maxTr];      // Values reported by participants
  int<lower=1> nObserved[nVis, nSub]; // Variable containing number of trials 
  int<lower=0> nOuts; // number of outcome types (should be 2 for sh/nosh)
}

parameters {
  real<lower=0, upper=1> alpha[nVis, nSub, nOuts]; // Learning rate for each participant, visit and outcome
  real<lower=0, upper=1> Vpred[nVis, nSub, maxTr]; // model prediction for each trial
  real<lower=0> sigma_V[nSub];                     // reporting noise, separate for each participant 

}

transformed parameters {
  real pe[nVis, nSub, maxTr]; // Prediction error
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
      // Prior on alpha
      for (o in 1:nOuts) {
        alpha[v,s,o] ~ beta(1,1);
      }
      // Prior on SD
      sigma_V[s] ~ cauchy(0, 2.5); 
      
      // Prior on first value
      Vpred[v, s, 1] ~ beta(1, 1);
      
      for (t in 1:nObserved[v,s]-1) {
          if (O[v,s,t] == 0) {
            V[v,s,t+1] ~ normal(V[v,s,t] + alpha[v,s,1]*pe[v,s,t], sigma_V[s]);
            // Also gather predictions
            Vpred[v,s,t+1] ~ normal(V[v,s,t] + alpha[v,s,1]*pe[v,s,t], sigma_V[s]);
          } else if (O[v,s,t] == 1) {
            V[v,s,t+1] ~ normal(V[v,s,t] + alpha[v,s,2]*pe[v,s,t], sigma_V[s]);
            // Also gather predictions
            Vpred[v,s,t+1] ~ normal(V[v,s,t] + alpha[v,s,2]*pe[v,s,t], sigma_V[s]);
          }
      }
    }
  }
}

generated quantities {
  real log_lik[nVis, nSub, maxTr];  // Log likelihood for each trial
  for (v in 1:nVis) {
    for (s in 1:nSub) {
      for (t in 1:nObserved[v,s]) {
        // Probability of the data given the model/parameters (P(D|M))
        log_lik[v,s,t] = normal_lpdf(V[v,s,t] | Vpred[v,s,t] , sigma_V[s]);
      }
    }
  };
  
}