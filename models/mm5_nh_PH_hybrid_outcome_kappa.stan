// A hierarchical pearce-hall - rescrola-wagner hybrid

// Author: Ondrej Zika, 25th August, 2023

data {
  int<lower=1> nSub;     // Number of participants
  int<lower=1> nVis;     // Number of visits
  int<lower=1> nOutc;     // Number of visits
  int<lower=1> maxTr;
  int<lower=0, upper=1> O[nVis, nSub, maxTr]; // Delivered outcomes
  real V[nVis, nSub, maxTr];      // Values reported by participants
  int<lower=1> nObserved[nVis, nSub]; // Variable containing number of trials
}

parameters {
  real<lower=0, upper=1> eta[nVis, nSub, nOutc];  // Eta voaltility param for each participant, visit and outcome
  real<lower=0, upper=1> alpha0[nVis, nSub];  // starting alpha
  real<lower=0, upper=1> kappa[nVis, nSub];  // starting alpha
  real<lower=0, upper=1> Vpred[nVis, nSub, maxTr]; // model prediction for each trial
  real<lower=0> sigma_V[nSub];       // reporting noise, separate for each participant 
}

transformed parameters {
  real pe[nVis, nSub, maxTr];
  real abspe[nVis, nSub, maxTr];
  real<lower=0, upper=1> alpha[nVis, nSub, maxTr]; // model prediction for each trial
  for (v in 1:nVis) {
    for (s in 1:nSub) {
      for (t in 1:nObserved[v,s]) {
        pe[v,s,t] = O[v,s,t] - V[v,s,t];
        abspe[v,s,t] = fabs(O[v,s,t] - V[v,s,t]);
        
      }
      alpha[v,s,1] = alpha0[v,s];
      for (t in 2:maxTr) {
            if (t <= nObserved[v,s]) { 
                alpha[v,s,t] = kappa[v,s]*(eta[v,s, O[v,s,t-1]+1]*alpha[v,s,t-1] + (1-eta[v,s, O[v,s,t-1]+1])*abspe[v,s,t-1]);
            } else {
              alpha[v,s,t] = 0;
            }
        }
    }
  }
}

model {

  for (v in 1:nVis) {

    for (s in 1:nSub) {
      // Prior on subject-level alpha
      for (o in 1:nOutc) {
        eta[v,s, o] ~ beta(1,1);
      }

      sigma_V[s] ~ cauchy(0, 2.5); 
      Vpred[v, s, 1] ~ beta(1, 1);
      alpha0[v,s] ~ beta(1, 1);
      kappa[v,s] ~ beta(1, 1);
      
      for (t in 1:nObserved[v,s]-1) {
        
        V[v,s,t+1] ~ normal(V[v,s,t] + alpha[v,s,t]*pe[v,s,t], sigma_V[s]);
        // Also gather predictions
        Vpred[v,s,t+1] ~ normal(V[v,s,t] + alpha[v,s,t]*pe[v,s,t], sigma_V[s]);
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