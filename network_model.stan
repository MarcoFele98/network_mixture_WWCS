data {
  int<lower=0> N; // number of data points
  int<lower=0> N_ind; // number of individuals
  int<lower=0> n_hours_same[N]; // number of hours in which both individuals are performing the same behaviour. 
  int<lower=0> id_1[N]; // first individual of edge
  int<lower=0> id_2[N]; // second individual of edge
  vector<lower=0,upper=1>[N] dominance_1; // dominance frst individual
  vector<lower=0,upper=1>[N] dominance_2; // dominance second individual
}

transformed data {
  // Centralize predictors to increase sampling performance. Response should not be centralized as it is a random variable (mean() would only give an estimate of the population mean), and will not improve sampling performance
  vector[N] dominance_1_c = dominance_1 - mean(dominance_1);
  vector[N] dominance_2_c = dominance_2 - mean(dominance_2);
}

parameters {
  // intercepts
  ordered[2] alpha; // funky type which indicates a vector such that alpha[1] < alpha[2]. This is needed for identifiability of mixture models, and is also useful for ordered logistic regression
  // fixed effect
  real beta; // effect of dominance
  // hyperparameters
  real<lower=0> sigma_id; // sd of individual identity effect
  // mixture mixing parameter
  real<lower=0,upper=1> theta; // probability the individuals are in the same location
  // random effects
  vector[N_ind] id_z; // individual identity effect
}

transformed parameters { // Back-transform centered parametrization. id ~ N(0, 1). 
  vector[N_ind] transformed_id = sigma_id * id_z + alpha[2];
}

model {
  // Priors
  // fixed effect
  alpha[1] ~ normal(-1, 5);
  beta ~ normal(0, 5);
  // Hyperparameters
  alpha[2] ~ normal(0, 5);
  sigma_id ~ exponential(0.1);
  // mixture mixing probability
  theta ~ beta(4, 4);
  // random effects 
  target += normal_lpdf(id_z | 0, 1); // target += log_probability and random_variable ~ distribution() are similar. The second calculates the log likelihood up to a scaling constant and is faster, but the likelihood (for example, used in model selection) will only be correct up to a scaling constant
  
  // Model
  for(i in 1:N) { // for every observation
    real p = transformed_id[id_1[i]] + transformed_id[id_2[i]] + beta * (dominance_1_c[i] + dominance_2_c[i]); 
    target += log_sum_exp( // marginalize likelihood. log_sum_exp() is the numerically stable way to calculate log(probability1 + probability2)
    log(theta) + binomial_logit_lpmf(n_hours_same[i] | 24, p), // "process" 1 : the pair is in the same location. binomial_logit_lpmf(k | n, p) = binomial_lpmf(k | n, inv_logit(p)) = binomial_lpmf(k | n, exp(p)/(1+exp(p))) 
    log1m(theta) + binomial_logit_lpmf(n_hours_same[i] | 24, alpha[1]) // "process" 2 : the pair is in a different location. log1m() is the numerically stable way to calculate 1 - probability
    );
  }
}

generated quantities {
  // Classify observation as same or different location
  vector[N] prob_same;
  for(i in 1:N) { // the mixture probability of a data point is the likelihood mulitplied by the mixture probability normalized over the mixtures
    real p = transformed_id[id_1[i]] + transformed_id[id_2[i]] + beta * (dominance_1_c[i] + dominance_2_c[i]); 
    real mixture_1 = binomial_logit_lpmf(n_hours_same[i] | 24, p) + log(theta);
    real mixture_2 = binomial_logit_lpmf(n_hours_same[i] | 24, alpha[1]) + log1m(theta);
    prob_same[i] = exp(mixture_1 - log_sum_exp(mixture_1, mixture_2)); 
  }
  // Generate edge posterior
  matrix[N_ind, N_ind] edges;
  for(ind_1 in 1:N_ind) {
    for(ind_2 in 1:N_ind) {
      edges[ind_1, ind_2] = inv_logit(transformed_id[ind_1] + transformed_id[ind_2] + beta * (dominance_1_c[ind_1] + dominance_2_c[ind_2]));
    }                                                      
  }
  // Generate centrality posterior
  vector[N_ind] centrality;
  for(ind_1 in 1:N_ind) {
    centrality[ind_1] = 0;
    for(ind_2 in 1:N_ind) {
      centrality[ind_1] += edges[ind_1, ind_2];
    }
  }
  centrality /= N_ind;
}
