functions {
  matrix evprob(real z, real alpha, real beta) {
    matrix[2,2] P;
    P[1,1] = (beta/(alpha+beta)) + (alpha/(alpha+beta)*exp(-(alpha+beta)*z));
    P[1,2] = (alpha/(alpha+beta)) - (alpha/(alpha+beta)*exp(-(alpha+beta)*z));
    P[2,1] = (beta/(alpha+beta)) - (beta/(alpha+beta)*exp(-(alpha+beta)*z));
    P[2,2] = (alpha/(alpha+beta)) + (beta/(alpha+beta)*exp(-(alpha+beta)*z));
    return P;
  }
  //compute likelihood via Felsenstein's Pruning Algorithm
  real pruning(int N, int B, int[] child, int[] parent, real[] brlen, matrix tiplik, real alpha, real beta) {
    matrix[N,2] lambda;                  //likelihoods at tips+nodes
    vector[2] pi;                         //stationary probability
    lambda = log(tiplik);
    for (b in 1:B) {
      matrix[2,2] P = evprob(brlen[b], alpha, beta); //via matrix exponentiation
      for (d in 1:2) {
        lambda[parent[b],d] += log(dot_product(P[d],exp(lambda[child[b]])));
      }
    }
    pi[1] = log(beta) - log(alpha+beta) + lambda[parent[B],1];
    pi[2] = log(alpha) - log(alpha+beta) + lambda[parent[B],2];
    return(log_sum_exp(pi));
  }
}
data {
  int<lower=1> T;                       //number of tips
  int<lower=1> N;                       //number of total nodes (incl. tips)
  int<lower=1> B;                       //number of branches
  int<lower=1> D;                       //number of features (i.e., lemma-pattern pairs)
  int<lower=1> L;                       //number of lemmas
  int<lower=1> J;                       //number of patterns
  int<lower=1> child[B];                //child of each branch
  int<lower=1> parent[B];               //parent of each branch
  real<lower=0> brlen[B];               //length of each branch
  matrix[N,D*2] tiplik;                 //likelihoods for data at tips in tree
  int<lower=1,upper=L> lemma_id[D];
  int<lower=1,upper=J> pattern_id[D];
}
parameters {
  real<lower=0> rho;
  real alpha_s;
  real alpha_p;
  real beta_s_lemma;
  real beta_p_lemma;
  simplex[L-1] simp_s_lemma;
  simplex[L-1] simp_p_lemma;
  vector[J-1] beta_s_pattern;
  vector[J-1] beta_p_pattern;
  real<lower=0> sigma_s;
  real<lower=0> sigma_p;
  real eps_s[D];
  real eps_p[D];
}
transformed parameters {
  real s[D];
  real p[D];
  vector[L] mo_s_lemma = rep_vector(0,L);
  vector[J] mo_s_pattern = rep_vector(0,J);
  vector[L] mo_p_lemma = rep_vector(0,L);
  vector[J] mo_p_pattern = rep_vector(0,J);
  mo_s_lemma[2:L] = cumulative_sum(simp_s_lemma)*beta_s_lemma;
  mo_s_pattern[2:J] = cumulative_sum(beta_s_pattern);
  mo_p_lemma[2:L] = cumulative_sum(simp_p_lemma)*beta_p_lemma;
  mo_p_pattern[2:J] = cumulative_sum(beta_p_pattern);
  for (d in 1:D) {
    s[d] = inv_logit(alpha_s + mo_s_lemma[lemma_id[d]] + mo_s_pattern[pattern_id[d]] + eps_s[d]*sigma_s)*rho;
    p[d] = inv_logit(alpha_p + mo_p_lemma[lemma_id[d]] + mo_p_pattern[pattern_id[d]] + eps_p[d]*sigma_p);
  }
}
model {
  rho ~ uniform(0,10);
  alpha_s ~ normal(0,1);
  alpha_p ~ normal(0,1);
  beta_s_lemma ~ normal(0,1);
  beta_p_lemma ~ normal(0,1);
  simp_s_lemma ~ dirichlet(rep_vector(1,L-1));
  simp_p_lemma ~ dirichlet(rep_vector(1,L-1));
  beta_s_pattern ~ normal(0,1);
  beta_p_pattern ~ normal(0,1);
  sigma_s ~ normal(0,1);
  sigma_p ~ normal(0,1);
  eps_s ~ normal(0,1);
  eps_p ~ normal(0,1);
  for (d in 1:D) {
    target += pruning(N,B,child,parent,brlen,tiplik[,((2*d)-1):(2*d)],p[d]*s[d],(1-p[d])*s[d]);
  }
}

