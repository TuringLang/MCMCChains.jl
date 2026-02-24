data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> y;
  int<lower=0> K;
  matrix[N, K] X;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
  array[2] real gamma;
}
model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 2);
  sigma ~ exponential(1);
  gamma ~ normal(0, 1);
  y ~ bernoulli_logit(alpha + X * beta);
}
generated quantities {
  matrix[2, 3] mat_param;
  for (i in 1:2)
    for (j in 1:3)
      mat_param[i, j] = normal_rng(0, 1);
}
