# Calculate Bayesian Multi-Chain MCMC Double Observer Estimate

Uses a multi-chain Collapsed Poisson-Gamma-Beta Gibbs Sampler (Data
Augmentation MCMC) to sample the joint posterior distribution of latent
undetected groups \\f\_{0,i}\\ and observer detection probabilities
across sectors. Supports both \`"independent"\` (two independent
observers) and \`"dependent"\` (2-pass primary/secondary removal)
protocols.

## Usage

``` r
bayesian_conjugate_estimate(
  both,
  obs1_only,
  obs2_only,
  obs_type = c("independent", "dependent"),
  iter = 2000,
  chains = 4,
  burnin = NULL,
  prior_alpha = 1,
  prior_beta = 0.5
)
```

## Arguments

- both:

  Integer vector of numbers of groups seen by both observers

- obs1_only:

  Integer vector of numbers of groups seen by Observer 1 only

- obs2_only:

  Integer vector of numbers of groups seen by Observer 2 only

- obs_type:

  Character string specifying observer protocol: \`"independent"\`
  (default) or \`"dependent"\`.

- iter:

  Integer number of total MCMC iterations per chain (default 2000)

- chains:

  Integer number of independent Markov chains (default 4)

- burnin:

  Integer number of warm-up/burn-in iterations per chain (default
  \`floor(iter / 4)\`)

- prior_alpha:

  Numeric Beta hyperprior shape1 parameter for detection probability
  (default 1)

- prior_beta:

  Numeric Beta hyperprior shape2 parameter for detection probability
  (default 0.5)

## Value

A list containing the posterior \`mean\`, \`median\`, \`var\`, 95
vectors (\`ci_lower\`, \`ci_upper\`), multi-chain convergence
diagnostics (\`rhat\`, \`ess\`), and full posterior sample arrays
(\`G_draws\`, \`p1_draws\`, \`mu_p1_draws\`).
