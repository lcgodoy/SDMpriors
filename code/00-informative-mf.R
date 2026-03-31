library(ggplot2)
library(cmdstanr)
library(bayesplot)

set.seed(42)

## Define the sampling design
N <- 100
## Sample temperatures uniformly between 5 and 35 degrees Celsius
x <- sort(runif(N, 5, 35)) 

## Define the TRUE biological parameters
true_T_opt <- 20
true_alpha_max <- 2.0 ## Logit of 2 is ~88% peak probability of occurrence

## If probability approaches 0 (logit ~ -4.6) at T = 10 (10 degrees away from T_opt):
## -4.6 = 2.0 - beta_width * (10 - 20)^2  ==>  beta_width = 0.066
true_beta_width <- 0.066

## Calculate the true biological mean function on the logit scale
true_mu <- true_alpha_max - true_beta_width * (x - true_T_opt)^2

## Simulate the Latent Gaussian Process (adding environmental noise)
true_alpha_gp <- 1.0 ## GP variance
true_rho_gp <- 4.0   ## GP length-scale (how fast the noise wiggles)

## Create the covariance matrix using the Squared Exponential kernel
dist_matrix <- as.matrix(dist(x))
K <- (true_alpha_gp^2) * exp(-0.5 * (dist_matrix / true_rho_gp)^2)
diag(K) <- diag(K) + 1e-10 # Jitter for numerical stability
L <- chol(K)

## Draw the true latent GP values (Mean function + GP noise)
set.seed(2026)
raw_f <- rnorm(N)
latent_f <- crossprod(L, raw_f) + true_mu 

# 4. Convert to probabilities and simulate presence/absence
p_occurrence <- 1 / (1 + exp(-latent_f)) # Inverse logit
y <- rbinom(N, size = 1, prob = p_occurrence)
plot(x, p_occurrence, type = "l")
points(x, y, pch = 19)

## Bundle the data for Stan
stan_data <- list(
  N = N,
  x = x,
  y = y,
  lit_T_opt = 20.5, ## what the literature says is the optimal temperature  
  lit_T_opt_sd = 1.5 ## How confident we are on the literature
)

## Compile the model
mod <- cmdstan_model("src/gp_inf_mean.stan")

## Sample from the posterior
fit <- mod$sample(
  data = stan_data,
  seed = 123,
  ## adapt_delta = 0.95 ## High adapt_delta helps GPs sample without divergent transitions
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000  
)

## Print a summary of the key biological parameters
fit$summary(variables = c("T_opt", "beta_width",
                          "alpha_max", "rho", "alpha"))

fit$draws(variables = c("T_opt", "beta_width",
                        "alpha_max", "rho", "alpha")) |>
  mcmc_combo(combo = c("dens_overlay", "trace"))

## Extract posterior samples of the latent GP (f)
posterior_f <- fit$draws("f", format = "matrix")

## Convert logit 'f' to probabilities
posterior_p <- 1 / (1 + exp(-posterior_f))

## Calculate median and 95% Credible Intervals for each temperature point
p_est <- apply(posterior_p, 2, median)
p_lower <- apply(posterior_p, 2, quantile, probs = 0.025)
p_upper <- apply(posterior_p, 2, quantile, probs = 0.975)

## Create a dataframe for plotting
plot_df <- data.frame(
  Temperature = x,
  Observed = y,
  True_Prob = p_occurrence,
  Est_Prob = p_est,
  Lower_CI = p_lower,
  Upper_CI = p_upper
)

## Plot everything together using ggplot2
ggplot(plot_df, aes(x = Temperature)) +
  ## 95% Credible Interval from the GP
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.3, fill = "blue") +
  ## Estimated mean probability from the GP
  geom_line(aes(y = Est_Prob), color = "blue", size = 1.2, linetype = "solid") +
  ## The true generating probability (what we want the GP to find!)
  geom_line(aes(y = True_Prob), color = "red", size = 1, linetype = "dashed") +
  ## The actual observed 0s and 1s
  geom_point(aes(y = Observed), alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(
    title = "Gaussian Process Recovery of Thermal Response",
    subtitle = "Blue = GP Posterior (Mean ± 95% CI), Red Dashed = True Underlying Probability",
    x = "Temperature (°C)",
    y = "Probability of Occurrence"
  )

##--- option 2 ----

## Assuming you have your true field data in vectors `x_field` and `y_field`
## Let's say the literature says T_min is ~10°C and T_max is ~30°C

## 1. Define the pseudo-data (building "walls" of zeros at the limits)
## We add points exactly at the limits, and a few degrees beyond them
x_pseudo <- c(8, 9, 10,   ## T_min wall
              30, 31, 32) ## T_max wall

y_pseudo <- rep(0, length(x_pseudo)) ## All are absences (0)

## 2. Combine your field data with the pseudo-data
x_combined <- c(x, x_pseudo)
y_combined <- c(y, y_pseudo)

## 3. Bundle the data for Stan
## Define the weight for your experimental pseudo-data
## A weight of 0.2 means 5 pseudo-points equal the "evidence" of 1 real point.
pseudo_weight <- 0.2 

## 1. Create the weights vectors
weights_field <- rep(1, length(x))
weights_pseudo <- rep(pseudo_weight, length(x_pseudo))

## 2. Combine them exactly like the x and y data
weights_combined <- c(weights_field, weights_pseudo)

## 3. Add to your Stan data list
stan_data_pseudo <- list(
  N = length(x_combined),
  x = x_combined,
  y = y_combined,
  weights = weights_combined  # Tell Stan who is who!
)

## Compile the model
mod_naive <- cmdstan_model("src/gp_augmented_naive.stan")

## Sample from the posterior
fit_naive <- mod_naive$sample(
  data = stan_data_pseudo,
  seed = 123,
  ## adapt_delta = 0.95 ## High adapt_delta helps GPs sample without divergent transitions
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000  
)

## Print a summary of the key biological parameters
fit_naive$summary(variables = c("rho", "alpha"))

fit_naive$draws(variables = c("rho", "alpha")) |>
  mcmc_combo(combo = c("dens_overlay", "trace"))

## Extract posterior samples of the latent GP (f)
posterior_f2 <- fit_naive$draws("f", format = "matrix")

## Convert logit 'f' to probabilities
posterior_p2 <- 1 / (1 + exp(-posterior_f2))

## Calculate median and 95% Credible Intervals for each temperature point
p_est2 <- apply(posterior_p2, 2, median)
p_lower2 <- apply(posterior_p2, 2, quantile, probs = 0.025)
p_upper2 <- apply(posterior_p2, 2, quantile, probs = 0.975)

## Create a dataframe for plotting
plot_df <- data.frame(
  Temperature = x,
  Observed = y,
  True_Prob = p_occurrence,
  Est_Prob = p_est2[seq_len(N)],
  Lower_CI = p_lower2[seq_len(N)],
  Upper_CI = p_upper2[seq_len(N)]
)

## Plot everything together using ggplot2
ggplot(plot_df, aes(x = Temperature)) +
  ## 95% Credible Interval from the GP
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.3, fill = "blue") +
  ## Estimated mean probability from the GP
  geom_line(aes(y = Est_Prob), color = "blue", size = 1.2, linetype = "solid") +
  ## The true generating probability (what we want the GP to find!)
  geom_line(aes(y = True_Prob), color = "red", size = 1, linetype = "dashed") +
  ## The actual observed 0s and 1s
  geom_point(aes(y = Observed), alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(
    title = "Gaussian Process Recovery of Thermal Response",
    subtitle = "Blue = GP Posterior (Mean ± 95% CI), Red Dashed = True Underlying Probability",
    x = "Temperature (°C)",
    y = "Probability of Occurrence"
  )

##--- option 3: conditional GPs ----

## 1. Define the Experimental (Literature) Data
## T_min = 10, T_max = 30
x_exp <- c(10.0, 30.0) 

## Set the exact latent logit values for these points
## logit(-5) is ~0.0067 probability (practically zero)
f_exp <- c(-5.0, -5.0) 

## 2. Simulate True Field Observational Data
N_obs <- 100
## Sort the temperatures so our line plots look clean later
x_obs <- sort(runif(N_obs, 12, 28)) 

## Simulate the true biological curve (peak at 20C)
true_T_opt <- 20
true_alpha_max <- 2.0 
true_beta_width <- 0.07
true_logit <- true_alpha_max - true_beta_width * (x_obs - true_T_opt)^2

## Convert to probabilities and simulate presence/absence
true_prob <- 1 / (1 + exp(-true_logit))
y_obs <- rbinom(N_obs, size = 1, prob = true_prob)

## 3. Bundle the data for the Exact Conditioned Stan Model
stan_data_cond <- list(
  N_obs = N_obs,
  x_obs = x_obs,
  y_obs = y_obs,
  N_exp = length(x_exp),
  x_exp = x_exp,
  f_exp = f_exp
)

## Compile the conditioned GP model
mod_cond <- cmdstan_model("src/gp_aug_conditioning.stan")

## Sample from the posterior
fit_cond <- mod_cond$sample(
  data = stan_data_cond,
  seed = 123,
  ## adapt_delta = 0.95,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

## Check the hyperparameters
fit_cond$summary(variables = c("rho", "alpha"))

## Extract the posterior of the latent GP at the observation points
posterior_f_obs <- fit_cond$draws("f_obs", format = "matrix")

## Convert from logit scale to probability scale (inverse logit)
posterior_p_obs <- 1 / (1 + exp(-posterior_f_obs))

## Calculate the median and 95% Credible Intervals
p_est <- apply(posterior_p_obs, 2, median)
p_lower <- apply(posterior_p_obs, 2, quantile, probs = 0.025)
p_upper <- apply(posterior_p_obs, 2, quantile, probs = 0.975)

## Create a dataframe for the observational predictions
plot_df <- data.frame(
  Temperature = x_obs,
  Observed_Y = y_obs,
  True_Prob = true_prob,
  Est_Prob = p_est,
  Lower_CI = p_lower,
  Upper_CI = p_upper
)

## Create a small dataframe for the experimental anchors to show on the plot
exp_prob <- 1 / (1 + exp(-f_exp)) ## Convert -5.0 logit to probability
anchor_df <- data.frame(
  Temperature = x_exp,
  Prob = exp_prob
)

## Plotting with ggplot2
ggplot(plot_df, aes(x = Temperature)) +
  ## 95% CI Ribbon
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.3, fill = "purple") +
  ## Estimated mean probability
  geom_line(aes(y = Est_Prob), color = "purple", size = 1.2) +
  ## True underlying probability (for comparison)
  geom_line(aes(y = True_Prob), color = "black", linetype = "dashed", size = 0.8) +
  ## Actual observed field data (0s and 1s)
  geom_point(aes(y = Observed_Y), alpha = 0.5, size = 2) +
  
  ## --- ADD THE EXPERIMENTAL ANCHORS ---
  geom_point(data = anchor_df, aes(x = Temperature, y = Prob), 
             color = "red", size = 4, shape = 4, stroke = 2) +
  geom_vline(xintercept = x_exp, linetype = "dotted", color = "red") +
  
  theme_minimal() +
  coord_cartesian(xlim = c(8, 32)) + ## Zoom out slightly to see the anchors
  labs(
    title = "Exact GP Conditioning on Experimental Limits",
    subtitle = "Purple = GP Posterior, Red 'X' = Literature Limits (Exact Anchors)",
    x = "Temperature (°C)",
    y = "Probability of Occurrence"
  )
