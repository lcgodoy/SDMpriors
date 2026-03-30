library(ggplot2)
library(cmdstanr)
library(dplyr)

set.seed(8675309)

## 1. Simulate the Field Observations
N_obs <- 80
x_obs <- sort(runif(N_obs, 12, 28)) 

## True biological curve (peak at 20C)
true_T_opt <- 20
true_alpha_max <- 2.0 
true_beta_width <- 0.07
true_logit <- true_alpha_max - true_beta_width * (x_obs - true_T_opt)^2
true_prob <- 1 / (1 + exp(-true_logit))
y_obs <- rbinom(N_obs, size = 1, prob = true_prob)

## 2. Define the Grids
## A continuous grid of virtual points to enforce the shape
x_deriv_grid <- seq(10, 30, length.out = 20) 
## A high-res grid for smooth plotting
x_pred_grid <- seq(8, 32, length.out = 100)  

## 3. Bundle Data for Stan
stan_data <- list(
    N_obs = N_obs,
    x_obs = x_obs,
    y_obs = y_obs,
    N_deriv = length(x_deriv_grid),
    x_deriv = x_deriv_grid,
    nu = 0.01,     ## Strictness of the derivative constraint
    kappa = 2.0,   ## Steepness of the tanh transition    
    N_pred = length(x_pred_grid),
    x_pred = x_pred_grid,
    jitter = 1e-10
)

## 4. Compile and Fit the Model
mod <- cmdstan_model("src/gp_diff.stan")

fit <- mod$sample(
               data = stan_data,
               seed = 42,
               ## adapt_delta = 0.99, ## High adapt_delta needed for strict constraints
               ## max_treedepth = 12,
               chains = 4,
               parallel_chains = 4,
               iter_warmup = 1000,
               iter_sampling = 1000,
           )

## Check the estimated Optimum!
fit$summary(variables = c("T_opt", "rho", "alpha"))
fit$draws(variables = c("T_opt", "rho", "alpha")) |>
  bayesplot::mcmc_combo(combo = c("dens_overlay", "trace"))

## --- 5. Generate Smooth Predictions ---
mod_pred <- cmdstan_model("src/gp_diff_pred.stan")

## We append the prediction grid to the data list
stan_data_pred <- stan_data
stan_data_pred$N_pred <- length(x_pred_grid)
stan_data_pred$x_pred <- x_pred_grid

## This bypasses HMC entirely and just runs the generated quantities block
## using the posterior draws from fit_base. It should take ~2 seconds.
fitted_pars <- fit$draws(variables = c("rho", "alpha", "baseline_logit",
                                       "T_opt", "eta_joint"))

fit_pred <-
  mod_pred$generate_quantities(fitted_params = fitted_pars[, -4, ],
                               data = stan_data_pred,
                               seed = 42)

## --- 3. Extract ---
## You now extract 'p_pred' from the predictive fit object!
posterior_p_pred <- fit_pred$draws("p_pred", format = "matrix")

## Calculate medians and 95% intervals
pred_est <- apply(posterior_p_pred, 2, median)
pred_lower <- apply(posterior_p_pred, 2, quantile, probs = 0.025)
pred_upper <- apply(posterior_p_pred, 2, quantile, probs = 0.975)

## Plotting dataframe for smooth curve
smooth_df <- data.frame(
    Temperature = x_pred_grid,
    Est_Prob = pred_est,
    Lower_CI = pred_lower,
    Upper_CI = pred_upper
)

## Plotting dataframe for raw observations
obs_df <- data.frame(
    Temperature = x_obs,
    Observed_Y = y_obs
)

## Generate Final Plot
ggplot() +
  ## 95% CI Ribbon
  geom_ribbon(data = smooth_df, aes(x = Temperature, ymin = Lower_CI, ymax = Upper_CI), 
              alpha = 0.3, fill = "dodgerblue") +
  ## Estimated mean probability
  geom_line(data = smooth_df, aes(x = Temperature, y = Est_Prob), 
            color = "dodgerblue", size = 1.2) +
  ## Actual observed field data (0s and 1s)
  geom_point(data = obs_df, aes(x = Temperature, y = Observed_Y), 
             alpha = 0.5, size = 2) +
  
  theme_minimal() +
  labs(
      title = "Dynamic Shape-Constrained Thermal Niche",
      subtitle = paste0("Model Estimated T_opt: ", 
                        round(median(fit$draws("T_opt")), 1), " °C"),
      x = "Temperature (°C)",
      y = "Probability of Occurrence"
  )
