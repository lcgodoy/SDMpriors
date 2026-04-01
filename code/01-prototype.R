library(ggplot2)
library(cmdstanr)
library(bayesplot)
library(terra)
library(dplyr)
library(cowplot)

species_list <- read.csv("out/presabs/SpeciesList_PresAbs.csv")

##--- clim data ----
##use july for max
temp <- rast("data/microclim/0_shade/TA1cm_soil_0_7.nc")
tmax_0 <- mean(temp) ##or max

##use july for max
temp <- rast("data/microclim/50_shade/TA1cm_soil_50_7.nc")
tmax_50 <- mean(temp)

##use july for max
temp <- rast("data/microclim/100_shade/TA1cm_soil_100_7.nc")
tmax_100 <- mean(temp)

##--- species presence ----

##---* choosing the species with least data ----
species_n <-
  species_list |>
  pull(spec) |>
  sapply(\(sp) {
    read.csv(sprintf("out/presabs/PresAbs_%s.csv",
                     janitor::make_clean_names(sp))) |>
      NROW()
  })

sp_id <- 1

pa <- species_list$spec[sp_id] |>
  janitor::make_clean_names() |>
  sprintf("out/presabs/PresAbs_%s.csv", ..2 = _) |>
  read.csv()

##----
##plot localities and temperature
##crop to observed range 
my_ext <- ext(c(range(pa$lon), range(pa$lat))) ## define the extent
## extent

my_ext[1] <- my_ext[1] - 10
my_ext[2] <- my_ext[2] + 10
my_ext[3] <- my_ext[3] - 10
my_ext[4] <- my_ext[4] + 10

##crop
tmax0 <- crop(tmax_0, my_ext)
tmax50 <- crop(tmax_50, my_ext)
tmax100 <-  crop(tmax_100, my_ext)

##--- model thermoregulation ----
CTmin1 <- species_list$tmin[sp_id]
CTmax1 <- species_list$tmax[sp_id]
## approximate Topt, but fix based on data
Topt <- CTmin1 + (CTmax1 - CTmin1) * 0.7

##--- setting up prior ----

## Normal distribution such that P(CTmin1 < x < CTmax1) = 0.99 and its mean is
## Topt:
prb <- 0.99
mu <- Topt
p_upper <- 1 - ((1 - prb) / 2)
z <- qnorm(p_upper)
sigma <- (CTmax1 - mu) / z

curve(dnorm(x, mu, sigma),
      from = CTmin1 - 2 * sigma,
      to = CTmax1 + 2 * sigma)
abline(v = CTmin1, lty = 2)
abline(v = CTmax1, lty = 2)
abline(v = Topt, lty = 2)

## same, but mean as the mid-point between CTmin1 and CTmax1
prb <- 0.99
mu <- (CTmin1 + CTmax1) / 2
p_upper <- 1 - ((1 - prb) / 2)
z <- qnorm(p_upper)
sigma <- (CTmax1 - mu) / z

curve(dnorm(x, mu, sigma),
      from = CTmin1 - 2 * sigma,
      to = CTmax1 + 2 * sigma)
abline(v = CTmin1, lty = 2)
abline(v = CTmax1, lty = 2)
abline(v = Topt, lty = 2)

## We can think of a skew normal too, since, for this species, the thermal range
## seems to be skewed

##---- GP model ----

## Define the Experimental (Literature) Data
x_exp <- seq(from = CTmin1 - 2 * sigma,
             to = CTmax1 + 2 * sigma,
             length.out = 5)
f_exp <- qlogis(dnorm(x_exp, mu, sigma))
## all these density values are pretty small, so I am setting it such that they
## lie between -5 and 5

lb <- -5
ub <- 5
f_exp <- (f_exp - min(f_exp)) / diff(range(f_exp))
f_exp <- f_exp * (ub - lb) + lb

## Set the exact latent logit values for these points
## logit(-5) is ~0.0067 probability (practically zero)
## logit(5) is ~0.9933 probability (practically one)


## 2. Field Observational Data

## 3. Bundle the data for the Exact Conditioned Stan Model
stan_data_cond <- list(
  N_obs = NROW(pa),
  x_obs = pa$trmax,
  y_obs = pa$pres,
  N_exp = length(x_exp),
  x_exp = x_exp,
  f_exp = f_exp,
  rho_thresh = sd(pa$trmax),
  p_rho = 0.4
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
fit_cond$summary(variables = c("rho", "alpha", "kappa"))

fit_cond$draws(variables = c("rho", "alpha", "kappa")) |>
  mcmc_combo(combo = c("dens_overlay", "trace"),
             facet_args = list(labeller = label_parsed))


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
  Temperature = stan_data_cond$x_obs,
  Observed_Y = stan_data_cond$y_obs,
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

##--- predictions ----

mod_pred <- cmdstan_model("src/gp_augc_pred.stan")

pred_x <- seq(from = min(x_exp), to = max(x_exp),
              length.out = 200)

dat_pred <- stan_data_cond
dat_pred$x_pred <- pred_x
dat_pred$N_pred <- NROW(pred_x)

my_preds <-
  mod_pred$generate_quantities(fitted_params =
                                 fit_cond$draws(variables = c("rho", "alpha", "kappa",
                                                              "f_obs")),
                               data = dat_pred,
                               seed = 2026)

pred_df <-
  my_preds$summary(variable = "y_pred_prob") |>
  mutate(Temperature = pred_x)

## Plotting with ggplot2
p1 <-
  ggplot(pred_df, aes(x = Temperature)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.3, fill = "purple") +
  geom_line(aes(y = median), color = "purple", size = 1.2) +
  geom_point(data = anchor_df,
             aes(x = Temperature, y = Prob), 
             color = "red", size = 4,
             shape = 4, stroke = 2) +
  geom_point(data = plot_df,
             aes(x = Temperature, y = Observed_Y), 
             size = 4, color = scales::alpha(1, .5)) +
  geom_vline(data = anchor_df,
             aes(xintercept = Temperature),
             linetype = "dotted",
             color = "red") +
  theme_bw() +
  labs(
      title = "GP Conditioning on Experimental Limits",
      x = "Temperature (°C)",
      y = "Probability of Occurrence"
  )

##--- No conditioning ----

## Compile the conditioned GP model
mod_uncond <- cmdstan_model("src/gp_vanilla.stan")

## Sample from the posterior
fit_uncond <- mod_uncond$sample(
  data = stan_data_cond,
  seed = 123,
  ## adapt_delta = 0.95,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

## Check the hyperparameters
fit_uncond$summary(variables = c("rho", "alpha"))

fit_cond$draws(variables = c("rho", "alpha")) |>
  mcmc_combo(combo = c("dens_overlay", "trace"),
             facet_args = list(labeller = label_parsed))


## Extract the posterior of the latent GP at the observation points
posterior_f_obs <- fit_uncond$draws("f_obs", format = "matrix")

## Convert from logit scale to probability scale (inverse logit)
posterior_p_obs <- 1 / (1 + exp(-posterior_f_obs))

## Calculate the median and 95% Credible Intervals
p_est <- apply(posterior_p_obs, 2, median)
p_lower <- apply(posterior_p_obs, 2, quantile, probs = 0.025)
p_upper <- apply(posterior_p_obs, 2, quantile, probs = 0.975)

## Create a dataframe for the observational predictions
plot_df <- data.frame(
  Temperature = stan_data_cond$x_obs,
  Observed_Y = stan_data_cond$y_obs,
  Est_Prob = p_est,
  Lower_CI = p_lower,
  Upper_CI = p_upper
)

##--- predictions ----

mod_pred2 <- cmdstan_model("src/gp_van_pred.stan")

pred_x <- seq(from = min(x_exp), to = max(x_exp),
              length.out = 200)

dat_pred <- stan_data_cond
dat_pred$x_pred <- pred_x
dat_pred$N_pred <- NROW(pred_x)

my_preds2 <-
  mod_pred2$generate_quantities(fitted_params =
                                  fit_uncond$draws(variables = c("rho", "alpha",
                                                                 "f_obs")),
                                data = dat_pred,
                                seed = 2026)

pred_df2 <-
  my_preds2$summary(variable = "y_pred_prob") |>
  mutate(Temperature = pred_x)

## Plotting with ggplot2
p2 <-
  ggplot(pred_df2, aes(x = Temperature)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.3, fill = "purple") +
  geom_line(aes(y = median), color = "purple", size = 1.2) +
  geom_point(data = plot_df,
             aes(x = Temperature, y = Observed_Y), 
             size = 4, color = scales::alpha(1, .5)) +
  geom_vline(data = anchor_df,
             aes(xintercept = Temperature),
             linetype = "dotted",
             color = "red") +
  theme_bw() +
  labs(
      title = "GP Vanilla",
      x = "Temperature (°C)",
      y = "Probability of Occurrence"
  )


plot_grid(p1, p2)
