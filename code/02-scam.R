library(ggplot2)
library(cmdstanr)
library(bayesplot)
library(terra)
library(dplyr)
library(cowplot)
library(scam)

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

error_invgamma <- function(comb, CTmin, CTmax){
  (abs(invgamma::pinvgamma(CTmin, comb[1], comb[2]) - 0.05) +
   abs(invgamma::pinvgamma(CTmax, comb[1], comb[2]) - 0.95)) / 2
}

param_irvgamma <- function(CTmin, CTmax){
  combs <- expand.grid(shape = seq(1, 50, 1),
                       rate = seq(1, 50, 1))
  error <- apply(combs,
                 MARGIN = 1,
                 error_invgamma,
                 CTmin = CTmin,
                 CTmax = CTmax)
  return(combs[which.min(error),])
}


params <- param_irvgamma(CTmin1, CTmax1)
s <- expm1(seq(log1p(0), log1p(50), len = 500))
plot(s, invgamma::dinvgamma(s, params$shape, params$rate), type = 'l')

## We can think of a skew normal too, since, for this species, the thermal range
## seems to be skewed

##---- SCAM model ----

aux_x <- splines2::cSpline(pa$trmax, knots = 10, intercept = FALSE)

aux_x2 <- splines2::cSpline(s, knots = 10, intercept = FALSE)

beta_init <- rnorm(NCOL(aux_x2))

objective <- invgamma::dinvgamma(s, params$shape, params$rate)
 
to_opt <- function(theta, x, y) {
  y_pred <- c(x %*% theta)
  y_pred <- y_pred / sum(y_pred)
  out <- y_pred - (y / sum(y))
  sum(out * out) |>
    log()
}

prior_betas <-
  optim(par = beta_init, fn = to_opt, method = "BFGS",
        x = aux_x2, y = objective)

plot(x = s, y = aux_x2 %*% prior_betas$par, type = "l")

## 3. Bundle the data for the Exact Conditioned Stan Model
stan_data_cond <- list(
  N_obs = NROW(pa),
  X = aux_x,
  K = NCOL(aux_x),
  y_obs = pa$pres,
  pri_beta = prior_betas$par
)

## Compile the conditioned GP model
mod_cond <- cmdstan_model("src/scam.stan")

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
fit_cond$summary(variables = c("alpha", "beta", "sigma"))

fit_cond$draws(variables = c("alpha", "beta", "sigma")) |>
  mcmc_combo(combo = c("dens_overlay", "trace"),
             facet_args = list(labeller = label_parsed))


## Extract the posterior of the latent GP at the observation points
posterior_f_obs <- fit_cond$draws("logit_mu", format = "matrix")

## Convert from logit scale to probability scale (inverse logit)
posterior_p_obs <- 1 / (1 + exp(-posterior_f_obs))

## Calculate the median and 95% Credible Intervals
p_est <- apply(posterior_p_obs, 2, median)
p_lower <- apply(posterior_p_obs, 2, quantile, probs = 0.025)
p_upper <- apply(posterior_p_obs, 2, quantile, probs = 0.975)

## Create a dataframe for the observational predictions
plot_df <- data.frame(
  Temperature = pa$trmax,
  Observed_Y = stan_data_cond$y_obs,
  Est_Prob = p_est,
  Lower_CI = p_lower,
  Upper_CI = p_upper
)

##--- predictions ----

pred_x <- s
beta_draws <- fit_cond$draws("beta", format = "matrix")
alpha_draws <- fit_cond$draws("alpha", format = "matrix")

preds <- tcrossprod(beta_draws, aux_x2)
preds <- apply(preds, 2, \(x, a) x + c(a),
               a = alpha_draws)
preds <- 1 / (1 + exp(-preds))

pred_dt <- posterior::summarise_draws(preds)

pred_df <-
  pred_dt |>
  mutate(Temperature = pred_x)

## Plotting with ggplot2
p1 <-
  ggplot(pred_df, aes(x = Temperature)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.3, fill = "purple") +
  geom_line(aes(y = median), color = "purple", size = 1.2) +
  geom_point(data = plot_df,
             aes(x = Temperature, y = Observed_Y), 
             size = 4, color = scales::alpha(1, .5)) +
  theme_bw() +
  labs(
      title = "GP Conditioning on Experimental Limits",
      x = "Temperature (°C)",
      y = "Probability of Occurrence"
  )

p1
