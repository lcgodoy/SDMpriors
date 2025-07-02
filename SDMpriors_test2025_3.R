#Gaussian Process Species Distribution Models with physiological priors, following the approach described in Golding & Purse (2016) but without relying on the GRaF package, https://github.com/goldingn/GRaF

#Gaussian process regression implementation: https://cran.r-project.org/web/packages/GauPro/vignettes/GauPro.html

#' GP-SDM Wrapper Class (Fixed)
#' 
#' Creates a wrapper around GauPro objects to store additional metadata
#' without modifying the locked GauPro environment.
GP_SDM <- R6::R6Class(
  "GP_SDM",
  public = list(
    #' @field gp_model The underlying GauPro model object
    gp_model = NULL,
    
    #' @field mean_function The physiological prior function
    mean_function = NULL,
    
    #' @field env_vars Character vector of environmental variable names
    env_vars = NULL,
    
    #' @field response_col Name of the response column
    response_col = NULL,
    
    #' @field prior_specs List of prior specifications
    prior_specs = NULL,
    
    #' @field training_data Training data used to fit the model
    training_data = NULL,
    
    #' @field combination_method Method for combining priors
    combination_method = NULL,
    
    #' @field kernel_type Type of kernel used
    kernel_type = NULL,
    
    #' @field n_training Number of training observations
    n_training = NULL,
    
    #' @field presence_rate Presence rate in training data
    presence_rate = NULL,
    
    #' Initialize the GP-SDM wrapper
    #'
    #' @param gp_model Fitted GauPro model
    #' @param mean_function Prior function
    #' @param env_vars Environmental variables
    #' @param response_col Response column name
    #' @param prior_specs Prior specifications
    #' @param training_data Training dataset
    #' @param combination_method Prior combination method
    #' @param kernel_type Kernel type used
    #' @param n_training Number of training observations
    #' @param presence_rate Presence rate
    initialize = function(gp_model, mean_function = NULL, env_vars = NULL,
                          response_col = NULL, prior_specs = NULL, 
                          training_data = NULL, combination_method = NULL,
                          kernel_type = NULL, n_training = NULL, 
                          presence_rate = NULL) {
      self$gp_model <- gp_model
      self$mean_function <- mean_function
      self$env_vars <- env_vars
      self$response_col <- response_col
      self$prior_specs <- prior_specs
      self$training_data <- training_data
      self$combination_method <- combination_method
      self$kernel_type <- kernel_type
      self$n_training <- n_training
      self$presence_rate <- presence_rate
    },
    
    #' Predict method for GP-SDM (Fixed)
    #'
    #' @param newdata Data frame with environmental variables
    #' @param include_prior Whether to include physiological prior
    #' @param return_se Whether to return standard errors
    #' @param type Either "response" or "link"
    #' @return Data frame with predictions
    predict = function(newdata = NULL, include_prior = TRUE, 
                       return_se = TRUE, type = "response") {
      
      # Use training data if newdata not provided
      if (is.null(newdata)) {
        newdata <- self$training_data[, self$env_vars, drop = FALSE]
        message("No newdata provided - using training data for predictions")
      }
      
      # Validate newdata structure
      if (!is.data.frame(newdata)) {
        stop("newdata must be a data frame")
      }
      
      # Check for required environmental variables
      missing_vars <- setdiff(self$env_vars, colnames(newdata))
      if (length(missing_vars) > 0) {
        stop(paste("Missing environmental variables in newdata:", 
                   paste(missing_vars, collapse = ", ")))
      }
      
      # Extract environmental variables in correct order
      X_new <- newdata[, self$env_vars, drop = FALSE]
      n_pred <- nrow(X_new)
      
      # Handle missing values by replacing with training means
      for (col in colnames(X_new)) {
        col_missing <- is.na(X_new[[col]])
        if (any(col_missing)) {
          col_mean <- mean(self$training_data[[col]], na.rm = TRUE)
          X_new[col_missing, col] <- col_mean
          if (sum(col_missing) > 0) {
            message(paste("Replaced", sum(col_missing), "missing values in", 
                          col, "with mean:", round(col_mean, 3)))
          }
        }
      }
      
      # Convert to matrix for GauPro
      X_new_matrix <- as.matrix(X_new)
      
      # Generate GP predictions with robust error handling
      message(paste("Generating predictions for", n_pred, "locations"))
      
      # Initialize prediction components
      gp_mean <- NULL
      gp_se <- NULL
      
      # Try different prediction approaches
      if (return_se) {
        # Approach 1: Try with standard errors
        tryCatch({
          pred_result <- self$gp_model$predict(X_new_matrix, se.fit = TRUE)
          
          if (is.list(pred_result)) {
            if ("mean" %in% names(pred_result)) {
              gp_mean <- pred_result$mean
            }
            if ("se.fit" %in% names(pred_result)) {
              gp_se <- pred_result$se.fit
            }
          } else {
            gp_mean <- as.numeric(pred_result)
          }
        }, error = function(e) {
          message("Failed to get predictions with SE, trying without SE")
        })
      }
      
      # Approach 2: Fallback without standard errors
      if (is.null(gp_mean)) {
        tryCatch({
          pred_result <- self$gp_model$predict(X_new_matrix, se.fit = FALSE)
          
          if (is.list(pred_result) && "mean" %in% names(pred_result)) {
            gp_mean <- pred_result$mean
          } else {
            gp_mean <- as.numeric(pred_result)
          }
        }, error = function(e) {
          stop(paste("All GP prediction methods failed:", e$message))
        })
      }
      
      # Validate prediction length
      if (length(gp_mean) != n_pred) {
        stop(paste("GP prediction length mismatch: expected", n_pred, 
                   "got", length(gp_mean)))
      }
      
      # Validate SE length if available
      if (!is.null(gp_se) && length(gp_se) != n_pred) {
        warning(paste("SE length mismatch: expected", n_pred, 
                      "got", length(gp_se), "- ignoring SE"))
        gp_se <- NULL
      }
      
      # Initialize adjusted mean with GP predictions
      adjusted_mean <- gp_mean
      prior_contribution <- NULL
      
      # Add physiological prior if requested and available
      if (include_prior && !is.null(self$mean_function)) {
        message("Including physiological prior in predictions")
        
        tryCatch({
          # Calculate prior values
          prior_values <- self$mean_function(X_new)
          
          # Validate prior values length
          if (length(prior_values) != n_pred) {
            warning(paste("Prior function length mismatch: expected", n_pred, 
                          "got", length(prior_values), "- skipping prior"))
          } else {
            # Convert prior to logit scale to match GP scale
            prior_clamped <- pmax(pmin(prior_values, 0.999), 0.001)
            prior_logit <- qlogis(prior_clamped)
            
            # Add to GP mean (both on logit scale)
            adjusted_mean <- gp_mean + prior_logit
            prior_contribution <- prior_clamped
            
            # Check for numerical issues
            if (any(!is.finite(adjusted_mean))) {
              warning("Non-finite values in adjusted predictions - using GP predictions only")
              adjusted_mean <- gp_mean
              prior_contribution <- NULL
            }
          }
        }, error = function(e) {
          warning(paste("Error applying prior function:", e$message, "- using GP predictions only"))
          adjusted_mean <- gp_mean
        })
      }
      
      # Create result data frame with safe assignment
      if (type == "response") {
        # Transform to probability scale
        result <- data.frame(
          fit = plogis(adjusted_mean),
          stringsAsFactors = FALSE
        )
      } else {  # type == "link"
        # Keep on logit scale
        result <- data.frame(
          fit = adjusted_mean,
          stringsAsFactors = FALSE
        )
      }
      
      # Add standard errors and confidence intervals if available
      if (return_se) {
        if (!is.null(gp_se)) {
          result$se.fit <- gp_se
          
          if (type == "response") {
            # Calculate confidence intervals on logit scale, then transform
            lower_logit <- adjusted_mean - 1.96 * gp_se
            upper_logit <- adjusted_mean + 1.96 * gp_se
            result$lower <- plogis(lower_logit)
            result$upper <- plogis(upper_logit)
          }
        } else {
          # Add NA columns if SE requested but not available
          result$se.fit <- rep(NA_real_, n_pred)
          if (type == "response") {
            result$lower <- rep(NA_real_, n_pred)
            result$upper <- rep(NA_real_, n_pred)
          }
          message("Standard errors not available from GP model")
        }
      }
      
      # Add prior contribution if calculated
      if (!is.null(prior_contribution)) {
        result$prior_contribution <- prior_contribution
      }
      
      # Ensure result has correct dimensions
      if (nrow(result) != n_pred) {
        stop(paste("Result data frame has wrong number of rows:", 
                   nrow(result), "expected", n_pred))
      }
      
      message("Prediction completed successfully")
      return(result)
    },
    
    #' Print method for GP-SDM
    print = function() {
      cat("Gaussian Process Species Distribution Model\n")
      cat("==========================================\n\n")
      
      cat("Model Details:\n")
      cat(paste("  Kernel type:", self$kernel_type, "\n"))
      cat(paste("  Training observations:", self$n_training, "\n"))
      cat(paste("  Presence rate:", round(self$presence_rate, 3), "\n"))
      cat(paste("  Environmental variables:", length(self$env_vars), "\n"))
      
      if (length(self$env_vars) <= 10) {
        cat(paste("    -", paste(self$env_vars, collapse = ", "), "\n"))
      } else {
        cat(paste("    -", paste(self$env_vars[1:10], collapse = ", "), "... (truncated)\n"))
      }
      
      cat("\nPrior Information:\n")
      if (is.null(self$prior_specs)) {
        cat("  No physiological priors specified\n")
      } else {
        n_priors <- sum(!sapply(self$prior_specs, is.null))
        cat(paste("  Active priors:", n_priors, "of", length(self$prior_specs), "\n"))
        cat(paste("  Combination method:", self$combination_method, "\n"))
        
        for (i in seq_along(self$prior_specs)) {
          if (!is.null(self$prior_specs[[i]])) {
            spec <- self$prior_specs[[i]]
            cat(paste("    -", self$env_vars[i], ":", spec$type, "prior\n"))
          }
        }
      }
      
      cat("\nUse object$predict() to generate habitat suitability predictions\n")
      invisible(self)
    }
  )
)

#' S3 Predict Method for GP_SDM Objects
#'
#' @param object A fitted GP_SDM object
#' @param newdata Data frame for predictions
#' @param include_prior Whether to include physiological prior
#' @param return_se Whether to return standard errors
#' @param type Either "response" or "link"
#' @param ... Additional arguments
#'
#' @return Data frame with predictions
#'
#' @method predict GP_SDM
#' @export
predict.GP_SDM <- function(object, newdata = NULL, include_prior = TRUE, 
                           return_se = TRUE, type = "response", ...) {
  # Call the R6 object's predict method
  object$predict(newdata = newdata, include_prior = include_prior,
                 return_se = return_se, type = type)
}

#' S3 Print Method for GP_SDM Objects
#'
#' @param x A fitted GP_SDM object
#' @param ... Additional arguments
#'
#' @method print GP_SDM
#' @export
print.GP_SDM <- function(x, ...) {
  x$print()
  invisible(x)
}

#==================================
# Load required libraries
if (!requireNamespace("R6", quietly = TRUE)) {
  install.packages("R6")
}
library(R6)
library(GauPro)
library(dplyr)

# Create example data
set.seed(123)
n_sites <- 100

example_data <- data.frame(
  presence = rbinom(n_sites, 1, 0.3),
  temperature = runif(n_sites, 0, 40),
  precipitation = runif(n_sites, 0, 2000),
  elevation = runif(n_sites, 0, 3000)
)




#--------------------------------
# Define prior specifications
prior_specs <- list(
  # Temperature: Gaussian prior
  list(
    type = "gaussian", 
    optimal_value = 22, 
    sigma = 6, 
    lower_bound = 10, 
    upper_bound = 35
  ),
  # Precipitation: Polynomial prior
  list(
    type = "polynomial", 
    coefficients = c(0, 4, -4), 
    lower_bound = 500, 
    upper_bound = 1500, 
    normalize = TRUE
  ),
  # Elevation: Sigmoid prior
  list(
    type = "sigmoid", 
    midpoint = 1500, 
    steepness = 0.002, 
    lower_bound = 0, 
    upper_bound = 3000, 
    ascending = FALSE
  )
)

# Fit model with priors
cat("Fitting GP-SDM with priors...\n")
gp_model_with_priors <- fit_gp_sdm(
  data = example_data,
  response_col = "presence",
  env_vars = c("temperature", "precipitation", "elevation"),
  prior_specifications = prior_specs,
  combination_method = "multiplicative",
  kernel_type = "Matern52",
  seed = 123
)

# Print model summary
print(gp_model_with_priors)

# Make predictions using S3 method
cat("\nMaking predictions...\n")
predictions <- predict.GP_SDM(gp_model_with_priors, 
                       newdata = example_data[1:10, ],
                       include_prior = TRUE,
                       return_se = TRUE,
                       type = "response")

print(head(predictions))

# Alternative: Call R6 method directly
cat("\nAlternative prediction method...\n")
predictions_direct <- gp_model_with_priors$predict(
  newdata = example_data[1:10, ],
  include_prior = TRUE,
  return_se = TRUE,
  type = "response"
)

print(head(predictions_direct))

# Fit model without priors for comparison
cat("\nFitting GP-SDM without priors...\n")
gp_model_no_priors <- fit_gp_sdm(
  data = example_data,
  response_col = "presence",
  env_vars = c("temperature", "precipitation", "elevation"),
  prior_specifications = NULL,
  kernel_type = "Matern52",
  seed = 123
)

# Compare predictions
predictions_no_priors <- predict.GP_SDM(gp_model_no_priors, 
                                 newdata = example_data[1:10, ])

cat("\nComparison of first 5 predictions:\n")
comparison <- data.frame(
  with_priors = predictions$fit[1:5],
  without_priors = predictions_no_priors$fit[1:5],
  difference = predictions$fit[1:5] - predictions_no_priors$fit[1:5]
)
print(round(comparison, 3))

#=========================================
#Model testing

# Calculate threshold-independent metrics
  
  # Area Under the ROC Curve (AUC)
  auc_value <- as.numeric(pROC::auc(obs, pred, quiet = TRUE))
  
  # Calculate threshold-dependent metrics
  pred_binary <- ifelse(pred >= threshold, 1, 0)
  
  # Confusion matrix components
  tp <- sum(obs == 1 & pred_binary == 1)  # True positives
  tn <- sum(obs == 0 & pred_binary == 0)  # True negatives
  fp <- sum(obs == 0 & pred_binary == 1)  # False positives
  fn <- sum(obs == 1 & pred_binary == 0)  # False negatives
  
  # Basic metrics
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  sensitivity <- tp / (tp + fn)  # True positive rate
  specificity <- tn / (tn + fp)  # True negative rate
  precision <- tp / (tp + fp)    # Positive predictive value
  
  # Handle division by zero
  if (is.nan(sensitivity)) sensitivity <- 0
  if (is.nan(specificity)) specificity <- 0
  if (is.nan(precision)) precision <- 0
  
  # True Skill Statistic (TSS) = Sensitivity + Specificity - 1
  tss <- sensitivity + specificity - 1
  
  # Cohen's Kappa
  po <- accuracy  # Observed agreement
  pe <- ((tp + fn) * (tp + fp) + (tn + fp) * (tn + fn)) / (tp + tn + fp + fn)^2  # Expected agreement
  kappa <- (po - pe) / (1 - pe)
  if (is.nan(kappa)) kappa <- 0
  
  # Brier Score (lower is better)
  brier_score <- mean((pred - obs)^2)
  
  # Log-likelihood
  # Add small constant to avoid log(0)
  epsilon <- 1e-15
  pred_safe <- pmax(pmin(pred, 1 - epsilon), epsilon)
  log_likelihood <- sum(obs * log(pred_safe) + (1 - obs) * log(1 - pred_safe))
  
    # Confusion matrix
    confusion_matrix = matrix(c(tn, fp, fn, tp), nrow = 2, 
                              dimnames = list(Observed = c("0", "1"), 
                                              Predicted = c("0", "1")))
  
    #===========================
    #Example
    
    # load physiological priors from Sunday database
    if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/Shared Drives/TrEnCh/Projects/SDMpriors/")
    if(desktop=="n") setwd("/Users/lbuckley/Google Drive/Shared Drives/TrEnCh/Projects/SDMpriors/")
    
    dat= read.csv("out/presabs/SpeciesList_PresAbs.csv")
    
    #Load envi data
    #DATA CAN BE DOWNLOADED TO SHARED DRIVE FROM HERE: https://figshare.com/collections/microclim_Global_estimates_of_hourly_microclimate_based_on_long_term_monthly_climate_averages/878253
    #USED 1cm Air temperature
    
    if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/SDMpriors/")
    if(desktop=="n") setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/SDMpriors/")
    
    #load all clim data
    #use july for max
    temp= brick("data/microclim/0_shade/TA1cm_soil_0_7.nc")
    tmax_0= mean(temp) #or max
    names(tmax_0)<-"tmax0"
    
    #use july for max
    temp= brick("data/microclim/50_shade/TA1cm_soil_50_7.nc")
    tmax_50= mean(temp)
    names(tmax_50)<-"tmax50"
    
    #use july for max
    temp= brick("data/microclim/100_shade/TA1cm_soil_100_7.nc")
    tmax_100= mean(temp)
    names(tmax_100)<-"tmax100"
    
    #---------------
    #set up data storage
    #models<- array(NA, dim=c(nrow(dat), 3, 2), dimnames=list(dat$species, c("model","AUC","RMSE"),c("flat","prior") ))
    models<- matrix(NA, nrow=nrow(dat), ncol=6)
    rownames(models)<- dat$species
    models<- as.data.frame(models)
    colnames(models)<- c("mod1", "AUC1", "RMSE1", "mod2", "AUC2", "RMSE2")
    
    #load presence absence
    if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/Shared Drives/TrEnCh/Projects/SDMpriors/")
    if(desktop=="n") setwd("/Users/lbuckley/Google Drive/Shared Drives/TrEnCh/Projects/SDMpriors/")
    
    for(spec.k in 1:nrow(dat)){
      
      print(spec.k)
      #load presence absence
      pa= read.csv(paste("out/presabs/PresAbs_",dat$spec[spec.k],".csv",sep=""))
      
      #----
      #plot localities and temperature
      #crop to observed range 
      ext = extent(rbind(range(pa$lon), range(pa$lat))) # define the extent
      # extent
      ext[1]= ext[1]-10; ext[2]= ext[2]+10; ext[3]=ext[3]-10; ext[4]=ext[4]+10
      #crop
      tmax0=  crop(tmax_0, ext)
      tmax50=  crop(tmax_50, ext)
      tmax100=  crop(tmax_100, ext)
      
      #------
      #model thermoregulation
      #Set up prior
      CTmin1= dat$tmin[spec.k]
      CTmax1= dat$tmax[spec.k]
      #approximate Topt, but fix based on data
      Topt= CTmin1+ (CTmax1-CTmin1)*0.7
      
      #-----
      # sun to shade
      # thermoregulation scenario
      
      #max
      tmax0.dif= abs(tmax0 - Topt) 
      tmax50.dif= abs(tmax50 - Topt) 
      tmax100.dif= abs(tmax100 - Topt) 
      tmax.dif= stack(tmax0.dif, tmax50.dif, tmax100.dif)
      tr.ind= which.min(tmax.dif)
      
      tr<- tmax0
      tr[]<-NA
      tr[tr.ind==1]<- tmax0[tr.ind==1]
      tr[tr.ind==2]<- tmax50[tr.ind==2]
      tr[tr.ind==3]<- tmax100[tr.ind==3]
      trmax=tr
      names(trmax)<-"trmax"
      #----
      
      # Define prior specifications
      prior_specs <- list(
        # Temperature: Gaussian prior
        list(
          type = "gaussian", 
          optimal_value = 22, 
          sigma = 6, 
          lower_bound = 10, 
          upper_bound = 35
        ),
        # Precipitation: Polynomial prior
        list(
          type = "polynomial", 
          coefficients = c(0, 4, -4), 
          lower_bound = 500, 
          upper_bound = 1500, 
          normalize = TRUE
        ),
        # Elevation: Sigmoid prior
        list(
          type = "sigmoid", 
          midpoint = 1500, 
          steepness = 0.002, 
          lower_bound = 0, 
          upper_bound = 3000, 
          ascending = FALSE
        )
      )
      
      # Fit model with priors
      cat("Fitting GP-SDM with priors...\n")
      gp_model_with_priors <- fit_gp_sdm(
        data = example_data,
        response_col = "presence",
        env_vars = c("temperature", "precipitation", "elevation"),
        prior_specifications = prior_specs,
        combination_method = "multiplicative",
        kernel_type = "Matern52",
        seed = 123
      )
      
      # Print model summary
      print(gp_model_with_priors)
      
      # Make predictions using S3 method
      cat("\nMaking predictions...\n")
      predictions <- predict.GP_SDM(gp_model_with_priors, 
                                    newdata = example_data[1:10, ],
                                    include_prior = TRUE,
                                    return_se = TRUE,
                                    type = "response")
      
      print(head(predictions))
      
      
      
      # Split into training and testing sets
      train_idx <- sample(1:nrow(pa), 0.7 * nrow(pa))
      train_data <- pa[train_idx, ]
      test_data <- pa[-train_idx, ]
      
      #------------------------------
      # Define physiological prior
      my_prior <- function(data, temp_col, ctmin=dat$tmin[spec.k], ctmax=dat$tmax[spec.k]) {
        temperature_prior(
          data, 
          temp_col="trmax",
          ctmin = ctmin, 
          ctmax = ctmax,
          type= "tpc"
        )
      }
      
      #"beta", "gaussian", "poly", "threshold", "sigmoid", "tpc"
      
      my_prior2 <- function(data, temp_col, ctmin=dat$tmin[spec.k], ctmax=dat$tmax[spec.k]) {
        temperature_prior(
          data, 
          temp_col="trmax",
          ctmin = ctmin, 
          ctmax = ctmax,
          type= "sigmoid"
        )
      }
      
      #-----------------
      #make predictions
      pred= temperature_prior(pa, temp_col="trmax", ctmin=dat$tmin[spec.k], ctmax=dat$tmax[spec.k], type="tpc")
      #"beta", "gaussian", "poly", "threshold", "sigmoid", "tpc"
      
      #convert back from logit scale
      pa$pred= 1 / (1 + exp(-pred))
      
      #plot prior
      plot.prior<- ggplot() +
        geom_point(aes(trmax, pres), data=pa)+
        #add prior
        geom_line(aes(trmax, pred), data=pa)+
        xlim(0, 42)
      
      #--------------------------------
      # Optimize lengthscales
      
      # opt_result <- optimize_lengthscales(
      #   pres ~ trmax,
      #   data = train_data,
      #   mean_function = my_prior
      # )
      
      #roughly estimate length scales
      
      # Fit models with optimized lengthscales
      # Fit GP model without physiological prior
      gp_flat <- gp_sdm(
        pres ~ trmax,
        data = train_data,
        mean_function = my_prior2,
        lengthscales = 5 #opt_result$lengthscales #5
      )
      
      # Fit GP model with physiological prior
      gp_prior <- gp_sdm(
        pres ~ trmax,
        data = train_data,
        mean_function = my_prior,
        lengthscales = 5 #opt_result$lengthscales
      )
      
      # Make predictions
      pred_flat <- predict(gp_flat, test_data)
      pred_prior <- predict(gp_prior, test_data)
      
      # Compare models
      eval_flat <- evaluate_model(pred_flat, test_data$pres)
      eval_prior <- evaluate_model(pred_prior, test_data$pres)
      
      results <- data.frame(
        Model = c("GP-Flat", "GP-Prior"),
        AUC = c(eval_flat$AUC, eval_prior$AUC),
        RMSE = c(eval_flat$RMSE, eval_prior$RMSE)
      )
      
      models[spec.k,1:3]<- c(results[1,])
      models[spec.k,4:6]<- c(results[2,])
      
      # # Cross validate
      # mean_functions <- list(
      #   "Flat" = function(data) rep(0, nrow(data)),
      #   # "Temperature_Only" = function(data) {
      #   #   temp_prob <- temp_response(data$trmax, optimal_temp = 22, tolerance = 6)
      #   #   return(log(temp_prob / (1 - temp_prob)))
      #   # },
      #   "Full_Physiological" = my_prior
      # )
      # 
      # # Run cross-validation
      # cv_results <- cross_validate_gp(
      #   pres ~ trmax,
      #   data = train_data,
      #   mean_functions = mean_functions,
      #   k = 5,
      #   lengthscales = opt_result$lengthscales
      # )
      # 
      # # Print summary
      # print(cv_results$summary)
      
      #---------------------------
      #point based approach
      #extract values
      pts= rasterToPoints(trmax)
      colnames(pts)=c("lon","lat","trmax")
      pts= as.data.frame(pts)
      pts$pres<- 0
      
      #predictions
      pred_flat <- predict(gp_flat, newdata=pts, type="response")
      pred_prior <- predict(gp_prior, newdata=pts, type="response")
      
      #combine predictions
      pts= cbind(pts, pred_flat[,1], pred_prior[,1])
      colnames(pts)[(ncol(pts)-1):ncol(pts)]=c("occ_pp","occ_np")
      
      #plot
      #to long format
      pts.l <- melt(pts[,which(names(pts)!="pres")], id=c("lon","lat","trmax"))
      #presence points
      pres= subset(pa, pa$pres=="1")
      pres$variable="occ_pp"
      #replicate pres for plotting
      pres2<- pres
      pres2$variable="occ_np"
      pres.all= rbind(pres, pres2)
      
      #trmax
      tr.plot=ggplot(pts, aes(lon, lat, fill= trmax)) + 
        geom_tile()+scale_fill_viridis(na.value = 'grey')
      
      #predictions
      occ.plot= ggplot(pts.l, aes(lon, lat)) + 
        geom_tile(aes(fill= value))+scale_fill_viridis(na.value = 'grey') +facet_wrap(~variable, nrow=1)+
        ggtitle(dat$species[spec.k])
      #add localities
      occ.plot= occ.plot +geom_point(pres.all, mapping=aes(lon, lat, color="red"))
      
      #combine and write out
      design <- "ABBB"
      
      #save figure 
      pdfname <- paste("./figures/",dat$species[spec.k], ".pdf", sep="")
      pdf(pdfname,height = 6, width = 10)
      print(plot.prior + occ.plot + plot_layout(design = design))
      dev.off()
      
    } #end loop species
    
    #-----------
    #examine results 
    plot(models$AUC1, models$AUC2)
    abline(a=0, b=1)
    
    plot(models$RMSE1, models$RMSE2)
    abline(a=0, b=1)
    
    #write out summary
    write.csv(models, "./out/model_stats.csv")
    
    
    