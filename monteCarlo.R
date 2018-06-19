# Copyright:    (C) 2017-2018 Sachs Undergraduate Research Apprentice Program
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3.
# Filename:     monteCarlo.R
# Purpose:      Concerns radiogenic mouse Harderian gland tumorigenesis. 
#               Contains the code to run Monte Carlo sampling and generate 
#               confidence intervals for dose-effect relationship models. It is 
#               part of the source code for the NASAmouseHG project.
# Contact:      Rainer K. Sachs 
# Website:      https://github.com/sachsURAP/NASAmouseHG
# Mod history:  18 Jun 2018
# Details:      See dataAndInfo.R for further licensing, attribution, 
#               references, and abbreviation information.

source("synergyTheory.R") # Load in data and models.

library(mvtnorm) # Random sampling.

#======================= MONTE CARLO SIMULATION FUNCTION ======================#

#' @description Runs the Monte Carlo method on a baseline no-synergy/antagonism 
#'              mixture DER.
#' 
#' @param n Numeric integer of the number of samples to be drawn.
#' @param dose Numeric vector of all total dose values to be evaluated. 
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of dose proportions applied on component DERs.
#' @param model String value corresponding to the model to be used. 
#' @param vcov Boolean for assessing inter-parameter correlation.
#' @param interval_length Numeric double of the confidence interval width.
#' @param seed Numeric value for pseudorandom generators.
           
#' @details Corresponding elements of ratios, LET should be associated with the
#'          same DER.
#'          
#' @return Named list representing lower and upper bounds for a Monte Carlo
#'         confidence interval for a mixture DER over an interval of doses.
#'         
#' @examples
#' ratios <- c(1/2, 1/2)
#' LET_vals <- c(195, 70)
#' simulate_monte_carlo(0:100, LET_vals, ratios)

simulate_monte_carlo <- function(n = 200, dose, LET, ratios, model = "NTE", 
                                 vcov = TRUE, interval_length = 0.95,
                                 seed = 100) {
  # Set the pseudorandom seed
  set.seed(seed)
  # Generate N randomly generated samples of parameters of HZE model.
  curve_list <- .generate_samples(n, dose, LET, ratios, model, vcov)
  monte_carlo_ci <- matrix(nrow = 2, ncol = length(dose))
  
  # Calculate CI for each dose point
  for (i in 1 : length(dose)) { # EGH: Possible vectorization opportunity
    monte_carlo_ci[, i] <- .generate_ci(n, i, curve_list, interval_length)
  }
  return(list(monte_carlo = monte_carlo_ci))
}


#========================= MONTE CARLO HIDDEN FUNCTIONS =======================#

#================== SAMPLING ===================#

#' @description Generates mixture no-synergy/antagonism DER samples with
#'              parameters drawn from a Gaussian distribution.
#'              
#' @param n Numeric integer of the number of samples to be drawn.
#' @param dose Numeric vector of all total dose values to be evaluated. 
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of dose proportions applied on component DERs.
#' @param model String value corresponding to the model to be used.
#' @param vcov Boolean for assessing inter-parameter correlation. 
#' @param nte_model The input HZE nontargeted and targeted effects model.
#' @param te_model The input HZE targeted effects only model.
#' @param low_model The input low LET model.
#'          
#' @return Numeric vector of sample mixture baseline DERs evaluated at the given doses. 

.generate_samples <- function(n, dose, LET, ratios, model, vcov, 
                              nte_model = HZE_nte_model, 
                              te_model = HZE_te_model, 
                              low_model = low_LET_model) {
  low_LET_samples <- rmvnorm(n, coef(low_model), vcov(low_model)) 
  if (model == "NTE") {
    ion_model <- nte_model
    num_coef <- 3
  } else if (model == "TE") {
    ion_model <- te_model
    num_coef <- 2
  }
  if (vcov) {
    samples <- rmvnorm(n, mean = coef(ion_model), sigma = vcov(ion_model))
  } else {
    samples <- mapply(rnorm, rep(n, num_coef), coef(ion_model), 
                      summary(ion_model)$coefficients[, "Std. Error"])
  }
  curve_list <- list(0)
  for (i in 1:n) {
    if (model == "NTE") {
      curve_list[[i]] <- calculate_id(dose, LET, ratios, 
                                      coef = list(NTE = samples[i, ],
                                      lowLET = low_LET_samples[i]),
                                      model = "NTE")
    } else if (model == "TE") {
      curve_list[[i]] <- calculate_id(dose, LET, ratios,
                                      coef = list(TE = samples[i, ],
                                      lowLET = low_LET_samples[i]),
                                      model = "TE")
    }
    cat(paste("  Currently at Monte Carlo step:", toString(i), "of", 
              toString(n)), sprintf('\r'))
  }
  return(curve_list)
}

#============= INTERVAL CONSTRUCTION ===========#

#' @description Generates confidence intervals for mixture baeline DER samples 
#'              at a dose.
#'              
#' @param n Numeric integer of the number of samples to be drawn.
#' @param dose_index Numeric integer of doses
#' @param sample_curves Numeric list of sampled mixture DER values.
#' @param interval_length Numeric double of the confidence interval width.
#' 
#' @return Numeric length-two vector of an upper and lower bound for the 
#'         confidence interval of a dose.

.generate_ci <- function(n, dose_index, sample_curves, interval_length) {
  # For each sample curve, evalute them at input dose, and sort.
  sample_values <- sort(sapply(sample_curves, function(x) x[, 2][dose_index]))
  # Returning resulting CI
  return(c(sample_values[ceiling((1 - interval_length) / 2 * n)], 
           sample_values[(interval_length + (1 - interval_length) / 2) * n]))
}
