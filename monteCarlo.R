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
# Mod history:  26 Apr 2018
# Details:      See hgData.R for further licensing, attribution, references, 
#               and abbreviation information.

source("synergyTheory.R") # load in data and models

library(mvtnorm) # Sampling

#======================= MONTE CARLO SIMULATION FUNCTION ======================#

#' @description Runs the Monte Carlo method on a MIXDER.
#' 
#' @param sample_num Numeric integer of the number of samples to be drawn.
#' @param d Numeric vector of all total dose values to be evaluated. 
#' @param r Numeric vector of the dose proportions applied to component IDERs.
#' @param l Numeric vector of all LET values, must be length n.
#' @param model String value corresponding to the model to be used. 
#' @param seed Numeric value for pseudorandom generators.
#' @param vcov Boolean for assessing inter-parameter correlation.
#'             
#' @details Corresponding elements of ratios, LET should be associated with the
#'          same IDER.
#'          
#' @return Numeric vector representing lower and upper bounds for a Monte Carlo
#'         confidence interval for the given MIXDER at the given dosages.
#'         
#' @examples

simulate_monte_carlo <- function(sample_num = 200, d, r, L,
                                 model = "NTE", seed = 100, vcov = TRUE) {
  # Set the pseudorandom seed
  set.seed(seed)
  # Generate N randomly generated samples of parameters of HZE model.
  curve_list <- .generate_samples(sample_num, d, L, r, model, vcov)
  numDosePoints <- length(d)
  monte_carlo_ci <- matrix(nrow = 2, ncol = numDosePoints)
  
  # Calculate CI for each dose point
  for (i in 1 : numDosePoints) { #EGH: Possible vectorization opportunity
    monte_carlo_ci[, i] <- .generate_ci(N = sample_num, dose_index = i, 
                                        sample_curves = curve_list)
  }
  return(list(monte_carlo = monte_carlo_ci))
}

#========================= MONTE CARLO HIDDEN FUNCTIONS =======================#

#================== SAMPLING ===================#

#' @description Generates MIXDER samples using parameters drawn from a Gaussian 
#'              distribution.
#'              
#' @param N Numeric integer of the number of samples to be drawn.
#' @param d Numeric vector of all total dose values to be evaluated. 
#' @param r Numeric vector of the dose proportions applied to component IDERs.
#' @param l Numeric vector of all LET values, must be length n.
#' @param model String value corresponding to the model to be used.
#' @param vcov Boolean for assessing inter-parameter correlation. 
#'             TRUE iff correlations are being taken into account. #RKS to EGH. Is this correct? #EGH Apr 16: Correct. I will update function documentation in the next few commmits.
#' @param HINmodel The input HIN model
#' @param HITmodel The input HIT model # RKS to EGH. Corrected misprint. #EGH Apr 16: Thanks.
#' @param LOWmodel The input LOW model
#' 
#' @details Corresponding elements of ratios, LET should be associated with the
#'          same IDER.
#'          
#' @return Numeric vector of sample MIXDERs evaluated at the given dosages. 
#' 
#' @examples

.generate_samples <- function(N = 200, d, L, r, model, vcov, 
                              HINmodel = HZE_nte_model, HITmodel = HZE_te_model, 
                              LOWmodel = low_LET_model) {
  low_LET_samples <- rmvnorm(N, coef(LOWmodel), vcov(LOWmodel)) # RKS to RKS and EGH. Is vcov a 1x1 matrix here? #EGH Apr 16: Correct.
  # RKS to EGH. I worried that your Monte Carlo was fast whereas that of Dae and the CA pod is much slower.
  # So I tried to understand the next chunk that defines curve_list, but got lost. #EGH Apr 16: That's fair. Funnily enough, I was worried that this Monte Carlo seemed too slow, and I was looking into ways to speed it up.
  # How should we handle this kind of situation to not lose information but not clutter up the script? #EGH Apr 16: Are you familiar with the browser() function? It allows you to freeze a function mid-evaluation and examine its environment and make other diagnostic function calls. 
  if (model == "NTE") {
    ion_model <- HINmodel
    num_coef <- 3
  } else if (model == "TE") {
    ion_model <- HITmodel
    num_coef <- 2
  }
  if (vcov) {
    samples <- rmvnorm(n = N, mean = coef(ion_model), sigma = vcov(ion_model))
  } else {
    samples <- mapply(rnorm, rep(N, num_coef), coef(ion_model), 
                      summary(ion_model)$coefficients[, "Std. Error"])
  }
  curve_list <- list(0)
  
  for (i in 1:N) { #EGH: Possible vectorization opportunity. Next block constitutes the bulk of the Monte Carlo runtime.
    if (model == "NTE") {
      curve_list[[i]] <- calculate_id(d, L, r, coef = list(NTE = samples[i, ],
                                               lowLET = low_LET_samples[i]),
                                               model = "NTE", lowLET = TRUE)
    } else if (model == "TE") {
      curve_list[[i]] <- calculate_id(d, L, r, coef = list(TE = samples[i, ],
                                               lowLET = low_LET_samples[i]),
                                               model = "TE", lowLET = TRUE)
    }
    cat(paste("  Currently at Monte Carlo step:", toString(i), 
              "of", toString(N)), sprintf('\r'))
  }
  return(curve_list)
}

#============= INTERVAL CONSTRUCTION ===========#

#' @description Generates confidence interval bounds for given MIXDER samples 
#'              at a given dosage.
#'              
#' @param N Numeric integer of the number of samples to be drawn.
#' @param interval_length Numeric double of the confidence interval width.
#' @param dose_index Numeric integer of dosage.
#' @param sample_curves Numeric list of sampled MIXDER values.
#' 
#' @details 
#' 
#' @return Numeric length-two vector of an upper and lower bound for the 
#'         confidence interval at the given dosage.
#'         
#' @examples

.generate_ci <- function(N, interval_length = 0.95, dose_index, sample_curves) {
  # For each sample curve, evalute them at input dose, and sort.
  sample_values <- sort(sapply(sample_curves, function(x) x[, 2][dose_index]))
  # Returning resulting CI
  return(c(sample_values[ceiling((1 - interval_length) / 2 * N)], 
           sample_values[(interval_length + (1 - interval_length) / 2) * N]))
}
