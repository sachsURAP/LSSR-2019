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
# Mod history:  16 Apr 2018
# Details:      See hgData.R for further licensing, attribution, references, 
#               and abbreviation information.

source("synergyTheory.R") # load in data and models

library(mvtnorm) #  Sampling

# helper function to generate samples
.generate_samples <- function(N = 200, r, L, d, model, HINmodel = HZE_nte_model, 
                             HITmodel = HZE_te_model, LOWmodel = low_LET_model, 
                             calib = TRUE) {
  # Function to generate Monte Carlo samples for calculating CI
  # @params:   N              - numbers of sample
  #            model          - select HIN or HIT model
  #                             0 - HIT model
  #                             1 - HIN model
  #            HINmodel       - the input HIN model
  #            HITmodel       - the input HIN model
  #            LOWmodel       - the input LOW model
  monteCarloSamplesLow <- rmvnorm(n = N, mean = coef(LOWmodel), sigma = vcov(LOWmodel))
  curve_list <- list(0)
  if (model) {
    if (calib) {
      monteCarloSamplesHin <- rmvnorm(n = N, mean = coef(HINmodel), sigma = vcov(HINmodel))
    } else {
      monteCarloSamplesHin <- mapply(rnorm, rep(N, 3), coef(HINmodel), 
                                     summary(HINmodel)$coefficients[, "Std. Error"])
    }
  } else {
    if (calib) {
      monteCarloSamplesHit <- rmvnorm(n = N, mean = coef(HITmodel), sigma = vcov(HITmodel))
    } else {
      monteCarloSamplesHit <- mapply(rnorm, rep(N, 2), coef(HITmodel), 
                                     summary(HITmodel)$coefficients[, "Std. Error"])
    }
  }
  for (i in 1:N) {
    if (model) { 
      curve_list[[i]] <- calculate_complex_id(r = r, LET = L, d = d, 
                                              coef = list(NTE = monteCarloSamplesHin[i, ],
                                              lowLET = monteCarloSamplesLow[i]),
                                              model = "NTE", lowLET = TRUE)
    } else {
      curve_list[[i]] <- calculate_complex_id(r = r, LET = L, d = d, 
                                              coef = list(TE = monteCarloSamplesHit[i, ],
                                              lowLET = monteCarloSamplesLow[i]),
                                              model = "TE", lowLET = TRUE)
    }
    cat(paste("Currently at Monte Carlo step:", toString(i), "of", toString(N)), sprintf('\r'))
  }
  return(curve_list)
}

.generate_ci <- function(N = 200, intervalLength = 0.95, doseIndex, r, 
                        L, HINmodel = HZE_nte_model, HITmodel = HZE_te_model, 
                        LOWmodel = low_LET_model, sampleCurves, model = mod) {
  # Function to generate CI for the input dose.
  # @params:   N              - numbers of sample
  #            intervalLength - size of confidence interval
  #            d              - input dose
  #            r              - proportion of ion
  #            L              - LTE
  #            HINmodel       - the input HIN model
  #            HITmodel       - the input HIN model
  #            LOWmodel       - the input LOW model
  #            method         - select Naive or Monte Carlo Approach
  #                             0 - Naive
  #                             1 - Monte Carlo
  
  # For each sample curve, evalute them at input dose, and sort.
  valueArr <- vector(length = 0)
  for (i in 1:N) {
    valueArr <- c(valueArr, sampleCurves[[i]][, 2][doseIndex])
  }
  valueArr <- sort(valueArr)
  # Returning resulting CI
  return(c(valueArr[ceiling((1 - intervalLength) / 2 * N)], 
            valueArr[(intervalLength + (1 - intervalLength) / 2) * N]))
}

#############
# CI Helper #
#############
ci_helper <- function(sample_num = 200, d, r, L,
                      model = 1, seed = 100, calib = TRUE) {
  # Set the pseudorandom seed
  set.seed(seed)
  # Generate N randomly generated samples of parameters of HZE model.
  curve_list <- .generate_samples(N = sample_num, r = r, L = L, d = d, 
                                 model = model, calib = calib)
  numDosePoints <- length(d)
  monte_carlo_ci <- matrix(nrow = 2, ncol = numDosePoints)
  
  # Calculate CI for each dose point
  for (i in 1 : numDosePoints) {
    monte_carlo_ci[, i] <- .generate_ci(N = sample_num, doseIndex = i, r = r, L = L, 
                                       model = model, sampleCurves = curve_list)
    cat(paste("Iterating on dose points. Currently at step:", toString(i), "of", 
              toString(numDosePoints)), sprintf('\r'))
  }
  return(list(monte_carlo = monte_carlo_ci))
}
