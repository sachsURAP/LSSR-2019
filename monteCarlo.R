#   Filename: monteCarlo.R 
#   Purpose: Concerns radiogenic mouse HG tumorigenesis. Contains the code to 
#            run Monte Carlo and generate confidence intervals for our models.

source("synergyTheory.R") # load in data and models

library(mvtnorm)
#==============================================#
#==========Confidence Interval Part============#
#==============================================#

# Parameter initialization
# r <- c(0.05, 0.05, 0.05, 0.05, 0.8); L <- c(25, 70, 100, 195)
sample_num = 200
mod = 1              # 1 if HIN, 0 if HIT

# d <- c(seq(0, .01, by = 0.001),
       # seq(.02, .1, by=.01),
       # seq(.2, 1, by=.1),
       # seq(2, 100, by=1))

# Set the pseudorandom seed
set.seed(100)

# helper function to generate samples
generate_samples <- function(N = sample_num, model = mod, HINmodel = HZE_nte_model, HITmodel = HZE_te_model, LOWmodel = low_LET_model, r, L, d) {
  # Function to generate Monte Carlo samples for calculating CI
  # @params:   N              - numbers of sample
  #            model          - select HIN or HIT model
  #                             0 - HIT model
  #                             1 - HIN model
  #            HINmodel       - the input HIN model
  #            HITmodel       - the input HIN model
  #            LOWmodel       - the input LOW model
  monteCarloSamplesLow = rmvnorm(n = N, mean = coef(LOWmodel), sigma = vcov(LOWmodel))
  curve_list = list(0)
  if (model) {
    monteCarloSamplesHin = rmvnorm(n = N, mean = coef(HINmodel), sigma = vcov(HINmodel))
  } else {
    monteCarloSamplesHit = rmvnorm(n = N, mean = coef(HITmodel), sigma = vcov(HITmodel))
  }
  for (i in 1:N) {
    if (model) { 
      curve_list[[i]] <- calculate_complex_id(r, L,d, coef = list(NTE = monteCarloSamplesHin[i, ],
                                                                 lowLET = monteCarloSamplesLow[i]),
                                            model = "NTE", lowLET = TRUE)
    } else {
      curve_list[[i]] <- calculate_complex_id(r, L, d, coef = list(TE = monteCarloSamplesHit[i, ],
                                                                  lowLET = monteCarloSamplesLow[i]),
                                            model = "TE", lowLET = TRUE)
    }
    cat(paste("Currently at Monte Carlo step:", toString(i), "of", toString(N)), sprintf('\r'))
    cat(sprintf('\r'))
  }
  return (curve_list)
}

# Generate N randomly generated samples of parameters of HZE model.
# curve_list <- generate_samples(N = sample_num, r = r, L = L, d = d)

generate_ci <- function(N = sample_num, intervalLength = 0.95, d, doseIndex, r, L, HINmodel = HZE_nte_model, HITmodel = HZE_te_model, LOWmodel = low_LET_model, method = 0, sampleCurves, model = mod) {
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
  if (method) {
    # For each sample curve, evalute them at input dose, and sort.
    valueArr = vector(length = 0)
    for (i in 1:N) {
      valueArr = c(valueArr, sampleCurves[[i]][, 2][doseIndex])
    }
    valueArr = sort(valueArr)
    # Returning resulting CI
    return (c(valueArr[(1 - intervalLength) / 2 * N], valueArr[(intervalLength + (1 - intervalLength) / 2) * N]))
  } else {
    #========= Naive =========#
    stdErrArrLow = summary(LOWmodel)$coefficients[, "Std. Error"]
    meanArrLow = summary(LOWmodel)$coefficients[, "Estimate"]
    if (model) {
      stdErrArrHin = summary(HINmodel)$coefficients[, "Std. Error"]
      meanArrHin = summary(HINmodel)$coefficients[, "Estimate"]
      upper = calculate_complex_id(r, L, d = c(0, d),
                                   coef = list(NTE = c(
                                   meanArrHin["aa1"] + 2*stdErrArrHin["aa1"],
                                   meanArrHin["aa2"] + 2*stdErrArrHin["aa2"],
                                   meanArrHin["kk1"] + 2*stdErrArrHin["kk1"]),
                                   lowLET = meanArrLow + 2*stdErrArrLow),
                                   lowLET = TRUE)[, 2][2]
      lower = calculate_complex_id(r, L, d = c(0, d),
                                   coef = list(NTE = c(
                                   meanArrHin["aa1"] - 2*stdErrArrHin["aa1"],
                                   meanArrHin["aa2"] - 2*stdErrArrHin["aa2"],
                                   meanArrHin["kk1"] - 2*stdErrArrHin["kk1"]),
                                   lowLET = meanArrLow - 2*stdErrArrLow),
                                   lowLET = TRUE)[, 2][2]
    } else {
      stdErrArrHit = summary(HITmodel)$coefficients[, "Std. Error"]
      meanArrHit = summary(HITmodel)$coefficients[, "Estimate"]
      upper = calculate_complex_id(r, L, d = c(0, d),
                                   coef = list(TE = c(meanArrHit["aate1"] + 2*stdErrArrHit["aate1"],
                                                      meanArrHit["aate2"] + 2*stdErrArrHit["aate2"]),
                                               lowLET = meanArrLow + 2*stdErrArrLow),
                                   model = "TE", lowLET = TRUE)[, 2][2]
      lower = calculate_complex_id(r, L, d = c(0, d),
                                   coef = list(TE = c(meanArrHit["aate1"] - 2*stdErrArrHit["aate1"],
                                                      meanArrHit["aate2"] - 2*stdErrArrHit["aate2"]),
                                               lowLET = meanArrLow - 2*stdErrArrLow),
                                   model = "TE", lowLET = TRUE)[, 2][2]
    }
    return (c(lower, upper))
  }
}

# Parameter initialization
# if (mod) {
#   mixderCurve = calculate_complex_id(r, L, d = d, lowLET = TRUE)
# } else {
#   mixderCurve = calculate_complex_id(r, L, d = d, lowLET = TRUE, model = "TE")
# }
# fourIonMIXDER = data.frame(d = mixderCurve[, 1], CA = mixderCurve[, 2])
# numDosePoints = length(fourIonMIXDER$d)
# naive_ci = matrix(nrow = 2, ncol = numDosePoints)
# monte_carlo_ci = matrix(nrow = 2, ncol = numDosePoints)

# Calculate CI for each dose point
# for (i in 1 : numDosePoints) {
#   naive_ci[, i] = generate_ci(d = fourIonMIXDER$d[i], r = r,  L = L, sampleCurves = curve_list)
#   monte_carlo_ci[, i] = generate_ci(doseIndex = i, r = r,  L = L, method = 1, sampleCurves = curve_list)
#   cat(paste("Iterating on dose points. Currently at step:", toString(i), "of", toString(numDosePoints)), sprintf('\r'))
#   
# }

#############
# CI Helper #
#############
ci_helper = function(sample_num, intervalLength = 0.95, d, r, L, HINmodel = HZE_nte_model, HITmodel = HZE_te_model, LOWmodel = low_LET_model, model = mod, seed = 100) {
  # Set the pseudorandom seed
  set.seed(seed)
  
  # Generate N randomly generated samples of parameters of HZE model.
  curve_list <- generate_samples(N = sample_num, r = r, L = L, d = d)
  
  # Parameter initialization
  if (mod) {
    mixderCurve <- calculate_complex_id(r, L, d = d, lowLET = TRUE)
  } else {
    mixderCurve <- calculate_complex_id(r, L, d = d, lowLET = TRUE, model = "TE")
  }
  fourIonMIXDER = data.frame(d = mixderCurve[, 1], CA = mixderCurve[, 2])
  numDosePoints = length(fourIonMIXDER$d)
  naive_ci = matrix(nrow = 2, ncol = numDosePoints)
  monte_carlo_ci = matrix(nrow = 2, ncol = numDosePoints)
  
  # Calculate CI for each dose point
  for (i in 1 : numDosePoints) {
    naive_ci[, i] = generate_ci(d = fourIonMIXDER$d[i], r = r,  L = L, sampleCurves = curve_list)
    monte_carlo_ci[, i] = generate_ci(doseIndex = i, r = r,  L = L, method = 1, sampleCurves = curve_list)
    cat(paste("Iterating on dose points. Currently at step:", toString(i), "of", toString(numDosePoints)), sprintf('\r'))
  }
  return (list(naive = naive_ci, monte_carlo = monte_carlo_ci))
}

#==============================================#
#==========Confidence Interval End=============#
#==============================================#

# Plotting checks
# plot(y = .01, x = 0.1, xlim = c(0, 128), ylim = c(0, 0.3))
# lines(1:128, naive_ci[1, ], col = "red")
# lines(1:128, naive_ci[2, ], col = "red")
# lines(1:128, monte_carlo_ci[1, ], col = "blue")
# lines(1:128, monte_carlo_ci[2, ], col = "blue")

