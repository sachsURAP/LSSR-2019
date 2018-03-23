#   Filename: HGsynergyMain.R 
#   Purpose: Concerns radiogenic mouse HG tumorigenesis. Contains the code to 
#            run Monte Carlo and generate confidence intervals for our models.

source("HGsynergyMain.R") # load in the data
source("synergytheory.R") # load in models

#==============================================#
#==========Confidence Interval Part============#
#==============================================#

# Set the pseudorandom seed

set.seed(1)
#
Generate_CI <- function(N = 500, intervalLength = 0.95, d, r, L, HZEmodel = hi_nte_model, method = 0) {
  # Function to generate CI for the input dose.
  # @params:   N              - numbers of sample
  #            intervalLength - size of confidence interval
  #            d              - input dose
  #            r              - proportion of ion
  #            L              - LTE
  #            HZEmodel       - the input HZE model
  #            method         - select Naive or Monte Carlo Approach
  #                             0 - Naive
  #                             1 - Monte Carlo
  if (method) {
    #========= Monte Carlo =========#
    valueArr = vector(length = 0)
    # Generate N randomly generated samples of parameters of HZE model.
    monteCarloSamples = rmvnorm(n = N, mean = coef(HZEmodel), sigma = vcov(HZEmodel))
    
    # For each sample curve, evalute them at input dose, and sort.
    for (i in 1:500) {
      # print(calculate_complex_id(r = r, L = L, d = c(0, d), aa1 = monteCarloSamples[, 1][i], aa2 = monteCarloSamples[, 2][i], kk1 = monteCarloSamples[, 3][i]))[1]
      # browser()
      #   valueArr = c(valueArr, calculate_complex_id(r = r, L = L, d = c(0, d), aa1 = monteCarloSamples[, 1][i], aa2 = monteCarloSamples[, 2][i], kk1 = monteCarloSamples[, 3][i])[, 2][2])
      valueArr = c(valueArr, calculate_complex_id(r = r, L = L, d = c(0, d), coef = list(NTE = c(aa1 = monteCarloSamples[, 1][i],
                                                                                                 aa2 = monteCarloSamples[, 2][i],
                                                                                                 kk1 = monteCarloSamples[, 3][i])))[, 2][2])
    }
    valueArr = sort(valueArr)
    
    # Returning resulting CI
    return (c(valueArr[(1-intervalLength)/2*500], valueArr[(intervalLength + (1-intervalLength)/2)*500]))
  } else {
    #========= Naive =========#
    stdErrArr = summary(HZEmodel)$coefficients[, "Std. Error"]
    meanArr = summary(HZEmodel)$coefficients[, "Estimate"]
    # upper = calculate_complex_id(r = r, L = L, d = c(0, d), aa1 = meanArr["aa1"] + 2*stdErrArr["aa1"], aa2 = meanArr["aa2"] + 2*stdErrArr["aa2"], kk1 = meanArr["kk1"] + 2*stdErrArr["kk1"])[, 2][2]
    # lower = calculate_complex_id(r = r, L = L, d = c(0, d), aa1 = meanArr["aa1"] - 2*stdErrArr["aa1"], aa2 = meanArr["aa2"] - 2*stdErrArr["aa2"], kk1 = meanArr["kk1"] - 2*stdErrArr["kk1"])[, 2][2]
    upper = calculate_complex_id(r = r, L = L, d = c(0, d),
                                 coef = list(NTE = c(aa1 = meanArr["aa1"] + 2*stdErrArr["aa1"],
                                                     aa2 = meanArr["aa2"] + 2*stdErrArr["aa2"],
                                                     kk1 = meanArr["kk1"] + 2*stdErrArr["kk1"])))[, 2][2]
    lower = calculate_complex_id(r = r, L = L, d = c(0, d),
                                 coef = list(NTE = c(aa1 = meanArr["aa1"] - 2*stdErrArr["aa1"],
                                                     aa2 = meanArr["aa2"] - 2*stdErrArr["aa2"],
                                                     kk1 = meanArr["kk1"] - 2*stdErrArr["kk1"])))[, 2][2]
    return (c(lower, upper))
  }
}

# Parameter initialization
r <- c(1/3, 1/3, 1/3)
L <- c(25, 70, 250)
mixderCurve = calculate_complex_id(r, L, d = dose_vector)
threeIonMIXDER = data.frame(d = mixderCurve[, 1], CA = mixderCurve[, 2])
numDosePoints = length(threeIonMIXDER$d)
naiveCI = matrix(nrow = 2, ncol = numDosePoints)
monteCarloCI = matrix(nrow = 2, ncol = numDosePoints)

# Calculate CI for each dose point
for (i in 1 : numDosePoints) {
  naiveCI[, i] = Generate_CI(d = threeIonMIXDER$d[i], r = r,  L = L)
  monteCarloCI[, i] = Generate_CI(d = threeIonMIXDER$d[i], r = r,  L = L, method = 1)
  print(paste("Currently at step:", toString(i)))
}

# Plot
mixderGraphWithNaiveCI = ggplot(data = threeIonMIXDER, aes(x = d, y = CA)) + geom_line(aes(y = CA), col = "red", size = 1) + geom_ribbon(aes(ymin = naiveCI[1, ], ymax = naiveCI[2, ]), alpha = .2)
mixderGraphWithMonteCarloCI = ggplot(data = threeIonMIXDER, aes(x = d, y = CA)) + geom_line(aes(y = CA), col = "red", size = 1) + geom_ribbon(aes(ymin = monteCarloCI[1, ], ymax = monteCarloCI[2, ]), alpha = .4)
print(mixderGraphWithNaiveCI)
print(mixderGraphWithMonteCarloCI)

mixderGraphWithNaiveAndMonteCarloCI = ggplot(data = threeIonMIXDER, aes(x = d, y = CA)) + geom_line(aes(y = CA), col = "red", size = 1) + geom_ribbon(aes(ymin = monteCarloCI[1, ], ymax = monteCarloCI[2, ]), alpha = .4) + geom_line(aes(y = CA), col = "red", size = 1) + geom_ribbon(aes(ymin = naiveCI[1, ], ymax = naiveCI[2, ]), alpha = .2)
print(mixderGraphWithNaiveAndMonteCarloCI)

#========================================#
#===================End==================#
#========================================#

#==============================================#
#==========Confidence Interval Part============#
#==============================================#

# Parameter initialization
r <- c(0.05, 0.05, 0.05, 0.05, 0.8); L <- c(25, 70, 100, 195)
sampleNum = 200
mod = 1              # 1 if HIN, 0 if HIT

d <- c(seq(0, .00001, by = 0.000001), 
       seq(.00002, .0001, by=.00001),
       seq(.0002, .001, by=.0001),
       seq(.002, .01, by=.001),
       seq(.02, 1., by=.01))

# Set the pseudorandom seed
set.seed(100)

# helper function to generate samples
Generate_samples = function(N = sampleNum, model = mod, HINmodel = hinm, HITmodel = hitm, LOWmodel = LOW.m, r, L, d) {
  # Function to generate Monte Carlo samples for calculating CI
  # @params:   N              - numbers of sample
  #            model          - select HIN or HIT model
  #                             0 - HIT model
  #                             1 - HIN model
  #            HINmodel       - the input HIN model
  #            HITmodel       - the input HIN model
  #            LOWmodel       - the input LOW model
  monteCarloSamplesLow = rmvnorm(n = N, mean = coef(LOWmodel), sigma = vcov(LOWmodel))
  curveList = list(0)
  if (model) {
    monteCarloSamplesHin = rmvnorm(n = N, mean = coef(HINmodel), sigma = vcov(HINmodel))
  } else {
    monteCarloSamplesHit = rmvnorm(n = N, mean = coef(HITmodel), sigma = vcov(HITmodel))
  }
  for (i in 1:N) {
    if (model) {
      curveList[[i]] = calculateComplexId(r = r, L = L, d = d, aa1 = monteCarloSamplesHin[, 1][i], aa2 = monteCarloSamplesHin[, 2][i], kk1 = monteCarloSamplesHin[, 3][i], beta = monteCarloSamplesLow[, 1][i], lowLET = TRUE)
    } else {
      curveList[[i]] = calculateComplexId.te(r = r, L = L, d = d, aate1 = monteCarloSamplesHit[, 1][i], aate2 = monteCarloSamplesHit[, 2][i], beta = monteCarloSamplesLow[, 1][i], lowLET = TRUE)
    }
    print(paste("Currently at Monte Carlo step:", toString(i), "Total of", toString(N), "steps"))
  }
  return (curveList)
}

# Generate N randomly generated samples of parameters of HZE model.
curveList = Generate_samples(N = sampleNum, r = r, L = L, d = d)

Generate_CI = function(N = sampleNum, intervalLength = 0.95, d, doseIndex, r, L, HINmodel = hinm, HITmodel = hitm, LOWmodel = LOW.m, method = 0, sampleCurves, model = mod) {
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
    return (c(valueArr[(1-intervalLength)/2*N], valueArr[(intervalLength + (1-intervalLength)/2)*N]))
  } else {
    #========= Naive =========#
    stdErrArrLow = summary(LOWmodel)$coefficients[, "Std. Error"]
    meanArrLow = summary(LOWmodel)$coefficients[, "Estimate"]
    if (model) {
      stdErrArrHin = summary(HINmodel)$coefficients[, "Std. Error"]
      meanArrHin = summary(HINmodel)$coefficients[, "Estimate"]
      upper = calculateComplexId(r = r, L = L, d = c(0, d), aa1 = meanArrHin["aa1"] + 2*stdErrArrHin["aa1"], aa2 = meanArrHin["aa2"] + 2*stdErrArrHin["aa2"], kk1 = meanArrHin["kk1"] + 2*stdErrArrHin["kk1"], beta = meanArrLow + 2*stdErrArrLow, lowLET = TRUE)[, 2][2]
      lower = calculateComplexId(r = r, L = L, d = c(0, d), aa1 = meanArrHin["aa1"] - 2*stdErrArrHin["aa1"], aa2 = meanArrHin["aa2"] - 2*stdErrArrHin["aa2"], kk1 = meanArrHin["kk1"] - 2*stdErrArrHin["kk1"], beta = meanArrLow - 2*stdErrArrLow, lowLET = TRUE)[, 2][2]
    } else {
      stdErrArrHit = summary(HITmodel)$coefficients[, "Std. Error"]
      meanArrHit = summary(HITmodel)$coefficients[, "Estimate"]
      upper = calculateComplexId.te(r = r, L = L, d = c(0, d), aate1 = meanArrHit["aate1"] + 2*stdErrArrHit["aate1"], aate2 = meanArrHit["aate2"] + 2*stdErrArrHit["aate2"], beta = meanArrLow + 2*stdErrArrLow, lowLET = TRUE)[, 2][2]
      lower = calculateComplexId.te(r = r, L = L, d = c(0, d), aate1 = meanArrHit["aate1"] - 2*stdErrArrHit["aate1"], aate2 = meanArrHit["aate2"] - 2*stdErrArrHit["aate2"], beta = meanArrLow - 2*stdErrArrLow, lowLET = TRUE)[, 2][2]
    }
    return (c(lower, upper))
  }
}

# Parameter initialization
if (mod) {
  mixderCurve = calculateComplexId(r, L, d = d, lowLET = TRUE)
} else {
  mixderCurve = calculateComplexId.te(r, L, d = d, lowLET = TRUE)
}
fourIonMIXDER = data.frame(d = mixderCurve[, 1], CA = mixderCurve[, 2])
numDosePoints = length(fourIonMIXDER$d)
naiveCI = matrix(nrow = 2, ncol = numDosePoints)
monteCarloCI = matrix(nrow = 2, ncol = numDosePoints)

# Calculate CI for each dose point
for (i in 1 : numDosePoints) {
  naiveCI[, i] = Generate_CI(d = fourIonMIXDER$d[i], r = r,  L = L, sampleCurves = curveList)
  monteCarloCI[, i] = Generate_CI(doseIndex = i, r = r,  L = L, method = 1, sampleCurves = curveList)
  print(paste("Iterating on dose points. Currently at step:", toString(i), "Total of", toString(numDosePoints), "steps."))
}

#############
# CI Helper #
#############
CIHelper = function(sampleNum, intervalLength = 0.95, d, r, L, HINmodel = hinm, HITmodel = hitm, LOWmodel = LOW.m, model = mod, seed = 100) {
  # Set the pseudorandom seed
  set.seed(seed)
  
  # Generate N randomly generated samples of parameters of HZE model.
  curveList = Generate_samples(N = sampleNum, r = r, L = L, d = d)
  
  # Parameter initialization
  if (mod) {
    mixderCurve = calculateComplexId(r, L, d = d, lowLET = TRUE)
  } else {
    mixderCurve = calculateComplexId.te(r, L, d = d, lowLET = TRUE)
  }
  fourIonMIXDER = data.frame(d = mixderCurve[, 1], CA = mixderCurve[, 2])
  numDosePoints = length(fourIonMIXDER$d)
  naiveCI = matrix(nrow = 2, ncol = numDosePoints)
  monteCarloCI = matrix(nrow = 2, ncol = numDosePoints)
  
  # Calculate CI for each dose point
  for (i in 1 : numDosePoints) {
    naiveCI[, i] = Generate_CI(d = fourIonMIXDER$d[i], r = r,  L = L, sampleCurves = curveList)
    monteCarloCI[, i] = Generate_CI(doseIndex = i, r = r,  L = L, method = 1, sampleCurves = curveList)
    print(paste("Iterating on dose points. Currently at step:", toString(i), "Total of", toString(numDosePoints), "steps."))
  }
  return (list(naiveCI, monteCarloCI))
}

#==============================================#
#==========Confidence Interval End=============#
#==============================================#