#   Filename: HGsynergyMain.R 
#   Purpose: Concerns radiogenic mouse HG tumorigenesis.

#   Copyright: (C) 2017 Mark Ebert, Edward Huang, Dae Woong Ham, Yimin Lin, and Ray Sachs 

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License version 3 as published 
#   by the Free Software Foundation.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#   Attribution Information: This R script was developed at UC Berkeley.
#   Written by Dae Woong Ham Summer 2017. Additions, corrections, changes, 
#   quality control, reorganization by Edward Huang, Yimin Lin, Mark Ebert 
#   and Ray Sachs Fall UCB semester 2017.

#   Relevant references and abbreviations:
#   ".93Alp" = Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-particle radiations." Rad Res 136:382-391 (1993)
#   ".94Alp" = Alpen et al. "Fluence-based relative biological effectiveness for charged particle carcinogenesis in mouse Harderian gland." Adv Space Res 14(10): 573-581. (1994).  
#   "16Chang" = Chang et al. "Harderian Gland Tumorigenesis: Low-Dose and LET Response." Radiat Res 185(5): 449-460. (2016).  
#   "16Srn" = Siranart et al."Mixed Beam Murine Harderian Gland Tumorigenesis: Predicted Dose-Effect Relationships if neither Synergism nor Antagonism Occurs." Radiat Res 186(6): 577-591 (2016).  
#   "17Cuc" = Cucinotta & Cacao. "Non-Targeted Effects Models Predict Significantly Higher Mars Mission Cancer Risk than Targeted Effects Models." Sci Rep 7(1): 1832. (2017). PMC5431989


#=========================== DEPENDENCIES ==========================#
rm(list=ls()) #   To be removed when script is finalized
library(deSolve) #  Solving differential equations
library(ggplot2) #   Plotting
library(mvtnorm) #   Monte Carlo simulation


#=============================== DATA ==============================#
hg_data <- data.frame( #  Data used in 16Chang; includes data analyzed in .93Alp and .94Alp  
  dose = c(0.2,0.4,0.6,1.2,2.4,3.2,5.1,7,0.05,0.1,0.15,0.2,0.4,0.8,1.6,0.05,0.1,0.2,0.4,0,0.1,0.2,0.4,0.8,1.6,0.4,0.8,1.6,3.2,0.05,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.04,0.08,0.16,0.32,0.033,0.066,0.13,0.26,0.52,.2, .4, .6),
  HG = c(0.091,0.045,0.101,0.169,0.347,0.431,0.667,0.623,0.156,0.215,0.232,0.307,0.325,0.554,0.649,0.123,0.145,0.207,0.31,0.026,0.083,0.25,0.39,0.438,0.424,0.093,0.195,0.302,0.292,0.109,0.054,0.066,0.128,0.286,0.183,0.167,0.396,0.536,0.192,0.234,0.317,0.092,0.131,0.124,0.297,0.082,0.088,0.146,0.236,0.371,.154,.132,.333), #  HG prevalence as defined in 16Chang
  NWeight = 0.01 * c(520,2048,1145,584,313,232,293,221,1162,877,455,409,374,223,320,742,661,347,131,6081,1091,251,244,191,131,645,255,199,111,649,378,973,833,201,468,381,197,109,496,257,185,1902,1063,884,350,1767,1408,874,299,261,322,206,67), #  nominal weight for weighted least squaresregression; see .93Alp. The Lanthanum entries were obtained by measuring the main graph in 17Cuc 
  index = c(rep(1, 8),rep(0,17), rep(1, 4),  rep(0, 24)), #  Index = 0 for Z > 3 ions, 1 otherwise. Not needed in some models
  L = c(rep(1.6,8), rep(193, 7), rep(250, 4), rep(195, 6), rep(0.4, 4), rep(25, 5), rep(464, 4), rep(193, 3),rep(70, 4), rep(100, 5), rep(953, 3)), #  L = LET = LET_infinity = stopping power (keV/micron)
  Z = c(rep(2, 8), rep(26, 17), rep(1, 4), rep(10, 5), rep(43, 4), rep(26, 3), rep(14, 4), rep(22, 5), rep(57, 3)), #  Atomic number, charge in units of proton charge on fully ionized atomic nucleus, e.g. 2 for 2He4
  Zeff = c(rep("TBD", 53)), #  Effective ion charge according to the formula of W.H Barkas. Zeff <= Z. Calculated below. For this data, only very slightly less than Z.
  beta = c(rep("TBD", 53)), #  Ion speed, relative to speed of light, calculated below
  MeVperu = c(rep(228, 8), rep(600, 7), rep(300, 4), rep(600, 6), rep(250, 4), rep(670, 5), rep(600, 4), rep(600, 3), rep(260, 4), rep(1000, 5), rep(593, 3)), #  Kinetic energy in MeV, divided by atomic mass, e.g. divided by 4u = 4 x 931.5 MeV/c^2 for 2He4
  Katz = c(rep("TBD", 53)), #  For fully ionized nuclei, Katz's Z^2/beta^2, Calculated below. It is part of the Bethe Barkas Bloch equation for stopping power. Our calculations don't use Katz, but various similar calculations do.
  ion = c(rep("He4", 8), rep("Fe56", 17), rep("p", 4), rep("Ne20", 5), rep("Nb93", 4), rep("Fe56", 3), rep("Si28", 4), rep("Ti48", 5), rep("La139", 3)),
  comments = c(".93AlpLooksOK", rep("", 7), ".93AlplooksOK", rep("", 11), ".93Alp.no.iso", "not in 17Cuc (or 16Chang?)", rep("", 3), "16Chang all OK?", rep('', 24), ".94Alp","From graphs",'e.g. in 17Cuc')
) 

# Data for HG induced by photons from Cs-137 or Co-60 beta decay; from 16Chang 
beta_decay_data <- data.frame(
  dose = c(0, 0.4, 0.8, 1.6, 3.2, 7, 0, .4, .8, .12, 1.6),
  HG = c(.026, .048, .093, .137, .322, .462, .0497, .054, .067, .128, .202),
  NWeight = c(6081.2, 4989.5, 1896.8, 981.1, 522.2, 205.2, 7474.1, 2877.6, 1423.7, 689.9, 514.9),
  Nucleus = c(rep("Cobalt-60", 6), rep("Cesium-137", 5)),
  Comments = c(rep("TBD", 11))
)

GeVu <- 0.001 * hg_data[, "MeVperu"] #  Convert to GeV/u for convenience in a calculation
hg_data[, "Katz"] <- round(hg_data[, "Z"] ^ 2 * (2.57 * GeVu ^2 + 4.781 * GeVu + 2.233) / (2.57 * GeVu ^ 2 + 4.781 * GeVu), 2) #  Special relativistic calculation of Z^2/beta^2. The numerics include conversion from GeV to joules and from u to kg.
hg_data[, "beta"] <- round(hg_data[, "Z"] * sqrt(1 / hg_data[, "Katz"]), 3) #  i.e. Z * sqrt(beta ^ 2 / Z ^ 2) 
hg_data[, "Zeff"] <- round(hg_data[, "Z"] * (1 - exp( -125 * hg_data[, "Z"] ^ ( - 2.0 / 3))), 2) #  Barkas formula for Zeff; for us Zeff is almost Z

hg_data <- within(hg_data, L[L < 200 & ion == 'Fe56'] <- 185) # Set all Fe56 with L < 200 to L = 185 
clean_hg_data <- hg_data[c(1:19, 21:53), ] # Removes the zero dose case 
clean_HZE_data <- subset(clean_hg_data, Z > 3) #  Look only at HZE not at much lower Z and LET ions. 
clean_light_ion_data <- subset(clean_hg_data, Z <= 3) 
#  NOTE: "Light" refers to ionized atomic nuclei lighter than Beryllium. 
#  The data published to date has such a big gap between alpha particles (Z=2, LET~1.6) 
#  and Neon (Z=10, LET~25) that we will here use independent models for Z<3 and Z>3 


#===================== MISC. OBJECTS & VARIABLES ===================#
# In next line phi controls how fast NTE build up from zero; not really needed during calibration since phi * Dose >> 1 at every observed Dose !=0. phi needed for later synergy calculations.
phi <- 2000 #  even larger phi should give the same final results, but might cause extra problems with R. 

dose_vector <- c(
  seq(0, .00001, by = 0.000001), #  Look carefully near zero, but go out to 0.5 Gy
  seq(.00002, .0001, by = .00001),
  seq(.0002, .001, by = .0001),
  seq(.002, .01, by = .001),
  seq(.02, .5, by = .01))

#####  photon model #####
beta_decay_lm <- lm(HG ~ dose , data = beta_decay_data) #  Linear model fit on beta_decay_data dataset 
summary(beta_decay_lm, correlation = TRUE) 


#========================= HZE/NTE MODEL ===========================#
# (HZE = high charge and energy; NTE = non-targeted effects are included)

# Uses 3 adjustable parameters. 
hi_nte_model <- nls( #  Calibrating parameters in a model that modifies the hazard function NTE models in 17Cuc. 
  HG ~ .0275 + (1 - exp ( -0.01 * (aa1 * L * dose * exp( - aa2 * L) + (1 - exp( - phi * dose)) * kk1))), 
  data = clean_HZE_data, 
  weights = NWeight,
  start = list(aa1 = .9, aa2 = .01, kk1 = 6)) 

summary(hi_nte_model, correlation = T) #  Parameter values & accuracy
vcov(hi_nte_model) #  Variance-covariance matrix RKSB
hi_nte_model_coef <- coef(hi_nte_model) #  Calibrated central values of the 3 parameters. Next is the IDER, = 0 at dose 0

calib_nte_hazard_func <- function(dose, L) { #  Calibrated hazard function 
  0.01 * (hi_nte_model_coef[1] * L * dose * exp( - hi_nte_model_coef[2] * L) + (1 - exp( - phi * dose)) * hi_nte_model_coef[3])
} 

calib_HZE_nte_ider <- function(dose, L) { #  Calibrated HZE NTE IDER
  1 - exp( - calib_nte_hazard_func(dose, L)) 
}


#=========================== HZE/TE MODEL  ==========================#
#  (TE = targeted effects only)

hi_te_model <- nls( #  Calibrating parameters in a TE only model.
  HG ~ .0275 + (1 - exp ( - 0.01 * (aate1 * L * dose * exp( - aate2 * L)))), 
  data = clean_HZE_data,  
  weights = NWeight,
  start = list(aate1 = .9, aate2 = .01)) 

summary(hi_te_model, correlation = T) #  Parameter values & accuracy
vcov(hi_te_model) #  Variance-covariance matrix RKSB
hi_te_model_coef <- coef(hi_te_model) #  Calibrated central values of the 2 parameters. Next is the IDER, = 0 at dose 0

calib_te_hazard_func <- function(dose, L) { #  Calibrated hazard function
  0.01 * (hi_te_model_coef[1] * L * dose * exp( - hi_te_model_coef[2] * L))
} 

calib_HZE_te_ider <- function(dose, L) {
  1 - exp( - calib_te_hazard_func(dose, L)) #  Calibrated HZE TE IDER
}


#============= LIGHT ION, LOW Z (<= 3), LOW LET MODEL ==============#
low_LET_model <- nls(
  HG ~ .0275 + 1 - exp( - bet * dose),
  data = clean_light_ion_data,
  weights = NWeight,
  start = list(bet = .5))

summary(low_LET_model)
low_LET_model_coef <- coef(low_LET_model)  # Calibrated central values of the parameter

calib_low_LET_ider <- function(dose, L) { # Calibrated Low LET model. Use L=0, but maybe later will use L > 0 but small 
  return(1 - exp( - low_LET_model_coef[1] * dose))
}  

low_LET_slope <- function(dose, L) { # Slope dE/dd of the low LET, low Z model; looking at the next plot() it seems fine
  low_LET_model_coef * exp( - low_LET_model_coef * dose)  
}


#========================== VISUAL CHECKS ==========================#
# plot () chunks such as the following are visual check to see if our calibration is consistent with 16Chang, .93Alp, .94Alp
# and 17Cuc; (ggplot commands are Yinmin's and concern CI)
# Put various values in our calibrated model to check with numbers and graphs in these references
plot(c(0, 7), c(0, 1), col = 'red', ann = 'F') 
lines(0.01 * 0:700, calib_low_LET_ider(0.01 * 0:700, 0) + .0275)  #  calibrated lowLET IDER
points(clean_light_ion_data[1:8, "dose"], clean_light_ion_data[1:8, "HG"], pch = 19) #  RKS: Helium data points
points(clean_light_ion_data[9:12, "dose"], clean_light_ion_data[9:12, "HG"] )  #  proton data points 


#======================= INFORMATION CRITERION =====================#
info_crit_table <- cbind(AIC(hi_te_model, hi_nte_model), BIC(hi_te_model, hi_nte_model))
print(info_crit_table)

#======= SEA Calculator ======#
calculate_SEA <- function(total_dose, r, L) {
  total = 0
  i = 1
  if (length(r) == 1) {
    while (i < r + 1) {
      total = total + calib_HZE_nte_ider(total_dose * r, L[i])
      i = i + 1
    }
  } else {
    while (i < length(r) + 1) {
      total = total + calib_HZE_nte_ider(total_dose * r[i], L[i])
      i = i + 1
    }
  }
  return(total)
}


#======= I(d) CALCULATOR; HZE NTE/TE MODELS; OPTIONAL LOW-LET ======#
calculate_complex_id <- function(r, L, d, lowLET = FALSE, model = "NTE",
                                 coef = list(NTE = hi_nte_model_coef, TE = hi_te_model_coef, lowLET = low_LET_model_coef),
                                 iders = list(NTE = calib_HZE_nte_ider, TE = calib_HZE_te_ider, lowLET = calib_low_LET_ider),
                                 calculate_dI = c(NTE = .calculate_dI_nte, TE = .calculate_dI_te),
                                 phi = 2000) {
  dE <- function(yini, state, pars) { #  Constructing an ode from the IDERS
    with(as.list(c(state, pars)), {
      aa <- u <- dI <- vector(length = length(L))
      for (i in 1:length(L)) {
        aa[i] <- pars[1] * L[i] * exp(-pars[2] * L[i])
        u[i] <- uniroot(function(d) HZE_ider(d, L[i]) - I, 
                        interval = c(0, 200), 
                        extendInt = "yes", 
                        tol = 10 ^ - 10)$root
        dI[i] <- r[i] * calc_dI(aa[i], u[i], pars[3])
      }
      if (lowLET == TRUE) { # If low-LET IDER is present then include it at the end of the dI vector
        u[length(L) + 1] <- uniroot(function(d) calib_low_LET_ider(d, coef["lowLET"]) - I, 
                                    interval = c(0, 200), 
                                    extendInt = "yes", 
                                    tol = 10 ^ - 10)$root
        dI[length(L) + 1] <- r[length(r)] * low_LET_slope(d = u[length(L) + 1], L = 0)
      }
      return(list(sum(dI)))
    })
  }
  p <- list(pars = coef[[model]], 
            HZE_ider = iders[[model]], 
            calib_low_LET_ider = iders[["lowLET"]], 
            calc_dI = calculate_dI[[model]])
  return(ode(c(I = 0), times = d, dE, parms = p))
}


#=================== dI HIDDEN FUNCTIONS =====================#
.calculate_dI_nte <- function(aa, u, kk1) {
  return(0.01 * (aa + exp( - phi * u) * kk1 * phi) * exp( - 0.01 * (aa * u + (1 -exp( - phi * u)) * kk1)))
}

.calculate_dI_te <- function(aa, u, pars = NULL) {
  return(0.01 * aa * exp(-0.01 * aa * u))
}


#========================== PLOTS ============================#
# Plot 1 : one HZE one low-LET; NTE
d <- .01 * 0:300.
r1 <- .2
r <- c(r1, 1 - r1) #Proportions. Next plot IDERs and MIXDER
plot(x = d, y = calib_HZE_nte_ider(dose = d, L = 173), type = "l", xlab = "dose", ylab = "HG", bty = 'l', col = 'green', lwd = 2)
lines(x = d, y = calib_low_LET_ider(d, 0), col = 'green', lwd = 2)
lines(x = d, y = calculate_complex_id(r = r, L = 193, d = d, lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)

# Plot 2 : one HZE one low-LET, TE & NTE
d <- .01 * 0:300
r1 <- .2
r <- c(r1, 1 - r1) #Proportions. Next plot IDERs and MIXDER
plot(x = d, y = calib_HZE_te_ider(dose = d, L = 173), type = "l", xlab = "dose", ylab = "HG", bty = 'l', col = 'green', lwd = 2)
lines(x = d, y = calib_low_LET_ider(d, 0), col = 'green', lwd = 2)
lines(x = d, y = calculate_complex_id(r = r, L = 193, d = d, model = "TE", lowLET = TRUE)[, 2], col = "orange", lwd = 2) # I(d)
lines(x = d, y = calculate_complex_id(r = r, L = 193, d = d, model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2)

# Plot 3: four HZE; NTE
plot(calculate_complex_id(r = rep(0.25, 4), L = c(25, 70, 190, 250), d = dose_vector), type = 'l', col = 'red', bty = 'l', ann = 'F') #  I(d) plot
SEA <- function(dose) {
  return(calib_HZE_nte_ider(dose / 4, 25) + 
         calib_HZE_nte_ider(dose / 4, 70) + 
         calib_HZE_nte_ider(dose / 4, 190) + 
         calib_HZE_nte_ider(dose / 4, 250))
}
lines(dose_vector, SEA(dose_vector), lty = 2)
lines(dose_vector, calib_HZE_nte_ider(dose_vector, 190), col = 'green') # component 4
lines(dose_vector, calib_HZE_nte_ider(dose_vector, 250), col = 'green') # component 3
lines(dose_vector, calib_HZE_nte_ider(dose_vector, 70), col = 'green') # component 2
lines(dose_vector, calib_HZE_nte_ider(dose_vector, 25), col = 'green') # component 1


# Plot 4: two HZE; NTE; one low-LET
d <- seq(0, .01, .0005)
plot(x = d, y = calib_HZE_nte_ider(dose = d, L = 173), type = "l", xlab = "dose", ylab = "HG", bty = 'l', col = 'green', lwd = 2)
lines(x = d, y = calib_HZE_nte_ider(d, 70), col = 'green', lwd = 2) # component 3
lines(x = d, y = calib_low_LET_ider(d, 0), col = 'green', lwd = 2)
lines(x = d, y = calculate_complex_id(r = c(1/20, 1/20, 9/10), L = c(70, 173), d = d, lowLET = TRUE)[, 2], col = 'red', lwd = 2)


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
forty_cGy <- .01 * 0:40
sixty_cGy <- .01 * 0:60
seventy_cGy <- .01 * 0:70
hundred_cGy <- .01 * 0:100
forty_nine_cGy <- .01 * 0:49

# Fig. 3.2.1.1. - Fe56 (600 MeV/u), Si28, and corresponding IEA and SEA MIXDERS.
setEPS()
postscript("fe56_si28_nte.eps")
plot(x = forty_cGy, y = calculate_SEA(forty_cGy, 1/2, c(70, 195)), type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2, xaxs="i")
lines(x = forty_cGy, y = calib_HZE_nte_ider(dose = forty_cGy, L = 70), col = "cyan", lwd = 2)
lines(x = forty_cGy, y = calib_HZE_nte_ider(dose = forty_cGy, L = 195), col = "darkcyan", lwd = 2)
lines(x = forty_cGy, y = calculate_complex_id(r = c(0.5 , 0.5), L = c(70, 195), d = forty_cGy, model = "NTE", lowLET = FALSE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 0.4, lwd = 2)
legend(x = "topleft", legend = c("Fe56 (600 MeV/u), NTE NTE-TE IDER", "Si28 HZE NTE-TE IDER", "IEA MIXDER (50% Fe56, 50% Si28)", "SEA MIXDER (50% Fe56, 50% Si28)"),
       col = c("darkcyan", "cyan", "red", "black"), lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 2), inset = 0.05)
dev.off()


# Fig. for a mixture of 60 cGy H1 (protons) with 40 cGy Si;
setEPS()
postscript("h1_si28_nte.eps")
plot(x = hundred_cGy, y = calculate_SEA(hundred_cGy, c(0.6, .4), c(0.4, 70)), type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2, xaxs="i")
lines(x = hundred_cGy, y = calib_low_LET_ider(dose = hundred_cGy, L = 0.4), col = "cyan", lwd = 2)
lines(x = hundred_cGy, y = calib_HZE_nte_ider(dose = hundred_cGy, L = 70), col = "darkcyan", lwd = 2)
lines(x = hundred_cGy, y = calculate_complex_id(r = c(0.6 , 0.4), L = c(0.4, 70), d = hundred_cGy, model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 1, lwd = 2)
legend(x = "topleft", legend = c("Si28 HZE NTE-TE IDER", "H1 Low-LET NTE-TE IDER","IEA MIXDER (60% H1, 40% Si28)", "SEA MIXDER (60% H1, 40% Si28)"),
       col = c("darkcyan", "cyan", "red", "black"), lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 2), inset = 0.025)
dev.off()

# Fig. for a mixture of 40 cGy H1 with 30 cGy Fe56 at 600 MeV/u;
setEPS()
postscript("h1_fe56_nte.eps")
plot(x = seventy_cGy, y = calculate_SEA(seventy_cGy, r = c(4/7, 3/7), c(0.4, 195)), type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2, xaxs="i")
lines(x = seventy_cGy, y = calib_low_LET_ider(dose = seventy_cGy, L = 0.4), col = "cyan", lwd = 2)
lines(x = seventy_cGy, y = calib_HZE_nte_ider(dose = seventy_cGy, L = 195), col = "darkcyan", lwd = 2)
lines(x = seventy_cGy, y = calculate_complex_id(r = c(4/7 , 3/7), L = c(0.4, 195), d = seventy_cGy, model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 0.7, lwd = 2)
legend(x = "topleft", legend = c("Fe56 (600 MeV/u), HZE NTE-TE IDER", "H1 Low-LET NTE-TE IDER","IEA MIXDER (57% H1, 43% Si28)", "SEA MIXDER (57% H1, 43% Si28)"),
       col = c("darkcyan", "cyan", "red", "black"), lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 2), inset = 0.005)
dev.off()

# Fig. for a mixture of all 7 HZE ions; total dose 49 cGy, each ion gets 7 cGy;
setEPS()
postscript("all_hze_nte.eps")

# assume Hi HZE implies Z > 3
plot(x = forty_nine_cGy, y = calculate_SEA(forty_nine_cGy, r =  c(1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7), c(25, 70, 100, 195, 250, 464, 953)), 
     type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2, axes=FALSE)

lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 25), col = "pink", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 70), col = "orange", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 100), col = "yellow", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 195), col = "green", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 250), col = "blue", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 464), col = "purple", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 953), col = "violet", lwd = 2)

lines(x = forty_nine_cGy, y = calculate_complex_id(r = c(1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7), L = c(25, 70, 100, 195, 250, 464, 953),
                                                d = forty_nine_cGy, model = "NTE", lowLET = FALSE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 0.49, lwd = 1)
axis(2, c(-0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
axis(1, c(-.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.49), xaxs = "i")
legend(x = "topleft", legend = c("Ne20 NTE-TE IDER", "Si28 NTE-TE IDER", 
                                 "Ti48 NTE-TE IDER", "Fe56 (600 MeV/u) NTE-TE IDER", 
                                 "Fe56 (300 MeV/u) NTE-TE IDER", "Nb93 NTE-TE IDER",
                                 "La139 NTE-TE IDER",
                                 "IEA MIXDER (Equally Distributed)", "SEA MIXDER (Equally Distributed)"),
       col = c("pink", "orange", "yellow", "green", "blue", "purple", "violet", "red", "black"), 
      lwd = c(2, 2, 2, 2, 2, 2, 2, 2, 2), 
      lty = c(1, 1, 1, 1, 1, 1, 1, 1, 2), cex = 0.55, inset = 0.0125)
dev.off()



# h1, 0.4, 60
# he4, 1.6, 20
# o16, 25, 10
# si28, 70, 2.5
# ti28, 100, 2.5
# fe56, 195, 5
setEPS()
postscript("big_mix_nte.eps")
plot(x = forty_nine_cGy, y = calculate_SEA(forty_nine_cGy, r = c(.6, .2, .1, .025, .025, .5), L = c(0.4, 1.6, 25, 70, 100, 195)), 
     type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2, xaxs = "i")

lines(x = forty_nine_cGy, y = calib_low_LET_ider(dose = forty_nine_cGy, L = 0.4), col = "orange", lwd = 2)
lines(x = forty_nine_cGy, y = calib_low_LET_ider(dose = forty_nine_cGy, L = 1.6), col = "yellow", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 25), col = "green", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 70), col = "blue", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 100), col = "purple", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 195), col = "violet", lwd = 2)

lines(x = forty_nine_cGy, y = calculate_complex_id(r = c(.6, .2, .1, .025, .025, .5), L =  c(0.4, 1.6, 25, 70, 100, 195),
                                                   d = forty_nine_cGy, model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 1, lwd = 1)
legend(x = "right", legend = c("H1 Low-LET NTE-TE IDER", "He4 Low-LET NTE-TE IDER", "O16 NTE-TE IDER", 
                                 "Si28 NTE-TE IDER", "Ti48 NTE-TE IDER", "Fe56 (600 MeV/u) NTE-TE IDER", 
                                 "IEA MIXDER (Equally Distributed)", "SEA MIXDER (Equally Distributed)"),
       col = c("orange", "yellow", "green", "blue", "purple", "violet", "red", "black"), 
       lwd = c(2, 2, 2, 2, 2, 2, 2, 2), 
       lty = c(1, 1, 1, 1, 1, 1, 1, 2), cex = 0.7, inset = 0.0125)
dev.off()

