#   Filename: synergyTheory.R 
#   Purpose: Concerns radiogenic mouse HG tumorigenesis. Contains the NTE-TE,
#            TE, HZE, low-LET, MIXDER, and other relevant synergy theory models.

#   Copyright: (C) 2017 Mark Ebert, Edward Huang, Dae Woong Ham, Yimin Lin, 
#                       Yunzhi Zhang, and Ray Sachs 

source("hgData.R") # load in the data. Dose in Gy (RKS: ASAP we should convert all cGy items to Gy); LET usually in keV/micron; prevalence Prev always < 1 (i.e. not in %, which would mean prevalence < 100 but is strongly deprecated).

library(deSolve) #  Solving differential equations
#===================== MISC. OBJECTS & VARIABLES ===================#
# In next line phi controls how fast NTE build up from zero; not really needed 
# during calibration since phi * Dose >> 1 at every observed Dose !=0. 
# phi needed for later synergy calculations.

phi <- 2000 #  even larger phi should give the same final results, but might cause extra problems with R. 

#======= PHOTON MODEL ===========#
#  Linear model fit on beta_decay_data dataset. We will never recalculate this unless new data comes in but here it is just in case.
#beta_decay_lm <- lm(HG ~ dose + I(dose ^ 2), data = beta_decay_data) 
#summary(beta_decay_lm, correlation = TRUE)

#======================= HZE/NTE MODEL ===================#
# (HZE = high charge and energy; NTE = non-targeted effects are included)
HZE_data <- ion_data[13:47,] # Includes 1-ion data iff Z > 3
# Uses 3 adjustable parameters. 
HZE_nte_model <- nls( #  Calibrating parameters in a model that modifies the hazard function NTE models in 17Cuc. 
  Prev ~ .0275 + (1 - exp ( - (aa1 * LET * dose * exp( - aa2 * LET) + (1 - exp( - phi * dose)) * kk1))), 
  data = HZE_data, 
  weights = NWeight,
  start = list(aa1 = .00009, aa2 = .001, kk1 = .06)) # use extra argument trace=TRUE if you want to watch convergence. 

summary(HZE_nte_model, correlation = TRUE) #  Parameter values & accuracy
vcov(HZE_nte_model) #  Variance-covariance matrix RKSB
HZE_nte_model_coef <- coef(HZE_nte_model) #  Calibrated central values of the 3 parameters. Next is the IDER, = 0 at dose 0

calib_nte_hazard_func <- function(dose, LET, coef) { #  Calibrated hazard function 
 return(coef[1] * LET * dose * exp( - coef[2] * LET) + (1 - exp( - phi * dose)) * coef[3])
} 

calib_HZE_nte_ider <- function(dose, LET, coef = HZE_nte_model_coef) { #  Calibrated HZE NTE IDER
  return(1 - exp( - calib_nte_hazard_func(dose, LET, coef)))
}


#=================== HZE/TE MODEL  ====================#
#  (TE = targeted effects only). RKS. This chunk runs and gives good results. We will not use it in the minor paper, only for later papers.

HZE_te_model <- nls( #  Calibrating parameters in a TE only model.
  Prev ~ .0275 + (1 - exp ( - (aate1 * LET * dose * exp( - aate2 * LET)))),
  data = HZE_data,
  weights = NWeight,
  start = list(aate1 = .00009, aate2 = .01))

summary(HZE_te_model, correlation = TRUE) #  Parameter values & accuracy
vcov(HZE_te_model) #  Variance-covariance matrix RKSB
HZE_te_model_coef <- coef(HZE_te_model) #  Calibrated central values of the 2 parameters. Next is the IDER, = 0 at dose 0

calib_te_hazard_func <- function(dose, LET, coef) { #  Calibrated hazard function
  return(coef[1] * LET * dose * exp( - coef[2] * LET))
}

calib_HZE_te_ider <- function(dose, LET, coef = HZE_te_model_coef) {
  return(1 - exp( - calib_te_hazard_func(dose, LET, coef))) #  Calibrated HZE TE IDER
}


#====== LIGHT ION, LOW Z (<= 3), LOW LET MODEL =========#
low_LET_data = ion_data[1:12, ] #swift protons and alpha particles
low_LET_model <- nls(
  Prev~ .0275 + 1 - exp( - bet * dose),
  data = low_LET_data,
  weights = NWeight,
  start = list(bet = .005))

summary(low_LET_model, correlation = TRUE)
low_LET_model_coef <- coef(low_LET_model)  # Calibrated central values of the parameter

calib_low_LET_ider <- function(dose, LET, beta = low_LET_model_coef[1]) { # Calibrated Low LET model. Use L=0, but maybe later will use L > 0 but small 
  return(1 - exp( - beta * dose))
}  

low_LET_slope <- function(dose, LET) { # Slope dE/dd of the low LET, low Z model; looking at the next plot() it seems fine
  low_LET_model_coef * exp( - low_LET_model_coef * dose)  
}


#========================== VISUAL CHECKS ==========================#
# plot () chunks such as the following are visual check to see if our 
# calibration is consistent with 16Chang, .93Alp, .94Alp
# and 17Cuc; (ggplot commands are Yinmin's and concern CI)
# Put various values in our calibrated model to check with numbers and 
# graphs in these references
plot(c(0, 7), c(0, 1), col = 'red', ann = 'F') 
lines(0.01 * 0:700, calib_low_LET_ider(0:700, 0) + .0275)  #  calibrated lowLET IDER
points(low_LET_data[1:8, "dose"]/100, low_LET_data[1:8, "Prev"], pch = 19) #  RKS: Helium data points
points(low_LET_data[9:12, "dose"]/100, low_LET_data[9:12, "Prev"] )  #  proton data points 


#======================= INFORMATION CRITERION =====================#
info_crit_table <- cbind(AIC(HZE_te_model, HZE_nte_model), BIC(HZE_te_model, HZE_nte_model))
print(info_crit_table)


#' @description Applies Simple Effect Additivity to a MIXDER.
#' @param total_dose Numeric vector corresponding to the sum dose in cGy.
#' @param ratios Numeric vector of all dose ratios, must be length n.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param lowLET Boolean of whether an LET IDER should be included in the MIXDER
#' @param n Number of IDERs, optional argument used to check parameter validity.
#' @details Corresponding elements of ratios, LET should be associated with the
#'          same IDER.
#' @return Numeric vector representing the estimated Harderian Gland 
#'         prevalence from a SEA MIXDER constructed from the given IDER 
#'         parameters. 
#' @examples
#' calculate_SEA(.01 * 0:40, r = c(1/2, 1/2), c(70, 195), n = 2)
#' calculate_SEA(.01 * 0:70, r = c(4/7, 3/7), c(0.4, 195))
#'
calculate_SEA <- function(total_dose, ratios, LET, lowLET = FALSE, n = NULL) {
  if (!is.null(n) && (n != length(ratios) | n != length(LET))) {
    stop("Length of arguments do not match.") 
  } else if (sum(ratios) != 1) {
    stop("Sum of ratios do not add up to one.")
  } #  End error handling
  total = 0
  i = 1
  if (lowLET == TRUE) { #  First elements of ratos and LET should be the low-LET IDER
    total = total + calib_low_LET_ider(total_dose * ratios[i], LET[i])
    i = i + 1
  } 
  while (i < length(ratios) + 1) { #  Iterate over HZE ions in MIXDER
    total = total + calib_HZE_nte_ider(total_dose * ratios[i], LET[i])
    i = i + 1
  }
  return(total)
}


#' @description Applies Incremental Effect Additivity to a MIXDER.
#' @param r Numeric vector of all dose ratios, must be length n.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param d Numeric vector corresponding to the sum dose in cGy.
#' @param lowLET Boolean of whether an LET IDER should be included in the MIXDER.
#' @param model String value corresponding to the model to be used, either 
#'              "NTE" or "TE". 
#' @param coef Named list of numeric vectors containing coefficients for IDERs.
#' @param iders Named list of functions containing relevant IDER models.
#' @param calculate_dI Named vector of functions to calculate dI depending on 
#'                     the selected model.
#' @param phi Numeric value, used in NTE models.
#' @details Corresponding elements of ratios, LET should be associated with the
#'          same IDER.
#' @return Numeric vector representing the estimated Harderian Gland 
#'         prevalence from an IEA MIXDER constructed from the given IDER 
#'         parameters. 
#' @examples
#' calculate_complex_id(d = .01 * 0:40, r = c(1/2, 1/2), L = c(70, 195))
#' calculate_complex_id(d = .01 * 0:70, r = c(4/7, 3/7), L = c(0.4, 195), 
#'                      lowLET = TRUE, model = "TE")
#'
calculate_complex_id <- function(r, LET, d, lowLET = FALSE, model = "NTE",
                                 coef = list(NTE = HZE_nte_model_coef, # [1] = aa1, [2] = aa2, [3] == kk1
                                             TE = HZE_te_model_coef, 
                                             lowLET = low_LET_model_coef),
                                 iders = list(NTE = calib_HZE_nte_ider, 
                                              TE = calib_HZE_te_ider, 
                                              lowLET = calib_low_LET_ider),
                                 calculate_dI = c(NTE = .calculate_dI_nte, 
                                                  TE = .calculate_dI_te),
                                 phi = 2000) {
  dE <- function(yini, state, pars) { #  Constructing an ode from the IDERS
    with(as.list(c(state, pars)), {
      aa <- u <- dI <- vector(length = length(LET))
      for (i in 1:length(LET)) {
        aa[i] <- pars[1] * LET[i] * exp( - pars[2] * LET[i])
        u[i] <- uniroot(function(d) HZE_ider(d, LET[i], pars) - I, 
                        interval = c(0, 20000), 
                        extendInt = "yes",
                        tol = 10 ^ - 10)$root
        dI[i] <- r[i] * calc_dI(aa[i], u[i], pars[3])
      }
      if (lowLET == TRUE) { # If low-LET IDER is present then include it at the end of the dI vector
        u[length(LET) + 1] <- uniroot(function(d) calib_low_LET_ider(dose = d, LET = LET, 
                                      beta = coef[["lowLET"]]) - I, 
                                      interval = c(0, 20000), 
                                      extendInt = "yes", 
                                      tol = 10 ^ - 10)$root
        dI[length(LET) + 1] <- r[length(r)] * low_LET_slope(d = u[length(LET) + 1], LET = 0)
      }
      return(list(sum(dI)))
    })
  }
  p <- list(pars = coef[[model]], 
            HZE_ider = iders[[model]], 
            calib_low_LET_ider = iders[["lowLET"]], 
            calc_dI = calculate_dI[[model]])
  return(ode(c(I = 0), times = d, dE, parms = p, method = "radau"))
}

#=================== dI HIDDEN FUNCTIONS =====================#
.calculate_dI_nte <- function(aa, u, kk1) {
  return((aa + exp( - phi * u) * kk1 * phi) * exp( - (aa * u + (1 -exp( - phi * u)) * kk1)))
}

.calculate_dI_te <- function(aa, u, pars = NULL) {
  return(aa * exp(- aa * u))
}