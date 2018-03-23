#   Filename: synergyTheory.R 
#   Purpose: Concerns radiogenic mouse HG tumorigenesis. Contains the NTE-TE,
#            TE, HZE, low-LET, MIXDER, and other relevant synergy theory models.

#   Copyright: (C) 2017 Mark Ebert, Edward Huang, Dae Woong Ham, Yimin Lin, and Ray Sachs 

source("hgData.R") # load in the data


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
#' calculate_SEA(forty_cGy, r = (1/2, 1/2), c(70, 195), n = 2)
#' calculate_SEA(seventy_cGy, r = c(4/7, 3/7), c(0.4, 195))
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


#======= I(d) CALCULATOR; HZE NTE/TE MODELS; OPTIONAL LOW-LET ======#
#' Applies Incremental Effect Additivity 
#' 
#' @description 
#' @param d number corresponding to the sum dose in cGy.
#' @param r Vector of dose ratios, must be length n.
#' @param L Vector of LET values, must be length n.
#' @return The estimate Harderian Gland prevalence from a SEA MIXDER constructed
#'         from the given IDER parameters. 
#' @examples
#' 
#' 
calculate_complex_id <- function(r, L, d, lowLET = FALSE, model = "NTE",
                                 coef = list(NTE = hi_nte_model_coef,
                                             TE = hi_te_model_coef, 
                                             lowLET = low_LET_model_coef),
                                 iders = list(NTE = calib_HZE_nte_ider, 
                                              TE = calib_HZE_te_ider, 
                                              lowLET = calib_low_LET_ider),
                                 calculate_dI = c(NTE = .calculate_dI_nte, 
                                                  TE = .calculate_dI_te),
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
