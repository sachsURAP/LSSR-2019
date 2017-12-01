#   Filename: HGsynergyMain.R 
#   Purpose: Concerns radiogenic mouse HG tumorigenesis.

#   Copyright: (C) 2017 Mark Ebert, Edward Huang, Dae Woong Ham, Yimin Lin, and Ray Sachs #Edward: why do we need this? 
#   #Ray: The gnu website (https://www.gnu.org/licenses/gpl-howto.html) 
#         recommends a copyright notice and license notice at the top of our 
#         files. I do not feel particularly strongly about keeping them but I 
#         have seen license/copyright notices in other published scripts
#         and feel that it adds to the professional appearance of ours.

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
clean_hze_data <- subset(clean_hg_data, Z > 3) #  Look only at HZE not at much lower Z and LET ions. 
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
  data = clean_hze_data, 
  weights = NWeight,
  start = list(aa1 = .9, aa2 = .01, kk1 = 6)) 

summary(hi_nte_model, correlation = T) #  Parameter values & accuracy
vcov(hi_nte_model) #  Variance-covariance matrix RKSB
hi_nte_model_coef <- coef(hi_nte_model) #  Calibrated central values of the 3 parameters. Next is the IDER, = 0 at dose 0

calib_nte_hazard_func <- function(dose, L) { #  Calibrated hazard function 
  0.01 * (hi_nte_model_coef[1] * L * dose * exp( - hi_nte_model_coef[2] * L) + (1 - exp( - phi * dose)) * hi_nte_model_coef[3])
} 

calib_hze_nte_ider <- function(dose, L) { #  Calibrated HZE NTE IDER
  1 - exp( - calib_nte_hazard_func(dose, L)) 
}


#=========================== HZE/TE MODEL  ==========================#
#  (TE = targeted effects only)

hi_te_model <- nls( #  Calibrating parameters in a TE only model.
  HG ~ .0275 + (1 - exp ( - 0.01 * (aate1 * L * dose * exp( - aate2 * L)))), 
  data = clean_hze_data,  
  weights = NWeight,
  start = list(aate1 = .9, aate2 = .01)) 

summary(hi_te_model, correlation = T) #  Parameter values & accuracy
vcov(hi_te_model) #  Variance-covariance matrix RKSB
hi_te_model_coef <- coef(hi_te_model) #  Calibrated central values of the 2 parameters. Next is the IDER, = 0 at dose 0

calib_te_hazard_func <- function(dose, L) { #  Calibrated hazard function
  0.01 * (hi_te_model_coef[1] * L * dose * exp( - hi_te_model_coef[2] * L))
} 

calib_hze_te_ider <- function(dose, L) {
  1 - exp( - calib_te_hazard_func(dose, L)) #  Calibrated HZE TE IDER
}


#============= LIGHT ION, LOW Z (<= 3), LOW LET MODEL ==============#
low_let_model <- nls(
  HG ~ .0275 + 1 - exp( - bet * dose),
  data = clean_light_ion_data,
  weights = NWeight,
  start = list(bet = .5))

summary(low_let_model)
low_let_model_coef <- coef(low_let_model)  # Calibrated central values of the parameter

calib_low_let_ider <- function(dose, L) { # Calibrated Low LET model. Use L=0, but maybe later will use L > 0 but small 
  return(1 - exp( - low_let_model_coef[1] * dose))
}  

low_let_slope <- function(dose, L) { # Slope dE/dd of the low LET, low Z model; looking at the next plot() it seems fine
  low_let_model_coef * exp( - low_let_model_coef * dose)  
}

#======================= INFORMATION CRITERION =====================#

AIC_function = function(RSS, k, n) {
  n + n*log(2*pi) + n*log(RSS/n) + 2*(k+1)
}

BIC_function = function(RSS, k, n) {
  n + n*log(2*pi) + n*log(RSS/n) + log(n)*(k+1)
}

calculate_AIC_BIC = function(numDataHZE = length(clean_hze_data$HG), backgroundConst = .0275, NTE_coeff = summary(hi_nte_model)$coefficients[, "Estimate"], TE_coeff  = summary(hi_te_model)$coefficients[, "Estimate"]) {
  
  NTE_function = function(d, L, aa1, aa2, kk1, y0 = backgroundConst) {
    y0 + (1 - exp ( -0.01 * (aa1 * L * d * exp( - aa2 * L) + (1 - exp( - phi * d)) * kk1)))
  }
  
  TE_function = function(d, L, aate1, aate2, y0 = backgroundConst) {
    y0 + (1 - exp ( - 0.01 * (aate1 * L * d * exp( - aate2 * L))))
  }
  
  residualSquareNTE = vector(length = 0)
  for (i in 1:numDataHZE) {
    residualSquareNTE = c(residualSquareNTE, clean_hze_data$NWeight[i] * (clean_hze_data$HG[i] - backgroundConst - NTE_function(d   = clean_hze_data$d[i],
                                                                                                                                L   = clean_hze_data$L[i],
                                                                                                                                aa1 = NTE_coeff[["aa1"]],
                                                                                                                                aa2 = NTE_coeff[["aa2"]],
                                                                                                                                kk1 = NTE_coeff[["kk1"]]))^2)
  }
  
  residualSquareTE = vector(length = 0)
  for (i in 1:numDataHZE) {
    residualSquareTE = c(residualSquareTE, clean_hze_data$NWeight[i] * (clean_hze_data$HG[i] - backgroundConst - TE_function(d     = clean_hze_data$d[i],
                                                                                                                             L     = clean_hze_data$L[i],
                                                                                                                             aate1 = TE_coeff[["aate1"]],
                                                                                                                             aate2 = TE_coeff[["aate2"]] ))^2)
  }
  
  weightedResidualSquareSumNTE = sum(residualSquareNTE)
  weightedResidualSquareSumTE  = sum(residualSquareTE)
  
  NTE_AIC = AIC_function(RSS = weightedResidualSquareSumNTE, k = 3, n = numDataHZE)
  NTE_BIC = BIC_function(RSS = weightedResidualSquareSumNTE, k = 3, n = numDataHZE)
  TE_AIC  = AIC_function(RSS = weightedResidualSquareSumTE , k = 2, n = numDataHZE)
  TE_BIC  = BIC_function(RSS = weightedResidualSquareSumTE , k = 2, n = numDataHZE)
  return (data.frame(AIC_NTE = NTE_AIC, AIC_TE = TE_AIC, BIC_NTE = NTE_BIC, BIC_TE = TE_BIC))
}

print(calculate_AIC_BIC(backgroundConst = .0275))
print(calculate_AIC_BIC(backgroundConst = .05))
print(calculate_AIC_BIC(backgroundConst = .01))

for (i in 1:100) {
  print(.001*i)
  print(calculate_AIC_BIC(backgroundConst = .001*i))
}

