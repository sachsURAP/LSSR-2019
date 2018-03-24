#   Filename: HGsynergyMain_merge2.R 
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
#   < authors and contributions to be added later > 

#   Relevant references and abbreviations:
#   ".93Alp" = Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-particle radiations." Rad Res 136:382-391 (1993)
#   ".94Alp" = Alpen et al. "Fluence-based relative biological effectiveness for charged particle carcinogenesis in mouse Harderian gland." Adv Space Res 14(10): 573-581. (1994).  
#   "16Chang" = Chang et al. "Harderian Gland Tumorigenesis: Low-Dose and LET Response." Radiat Res 185(5): 449-460. (2016).  
#   "16Srn" = Siranart et al."Mixed Beam Murine Harderian Gland Tumorigenesis: Predicted Dose-Effect Relationships if neither Synergism nor Antagonism Occurs." Radiat Res 186(6): 577-591 (2016).  
#   "17Cuc" = Cucinotta & Cacao. "Non-Targeted Effects Models Predict Significantly Higher Mars Mission Cancer Risk than Targeted Effects Models." Sci Rep 7(1): 1832. (2017). PMC5431989

library(deSolve) #  solving differential equations
library(ggplot2) #  plotting
library(mvtnorm) #  Monte Carlo simulation
library(minpack.lm) #  non-linear regression 
rm(list=ls())
#=========================== DATA START ===========================#
dfr <- data.frame( #  data used in 16Chang; includes data analyzed in .93Alp and .94Alp  
  dose.1 = c(0.2,0.4,0.6,1.2,2.4,3.2,5.1,7,0.05,0.1,0.15,0.2,0.4,0.8,1.6,0.05,0.1,0.2,0.4,0,0.1,0.2,0.4,0.8,1.6,0.4,0.8,1.6,3.2,0.05,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.04,0.08,0.16,0.32,0.033,0.066,0.13,0.26,0.52,.2, .4, .6),
  HG = c(0.091,0.045,0.101,0.169,0.347,0.431,0.667,0.623,0.156,0.215,0.232,0.307,0.325,0.554,0.649,0.123,0.145,0.207,0.31,0.026,0.083,0.25,0.39,0.438,0.424,0.093,0.195,0.302,0.292,0.109,0.054,0.066,0.128,0.286,0.183,0.167,0.396,0.536,0.192,0.234,0.317,0.092,0.131,0.124,0.297,0.082,0.088,0.146,0.236,0.371,.154,.132,.333), #  HG prevalence as defined in 16Chang
  NWeight = c(520,2048,1145,584,313,232,293,221,1162,877,455,409,374,223,320,742,661,347,131,6081,1091,251,244,191,131,645,255,199,111,649,378,973,833,201,468,381,197,109,496,257,185,1902,1063,884,350,1767,1408,874,299,261,322,206,67), #  nominal weight for weighted least squaresregression; see .93Alp. The Lanthanum entries were obtained by measuring the main graph in 17Cuc 
  index=c(rep(1,8),rep(0,17), rep(1,4),  rep(0,24)), #  index=0 for Z>3 ions, 1 otherwise. Not needed in some models
  L = c(rep(1.6,8), rep(193, 7), rep(250, 4), rep(195, 6), rep(0.4, 4), rep(25, 5), rep(464, 4), rep(193, 3),rep(70, 4), rep(100, 5), rep(953, 3)), #  L = LET = LET_infinity = stopping power (keV/micron)
  Z = c(rep(2, 8), rep(26, 17), rep(1, 4), rep(10, 5), rep(43, 4), rep(26, 3), rep(14, 4), rep(22, 5), rep(57, 3)), #  atomic number, charge in units of proton charge on fully ionized atomic nucleus, e.g. 2 for 2He4
  Zeff = c(rep("TBD", 53)), #  effective ion charge according to the formula of W.H Barkas. Zeff <= Z. Calculated below. For this data, only very slightly less than Z.
  beta = c(rep("TBD", 53)), #  ion speed, relative to speed of light, calculated below
  MeVperu = c(rep(228, 8), rep(600, 7), rep(300, 4), rep(600, 6), rep(250, 4), rep(670, 5), rep(600, 4), rep(600, 3), rep(260, 4), rep(1000, 5), rep(593, 3)), #  Kinetic energy in MeV, divided by atomic mass, e.g. divided by 4u=4x931.5 MeV/c^2 for 2He4
  Katz = c(rep("TBD", 53)), #  for fully ionized nuclei, Katz's Z^2/beta^2, Calculated below. It is part of the Bethe Barkas Bloch equation for stopping power. Our calculations don't use Katz, but various similar calculations do.
  ion = c(rep("He4", 8), rep("Fe56", 17), rep("p", 4), rep("Ne20", 5), rep("Nb93", 4), rep("Fe56", 3), rep("Si28", 4), rep("Ti48", 5), rep("La139", 3)),
  comments = c(".93AlpLooksOK", rep("", 7), ".93AlplooksOK", rep("", 11), ".93Alp.no.iso", "not in 17Cuc (or 16Chang?)", rep("", 3), "16Chang all OK?", rep('', 24), ".94Alp","From graphs",'e.g. in 17Cuc')
) 

# Data for HG induced by photons from Cs-137 or Co-60 beta decay; from 16Chang (and calibration of LQ model)
ddd <- data.frame(dose.1 = c(0, 0.4, 0.8, 1.6, 3.2, 7, 0, .4, .8, .12, 1.6),
                  HG = c(.026, .048, .093, .137, .322, .462, .0497, .054, .067, .128, .202),
                  NWeight = c(6081.2, 4989.5, 1896.8, 981.1, 522.2, 205.2, 7474.1, 2877.6, 1423.7, 689.9, 514.9),
                  Nucleus = c(rep("Cobalt-60", 6), rep("Cesium-137", 5)),
                  Comments = c(rep("TBD", 11))
)
GeVu <- 0.001 * dfr[, "MeVperu"] #  convert to GeV/u for convenience in a calculation
dfr[, "Katz"] <- round(dfr[, "Z"] ^2 * (2.57 * GeVu ^2 + 4.781 * GeVu + 2.233) / (2.57 * GeVu ^2 + 4.781 * GeVu), 2) #  special relativistic calculation of Z^2/beta^2. The numerics include conversion from GeV to joules and from u to kg.
dfr[, "beta"] <- round(dfr[, "Z"] * sqrt(1 / dfr[, "Katz"]), 3) #  i.e. Z*sqrt(beta^2/Z^2) 
dfr[, "Zeff"] <- round(dfr[, "Z"] * (1 - exp( -125 * dfr[, "Z"] ^ (-2.0 / 3))), 2) #  Barkas formula for Zeff; for us Zeff is almost Z
dfr <- within(dfr, L[L < 200 & ion == 'Fe56'] <- 195) # Set all Fe56 with L < 200 to L = 185 
dfra <- dfr[c(1:19, 26:53), ] #  removes the zero dose case and the no isograft data
#=========================== DATA END ===========================#

#####  photon model #####
LQ <- lm(HG ~ dose.1 + I(dose.1 ^ 2), data = ddd) # linear model fit on ddd dataset
summary(LQ, correlation = T) 

#===================== HZE/NTE MODEL, abbreviated  "hin" for "high non-targeted" =====================# 
# Uses 3 adjustable parameters. There is also an HZE/TE model, abbreviated "hit" for "high targeted" and a "LOW"
# model for Z <= 3. Both hin and hit are for Z>3 in principle and here have data for Z >= 8.  

dfrHZE <- subset(dfra, Z > 3) # look only at HZE not at much lower Z and LET ions. # In next line phi controls how fast NTE build up from zero; not really needed during calibration since phi*Dose>>1 at every observed Dose !=0. phi needed for later synergy calculations.

phi <- 2000#  even larger phi should give the same final results, but might cause extra problems with R. 
hinm <- nls(HG ~ .0275 + (1 - exp ( -0.01 * (aa1 * L * dose.1 * exp( -aa2 * L) + (1 - exp( - phi * dose.1)) * kk1))), #  calibrating parameters in a model that modifies the hazard function NTE models in 17Cuc. "hinm" is for hin model
            data = dfrHZE, 
            weights = NWeight,
            start = list(aa1 = .9, aa2 = .01, kk1 = 6)) 
summary(hinm, correlation = T); vcov(hinm) #  parameter values & accuracy; variance-covariance matrix RKSB
hin.c <- coef(hinm) #  calibrated central values of the 3 parameters. Next is the IDER, = 0 at dose 0
hanC <- function(dose.1,L) { #  calibrated hazard function "hanC" is for hazard non-targeted calibrated
  0.01 * (hin.c[1] * L * dose.1 * exp(-hin.c[2] * L) + (1 - exp(- phi * dose.1)) * hin.c[3])
} 
Calculate.hinC <- function(dose.1, L) {
  1 - exp(-hanC(dose.1, L)) #  Calibrated HZE NTE IDER
}
######### TE model #########
hitm <- nls(HG ~ .0275 + (1 - exp ( -0.01 * (aate1 * L * dose.1 * exp( -aate2 * L) ))), #  calibrating parameters in a TE only model.
            data = dfrHZE,  
            weights = NWeight,
            start = list(aate1 = .9, aate2 = .01)) 
summary(hitm, correlation = T); vcov(hitm) #  parameter values & accuracy; variance-covariance matrix RKSB
hit.c <- coef(hitm) #  calibrated central values of the 2 parameters. Next is the IDER, = 0 at dose 0
hatC <- function(dose.1,L) { #  calibrated hazard function
  0.01 * (hit.c[1] * L * dose.1 * exp(-hit.c[2] * L))
} 
Calculate.hitC <- function(dose.1, L) {
  1 - exp(-hatC(dose.1, L)) #  Calibrated HZE TE IDER
}
IC<-cbind(AIC(hitm,hinm),BIC(hitm,hinm))
print(IC)
dose <- c(seq(0, .00001, by = 0.000001), #  look carefully near zero, but go out to 0.5 Gy
          seq(.00002, .0001, by=.00001),
          seq(.0002, .001, by=.0001),
          seq(.002, .01, by=.001),
          seq(.02, .5, by=.01))
# dose <- dose[1:30] #  this can be used to zoom in on the very low dose behavior in the graphs

####### calculate baseline MIXDER I(d) for mixtures of HZE components modeled by NTE IDER and then those by TE IDER #######
IntegratehinMIXDER <- function(r, L, d = dose, aa1 = hin.c[1], aa2 = hin.c[2], kk1 = hin.c[3]) {
  dE <- function(yini, State, Pars) {
    aa1 <- aa1; aa2 <- aa2; kk1 <- kk1
    with(as.list(c(State, Pars)), {
      aa = vector(length = length(L))
      u = vector(length = length(L))
      for (i in 1:length(L)) {
        aa[i] = aa1*L[i]*exp(-aa2*L[i])
        u[i] = uniroot(function(d) 1-exp(-0.01*(aa[i]*d+(1-exp(-phi*d))*kk1)) - I, lower = 0, upper = 20, tol = 10^-10)$root
      }
      dI = vector(length = length(L))
      for (i in 1:length(L)) {
        dI[i] = r[i]*0.01*(aa[i]+exp(-phi*u[i])*kk1*phi)*exp(-0.01*(aa[i]*u[i]+(1-exp(-phi*u[i]))*kk1))
        
      }
      dI = sum(dI)
      return(list(c(dI)))
    })
  }
  pars = NULL; yini = c(I= 0); d = d
  out = ode(yini, times = d, dE, pars, method = "radau")
  return(out)
} 
#Now hit instead of hin
Integrate_hiteMIXDER <- function(r, L, d = dose, aate1 = hit.c[1], aate2 = hit.c[2]) {
  dE <- function(yini, State, Pars) {
    aate1 <- aate1; aate2 <- aate2; kk1 <- kk1
    with(as.list(c(State, Pars)), {
      aate = vector(length = length(L))
      u = vector(length = length(L))
      for (i in 1:length(L)) {
        aate[i] = aate1*L[i]*exp(-aate2*L[i])
        u[i] = uniroot(function(d) 1-exp(-0.01*(aate[i]*d)) - I, lower = 0, upper = 20, tol = 10^-10)$root
      }
      dI = vector(length = length(L))
      for (i in 1:length(L)) {
        dI[i] = r[i]*0.01*aate[i]*exp(-0.01*(aate[i]*u[i]))
        
      }
      dI = sum(dI)
      return(list(c(dI)))
    })
  }
  pars = NULL; yini = c(I= 0); d = d
  out = ode(yini, times = d, dE, pars, method = "radau")
  return(out)
} 
########### Light ion, low Z (<= 3), low LET model ######### 
dfrL <- subset(dfra, Z <= 3) #  for Light ions
LOW.m <- nls(HG ~ .0275 + 1-exp(-bet * dose.1),
             data = dfrL,
             weights = NWeight,
             start = list(bet = .5))
summary(LOW.m)
LOW.c <- coef(LOW.m)  # calibrated central values of the parameter
CalculateLOW.C <- function(dose.1, L) { # Calibrated Low LET model. Use L=0, but maybe later will use L >0 but small 
  return(1 - exp(-LOW.c[1] * dose.1))
}  

dE_2 <- function(dose,L) { # Slope dE/dd of the low LET, low Z model; looking at the next plot() it seems fine
  LOW.c*exp(-LOW.c*dose)  
}

# plot () chunks such as the following are visual check to see if our calibration is consistent with 16Chang, .93Alp, .94Alp
# and 17Cuc; (ggplot commands are Yinmin's and concern CI)
# Put various values in our calibrated model to check with numbers and graphs in these references
# plot(c(0, 7), c(0, 1), col = 'red', ann = 'F') 
# ddose <- 0.01 * 0:700; lines(ddose, CalculateLOW.C(ddose, 0) + .0275)  #  calibrated lowLET IDER
# points(dfrL[1:8, "dose.1"], dfrL[1:8,"HG"],pch=19) #  RKS: Helium data points
# points(dfrL[9:12, "dose.1"], dfrL[9:12, "HG"] )  #  proton data points 

################## I(d) calculator for high Z, high E, NTE model hinm plus optionally LOW. START ##################
calculateComplexId <- function(r, L, d, aa1 = hin.c[1], aa2 = hin.c[2], kk1 = hin.c[3], phi = 2000, beta = LOW.c, lowLET = FALSE) {
  # Calculates incremental effect additivity function I(d) for mixture of N >= 1 HZE NTE IDERs and optionally one low-LET IDER 
  # new argument: lowLET (FALSE by default, TRUE when one IDER is low-LET)
  dE <- function(yini, State, Pars) { #  Constructing an ode from the IDERS
    aa1 <- aa1; aa2 <- aa2; kk1 <- kk1; beta <- beta; phi <- phi; L <- L
    with(as.list(c(State, Pars)), {
      aa <- vector(length = length(L))  
      u <- vector(length = length(L))  
      for (i in 1:length(L)) {
        aa[i] <- aa1 * L[i] * exp(-aa2 * L[i])
        u[i] <- uniroot(function(d) 1-exp(-0.01*(aa1*L[i]*d*exp(-aa2*L[i])+(1-exp(-phi*d))*kk1)) - I, lower = 0, upper = 200, extendInt = "yes", tol = 10^-10)$root #egh this is used in the single HZE and lowLET example
      }
      dI <- vector(length = length(L))
      for (i in 1:length(L)) {
        dI[i] <- r[i] * 0.01*(aa[i]+exp(-phi*u[i])*kk1*phi)*exp(-0.01*(aa[i]*u[i]+(1-exp(-phi*u[i]))*kk1))
      }
      if (lowLET == TRUE) { # If low-LET IDER is present then include it at the end of the dI vector
        u[length(L) + 1] <- uniroot(function(d) 1-exp(-beta*d) - I, lower = 0, upper = 200, extendInt = "yes", tol = 10^-10)$root
        dI[length(L) + 1] <- r[length(r)] * dE_2(d = u[length(L) + 1], L = 0)
      }
      dI <- sum(dI)
      return(list(c(dI)))
    })
  }
  return(ode(c(I = 0), times = d, dE, parms = NULL, method = "radau")) #  Finds solution I(d) of the differential equation
  # RKS to Yimin and Edward: I'm not convinced we need to or should add that the method is radau
}
################## I(d) calculator for high Z, high E, NTE model hinm plus optionally LOW. END ##################

###### RKS to Yimin and Edward: Next is the same for hitm; then some plots #####
calculateComplexId.te <- function(r, L, d, aate1 = hit.c[1], aate2 = hit.c[2], beta = LOW.c, lowLET = FALSE) {
  dE <- function(yini, State, Pars) { #  Constructing an ode from the IDERS
    aate1 <- aate1; aate2 <- aate2; beta <- beta; L <- L
    with(as.list(c(State, Pars)), {
      aate <- vector(length = length(L))
      u <- vector(length = length(L))
      for (i in 1:length(L)) {
        aate[i] <- aate1 * L[i] * exp(-aate2 * L[i])
        u[i] <- uniroot(function(d) 1-exp(-0.01*(aate1*L[i]*d*exp(-aate2*L[i]))) - I, lower = 0, upper = 200, extendInt = "yes", tol = 10^-10)$root
      }
      dI <- vector(length = length(L))
      for (i in 1:length(L)) {
        dI[i] <- r[i] * 0.01*aate[i]*exp(-0.01*aate[i]*u[i])
      }
      if (lowLET == TRUE) { # If low-LET IDER is present then include it at the end of the dI vector
        u[length(L) + 1] <- uniroot(function(d) 1-exp(-beta*d) - I, lower = 0, upper = 200, extendInt = "yes", tol = 10^-10)$root
        dI[length(L) + 1] <- r[length(r)] * dE_2(d = u[length(L) + 1], L = 0)
      }
      dI <- sum(dI)
      return(list(c(dI)))
    })
  }
  return(ode(c(I = 0), times = d, dE, parms = NULL, method = "radau")) #  Finds solution I(d) of the differential equation
}
# RKS to Yimin and Edward: I'm not convinced we need to or should add that the method is radau

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

#==============================================#
#=================  Plotting  =================#
#==============================================#

library(grid)

# https://stackoverflow.com/questions/27667017/removing-right-border-from-ggplot2-graph
element_grob.element_custom <- function(element, ...)  {
  
  segmentsGrob(c(1),
               c(1),
               c(1),
               c(0), gp=gpar(lwd=3))
}
border_custom <- function(...){
  structure(
    list(...), 
    class = c("element_custom","element_blank", "element")
  ) 
  
}

# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplotHelper <- function(plotArr, file_path, num_col, w, h) {
  if (length(plotArr) == 8) {
    multiplot(plotArr[[1]], plotArr[[2]], plotArr[[3]], plotArr[[4]], plotArr[[5]], plotArr[[6]], plotArr[[7]], plotArr[[8]], cols = num_col)
  } else if (length(plotArr) == 2) {
    multiplot(plotArr[[1]], plotArr[[2]], cols = num_col)
  }
  dev.copy2eps(file = file_path, width = w, height = h)
}
  

plottingHelper <- function(Lval, blackwhite, zoom = 0, save_path = "~/Desktop/plots/", has_axis_text = TRUE, is_multiplot = FALSE, output_to_file = TRUE) {
  # zoom = 0       no zoom
  #        1       zoom for Fe 0.4
  #        2       zoom for Fe 0.2
  
  #######################
  # Data initialization #
  #######################
  d = c()
  nWeight = c()
  hgArr = c()
  
  dProton = c()
  nWeightProton = c()
  protonArr = c()
  
  dfCurve = NULL
  dfData = NULL
  dfProton = NULL
  
  if (blackwhite) {
    color1 = "black"
    color2 = "black"
  } else {
    color1 = "red"
    color2 = "blue"
  }
  
  
  ############################
  # Read data from dataframe #
  ############################
  for (row in 1:nrow(dfr)) {
    if (Lval == 1) {                                            # LOW LET
      if (dfr[row, "L"] == 0.4) {
        dProton = c(dProton, dfr[row, "dose.1"])
        nWeightProton = c(nWeightProton, dfr[row, "NWeight"])
        protonArr = c(protonArr, dfr[row, "HG"])
      } else if (dfr[row, "L"] == 1.6) {
        d = c(d, dfr[row, "dose.1"])
        nWeight = c(nWeight, dfr[row, "NWeight"])
        hgArr = c(hgArr, dfr[row, "HG"])
      }
    } else if (Lval == 193 || Lval == 195) {                    # 193, 195
      if (!zoom) {
        if (dfr[row, "L"] == 193 || dfr[row, "L"] == 195) {
          d = c(d, dfr[row, "dose.1"])
          nWeight = c(nWeight, dfr[row, "NWeight"])
          hgArr = c(hgArr, dfr[row, "HG"])
        }
      } else {
        if (((zoom == 1) && (dfr[row, "L"] == 193 || dfr[row, "L"] == 195) && (dfr[row, "dose.1"] <= 0.4)) ||
            ((zoom == 2) && (dfr[row, "L"] == 193 || dfr[row, "L"] == 195) && (dfr[row, "dose.1"] <= 0.2))) {
          d = c(d, dfr[row, "dose.1"])
          nWeight = c(nWeight, dfr[row, "NWeight"])
          hgArr = c(hgArr, dfr[row, "HG"])
        }
      }
    } else if (dfr[row, "L"] == Lval) {                         # Other LET
      d = c(d, dfr[row, "dose.1"])
      nWeight = c(nWeight, dfr[row, "NWeight"])
      hgArr = c(hgArr, dfr[row, "HG"])
    }
  }
  
  if (Lval == 1) {
    endpoint = 7.2
  } else if (Lval == 195) {
    if (!zoom) {
      endpoint = 1.6
    } else if (zoom == 1) {
      endpoint = 0.4
    } else if (zoom == 2) {
      endpoint = 0.2
    }
  } else {
    endpoint = d[length(d)]
  }
  
  dose <- c(seq(0, .00001, by = 0.000001),
            seq(.00002, .0001, by=.00001),
            seq(.0002, .001, by=.0001),
            seq(.002, .01, by=.001),
            seq(.02, endpoint, by=.01))
  
  ###################
  # Actual plotting #
  ###################
  
  #################################
  # Plot parameter initialization #
  #################################
  
  p = ggplot() + theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme(axis.ticks.length=unit(0.4,"cm"))
  
  if (has_axis_text) {
    p = p + theme(axis.ticks.x = element_line(size = 1)) + 
      theme(axis.ticks.y = element_line(size = 1))
  } else {
    p = p + theme(axis.ticks.x = element_line(size = 1), axis.text.x=element_blank()) +
      theme(axis.ticks.y = element_line(size = 1), axis.text.y=element_blank())
  }
  
  if (is_multiplot) {
    if (Lval != 100 && Lval != 953) {
      if (!zoom) {
        p = last_plot() + theme(panel.border = border_custom())
      } else if (zoom == 1) {
        p = last_plot() + theme(panel.border = border_custom())
      }
    }
    
    if (!zoom) {
      p = p + theme(plot.margin = unit(c(1,1.5,1,1.5), "cm"))
    } else {
      p = p + theme(plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))
    }
  }
  
  #################
  # Data plotting #
  #################
  
  if (Lval == 1) {
    lowVal = CalculateLOW.C(dose, Lval)
    dfCurve = data.frame(d = dose, low = lowVal)
    dfData = data.frame(x = d, hg = hgArr, nw = nWeight)
    dfProton = data.frame(x = dProton, hg = protonArr, nw = nWeightProton)
  } else {
    if (Lval == 193 || Lval == 195) {
      hinVal = Calculate.hinC(dose, 180)
      hitVal = Calculate.hitC(dose, 180)
    } else {
      hinVal = Calculate.hinC(dose, Lval)
      hitVal = Calculate.hitC(dose, Lval)
    }
    dfCurve = data.frame(d = dose, hgNTE = hinVal, hgTE = hitVal)
    dfData = data.frame(x = d, hg = hgArr, nw = nWeight)
  }
  
  
  if (Lval == 1) {
    p = p + 
        geom_line(data = dfCurve, aes(x = dose, y = low), colour = color1, size = 1) +
        geom_segment(data = dfData, size = 0.8, aes(x = d, xend = d, y = hg - 1/sqrt(nWeight), yend = hg + 1/sqrt(nWeight)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) + 
        geom_segment(data = dfData, size = 0.8, aes(x = d, xend = d, y = hg, yend = hg - 1/sqrt(nWeight)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) +
        geom_segment(data = dfProton, size = 0.8, aes(x = dProton, xend = dProton, y = hg - 1/sqrt(nWeightProton), yend = hg + 1/sqrt(nWeightProton)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) + 
        geom_segment(data = dfProton, size = 0.8, aes(x = dProton, xend = dProton, y = hg, yend = hg - 1/sqrt(nWeightProton)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) +
        geom_point(data = dfData, aes(x = d, y = hg), size = 4, fill ="white", shape = 21) + 
        geom_point(data = dfProton, aes(x = dProton, y = hg), size = 4, fill ="black", shape = 21)
  } else {
    p = p + 
        geom_line(data = dfCurve, aes(x = dose, y = hgNTE), colour = color1, size = 1) +
        geom_line(data = dfCurve, aes(x = dose, y = hgTE), colour = color2, size = 1, linetype = "dashed") + 
        geom_segment(data = dfData, size = 0.8, aes(x = d, xend = d, y = hg - 1/sqrt(nWeight), yend = hg + 1/sqrt(nWeight)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) + 
        geom_segment(data = dfData, size = 0.8, aes(x = d, xend = d, y = hg, yend = hg - 1/sqrt(nWeight)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) +
        geom_point(data = dfData, aes(x = d, y = hg), size = 2, fill ="black", shape = 21)
  }
  
  ###############
  # Plot output #
  ###############
  if (output_to_file) {
    name = " "
    if (blackwhite) {
      if (zoom == 1) {
        name = paste(save_path, "L =", toString(Lval), "zoomed max dose 0.4", "Blackwhite.eps")
      } else if (zoom == 2) {
        name = paste(save_path, "L =", toString(Lval), "zoomed max dose 0.2", "Blackwhite.eps")
      } else {
        name = paste(save_path, "L =", toString(Lval), "Blackwhite.eps")
      }
      ggsave(filename = name, plot = p, device = "eps")
    } else {
      name = " "
        if (zoom == 1) {
          name = paste(save_path, "L =", toString(Lval), "zoomed max dose 0.4", "Color.eps")
        } else if (zoom == 2) {
          name = paste(save_path, "L =", toString(Lval), "zoomed max dose 0.2", "Color.eps")
        } else {
          name = paste(save_path, "L =", toString(Lval), "Color.eps")
        }
      ggsave(filename = name, plot = p, device = "eps")
    }
  }
  return (p)
}




##########################
# Plot IDER for each LET #
##########################

possibleLVal = c(1, 25, 70, 100, 195, 250, 464, 953)

for (val in possibleLVal) {
  plottingHelper(val, TRUE)
  plottingHelper(val, FALSE)
}

#################################
# Different kinds of multiplots #
#################################

multiplotArrBW = list(0)
multiplotArrColor = list(0)
twoPanelArrBW = list(0)
twoPanelArrColor = list(0)
count = 1

for (val in possibleLVal) {
  p1 = plottingHelper(val, TRUE, is_multiplot = TRUE, has_axis_text = FALSE, output_to_file = FALSE)
  p2 = plottingHelper(val, FALSE, is_multiplot = TRUE, has_axis_text = FALSE, output_to_file = FALSE)
  multiplotArrBW[[count]] = p1
  multiplotArrColor[[count]] = p2
  count = count + 1
}

# permute multiplot array to change the order of plots inside the multiplot
multiplotArrColor = multiplotArrColor[c(1,5,2,6,3,7,4,8)]
multiplotArrBW = multiplotArrBW[c(1,5,2,6,3,7,4,8)]

twoPanelArrBW[[1]] = multiplotArrBW[[1]]
twoPanelArrBW[[2]] = multiplotArrBW[[2]]
twoPanelArrColor[[1]] = multiplotArrColor[[1]]
twoPanelArrColor[[2]] = multiplotArrColor[[2]]


# blackwhite 2 panel
multiplotHelper(twoPanelArrBW, "~/Desktop/plots/1*2 Multiplot Blackwhite, LOW LET + 195.eps", num_col = 2, w = 8, h = 4)

# color 2 panel
multiplotHelper(twoPanelArrColor, "~/Desktop/plots/1*2 Multiplot Color, LOW LET + 195.eps", num_col = 2, w = 8, h = 4)

# blackwhite 8 panel
multiplotHelper(multiplotArrBW, "~/Desktop/plots/4*2 Multiplot Blackwhite.eps", num_col = 2, w = 10, h = 20)
multiplotHelper(multiplotArrBW, "~/Desktop/plots/2*4 Multiplot Blackwhite.eps", num_col = 4, w = 20, h = 10)

# Color 8 panel
multiplotHelper(multiplotArrColor, "~/Desktop/plots/4*2 Multiplot Color.eps", num_col = 2, w = 10, h = 20)
multiplotHelper(multiplotArrColor, "~/Desktop/plots/2*4 Multiplot Color.eps", num_col = 4, w = 20, h = 10)

##############################
# Different kinds IDER plots #
##############################

# Some parameter initializations
save_path =  "~/Desktop/plots/"
dose <- c(seq(0, .00001, by = 0.000001), 
          seq(.00002, .0001, by=.00001),
          seq(.0002, .001, by=.0001),
          seq(.002, .01, by=.001),
          seq(.02, 1., by=.01))

# 80% 195, 20% low
r <- c(0.8, 0.2); L <- c(195)
incremental = calculateComplexId(r, L, d = dose, lowLET = TRUE)
highVal = Calculate.hinC(dose, 195)
lowVal = CalculateLOW.C(dose, 1)

df1 = data.frame(d = dose, hg = incremental[, 2])
df2 = data.frame(d = dose, hg = highVal)
df3 = data.frame(d = dose, hg = lowVal)
p = ggplot() + theme_bw() + 
  geom_line(data = df2, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df3, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 2) + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm"))

print(p)

ggsave(filename = paste(save_path, "IDER 80percent_195, 20percent_LOW", ".eps"), plot = p, device = "eps")

# 20% 195, 80% low
r <- c(0.2, 0.8); L <- c(195)
incremental = calculateComplexId(r, L, d = dose, lowLET = TRUE)
highVal = Calculate.hinC(dose, 195)
lowVal = CalculateLOW.C(dose, 1)

df1 = data.frame(d = dose, hg = incremental[, 2])
df2 = data.frame(d = dose, hg = highVal)
df3 = data.frame(d = dose, hg = lowVal)
p = ggplot() + theme_bw() + 
  geom_line(data = df2, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df3, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 2) + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm"))

ggsave(filename = paste(save_path, "IDER 20percent_195, 80percent_LOW", ".eps"), plot = p, device = "eps")

# 80% low, 5% 25,  5% 70,  5% 100,  5% 195
r <- c(0.05, 0.05, 0.05, 0.05, 0.8); L <- c(25, 70, 100, 195)
incremental = calculateComplexId(r, L, d = dose, lowLET = TRUE)
Val195 = Calculate.hinC(dose, 195)
Val25 = Calculate.hinC(dose, 25)
Val70 = Calculate.hinC(dose, 70)
Val100 = Calculate.hinC(dose, 100)
lowVal = CalculateLOW.C(dose, 1)

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df70 = data.frame(d = dose, hg = Val70)
df100 = data.frame(d = dose, hg = Val100)
df195 = data.frame(d = dose, hg = Val195)
dflow = data.frame(d = dose, hg = lowVal)
p = ggplot() + theme_bw() + 
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 2) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm"))

ggsave(filename = paste(save_path, "IDER 80percent_LOW, 5percent_25,70,100,195 each", ".eps"), plot = p, device = "eps")

# 20% low, 20% 25,  20% 70,  20% 100,  20% 195
r <- c(0.2, 0.2, 0.2, 0.2, 0.2); L <- c(25, 70, 100, 195)
incremental = calculateComplexId(r, L, d = dose, lowLET = TRUE)
Val195 = Calculate.hinC(dose, 195)
Val25 = Calculate.hinC(dose, 25)
Val70 = Calculate.hinC(dose, 70)
Val100 = Calculate.hinC(dose, 100)
lowVal = CalculateLOW.C(dose, 1)

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df70 = data.frame(d = dose, hg = Val70)
df100 = data.frame(d = dose, hg = Val100)
df195 = data.frame(d = dose, hg = Val195)
dflow = data.frame(d = dose, hg = lowVal)

p = ggplot() + theme_bw() + 
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 2) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm"))

ggsave(filename = paste(save_path, "IDER 20percent_LOW,25,70,100,195 each", ".eps"), plot = p, device = "eps")

# 0.1 each
r <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1); L <- c(25,50,75,100,150,200,300,500,750,1000)
incremental = calculateComplexId(r, L, d = dose, lowLET = FALSE)
SEA <- function(dose.1) Calculate.hinC(dose.1/10, 25) + Calculate.hinC(dose.1/10, 50) + Calculate.hinC(dose.1/10, 75) + Calculate.hinC(dose.1/10, 100) + CalculateLOW.C(dose.1/10, 150) + Calculate.hinC(dose.1/10, 200) + CalculateLOW.C(dose.1/10, 300) + Calculate.hinC(dose.1/10, 500) + CalculateLOW.C(dose.1/10, 750)+ CalculateLOW.C(dose.1/10, 1000)
Val25 = Calculate.hinC(dose, 25)
Val50 = Calculate.hinC(dose, 50)
Val75 = Calculate.hinC(dose, 75)
Val100 = Calculate.hinC(dose, 100)
Val150 = Calculate.hinC(dose, 150)
Val200 = Calculate.hinC(dose, 200)
Val300 = Calculate.hinC(dose, 300)
Val500 = Calculate.hinC(dose, 500)
Val750 = Calculate.hinC(dose, 750)
Val1000 = Calculate.hinC(dose, 1000)


dfsea = data.frame(d = dose, hg = SEA(dose))

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df50 = data.frame(d = dose, hg = Val50)
df75 = data.frame(d = dose, hg = Val75)
df100 = data.frame(d = dose, hg = Val100)
df150 = data.frame(d = dose, hg = Val150)
df200 = data.frame(d = dose, hg = Val200)
df300 = data.frame(d = dose, hg = Val300)
df500 = data.frame(d = dose, hg = Val500)
df750 = data.frame(d = dose, hg = Val750)
df1000 = data.frame(d = dose, hg = Val1000)

p = ggplot() + theme_bw() + 
  geom_line(data = dfsea, aes(x = d, y = hg), colour = "black", size = 1, linetype = "dashed") +
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df50, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df75, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df150, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df200, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df300, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df500, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df750, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df1000, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 1) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave(filename = paste(save_path, "IDER 10percent_LOW,25,50,75,100,150,200,300,500,750,1000 each", ".eps"), plot = p, device = "eps")

#############################################
# Different kinds confidence interval plots #
#############################################

# Some parameter initializations
sampleN = 200
mod = 1              # 1 if HIN, 0 if HIT
dose <- c(seq(0, .00001, by = 0.000001), 
          seq(.00002, .0001, by=.00001),
          seq(.0002, .001, by=.0001),
          seq(.002, .01, by=.001),
          seq(.02, 1., by=.01))

dose <- c(seq(0.0, 1., by = 0.1))

# 5x0.2 Confidence Interval Monte
confidenceIntervals = CIHelper(sampleNum = sampleN, d = dose, r = c(0.2, 0.2, 0.2, 0.2, 0.2), L = c(25, 70, 100, 180), model = 0)
naiveCI = confidenceIntervals[[1]]
monteCarloCI = confidenceIntervals[[2]]
incremental = calculateComplexId(r = c(0.2, 0.2, 0.2, 0.2, 0.2), L = c(25, 70, 100, 180), d = dose, lowLET = TRUE)
SEA <- function(dose.1) Calculate.hinC(dose.1/5, 25) + Calculate.hinC(dose.1/5, 70) + Calculate.hinC(dose.1/5, 100) + Calculate.hinC(dose.1/5, 180) + CalculateLOW.C(dose.1/5, 1)
Val195 = Calculate.hinC(dose, 180)
Val25 = Calculate.hinC(dose, 25)
Val70 = Calculate.hinC(dose, 70)
Val100 = Calculate.hinC(dose, 100)
lowVal = CalculateLOW.C(dose, 1)

dfsea = data.frame(d = dose, hg = SEA(dose))

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df70 = data.frame(d = dose, hg = Val70)
df100 = data.frame(d = dose, hg = Val100)
df195 = data.frame(d = dose, hg = Val195)
dflow = data.frame(d = dose, hg = lowVal)
ciPlot1 = ggplot() + theme_bw() + 
  geom_ribbon(aes(x = dose, ymin = monteCarloCI[1, ], ymax = monteCarloCI[2, ]), fill = "yellow") +
  geom_line(data = dfsea, aes(x = d, y = hg), colour = "black", size = 1, linetype = "dashed") +
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 1) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

print(ciPlot1)

ggsave(filename = paste(save_path, "ConfidenceIntervalMonteCarlo 20percent_LOW,25,70,100,195 each", ".eps"), plot = ciPlot1, device = "eps")

# 5x0.2 Confidence Interval Naive
incremental = calculateComplexId(r = c(0.2, 0.2, 0.2, 0.2, 0.2), L = c(25, 70, 100, 180), d = dose, lowLET = TRUE)
Val195 = Calculate.hinC(dose, 180)
Val25 = Calculate.hinC(dose, 25)
Val70 = Calculate.hinC(dose, 70)
Val100 = Calculate.hinC(dose, 100)
lowVal = CalculateLOW.C(dose, 1)

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df70 = data.frame(d = dose, hg = Val70)
df100 = data.frame(d = dose, hg = Val100)
df195 = data.frame(d = dose, hg = Val195)
dflow = data.frame(d = dose, hg = lowVal)
ciPlot2 = ggplot() + theme_bw() + 
  geom_ribbon(aes(x = dose, ymin = naiveCI[1, ], ymax = naiveCI[2, ]), fill = "yellow") +
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 1) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave(filename = paste(save_path, "ConfidenceIntervalNaive 20percent_LOW,25,70,100,195 each", ".eps"), plot = ciPlot2, device = "eps")

twoPanelArr = list(0)
twoPanelArr[[1]] = ciPlot1
twoPanelArr[[2]] = ciPlot2
multiplotHelper(twoPanelArr, "~/Desktop/plots/1*2 ConfidenceInterval 20percent_LOW,25,70,100,195 each.eps", num_col = 2, w = 8, h = 4)


# 0.05 0.8 Confidence Interval Monte
r <- c(0.05, 0.05, 0.05, 0.05, 0.8);
confidenceIntervals = CIHelper(sampleNum = sampleN, d = dose, r = c(0.05, 0.05, 0.05, 0.05, 0.8), L = c(25, 70, 100, 180), model = 0)
naiveCI = confidenceIntervals[[1]]
monteCarloCI = confidenceIntervals[[2]]
incremental = calculateComplexId(r = c(0.05, 0.05, 0.05, 0.05, 0.8), L = c(25, 70, 100, 180), d = dose, lowLET = TRUE)
SEA <- function(dose.1) Calculate.hinC(dose.1/20, 25) + Calculate.hinC(dose.1/20, 70) + Calculate.hinC(dose.1/20, 100) + Calculate.hinC(dose.1/20, 180) + CalculateLOW.C(4*dose.1/5, 1)
Val195 = Calculate.hinC(dose, 180)
Val25 = Calculate.hinC(dose, 25)
Val70 = Calculate.hinC(dose, 70)
Val100 = Calculate.hinC(dose, 100)
lowVal = CalculateLOW.C(dose, 1)

dfsea = data.frame(d = dose, hg = SEA(dose))

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df70 = data.frame(d = dose, hg = Val70)
df100 = data.frame(d = dose, hg = Val100)
df195 = data.frame(d = dose, hg = Val195)
dflow = data.frame(d = dose, hg = lowVal)
ciPlot1 = ggplot() + theme_bw() + 
  geom_ribbon(aes(x = dose, ymin = monteCarloCI[1, ], ymax = monteCarloCI[2, ]), fill = "yellow") +
  geom_line(data = dfsea, aes(x = d, y = hg), colour = "black", size = 1, linetype = "dashed") +
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 1) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave(filename = paste(save_path, "ConfidenceIntervalMonteCarlo 80percent_LOW, 5percent_25,70,100,195 each", ".eps"), plot = ciPlot1, device = "eps")

# 0.05 0.8 Confidence Interval Naive
incremental = calculateComplexId(r = c(0.05, 0.05, 0.05, 0.05, 0.8), L = c(25, 70, 100, 180), d = dose, lowLET = TRUE)
Val195 = Calculate.hinC(dose, 180)
Val25 = Calculate.hinC(dose, 25)
Val70 = Calculate.hinC(dose, 70)
Val100 = Calculate.hinC(dose, 100)
lowVal = CalculateLOW.C(dose, 1)

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df70 = data.frame(d = dose, hg = Val70)
df100 = data.frame(d = dose, hg = Val100)
df195 = data.frame(d = dose, hg = Val195)
dflow = data.frame(d = dose, hg = lowVal)
ciPlot2 = ggplot() + theme_bw() + 
  geom_ribbon(aes(x = dose, ymin = naiveCI[1, ], ymax = naiveCI[2, ]), fill = "yellow") +
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 1) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave(filename = paste(save_path, "ConfidenceIntervalNaive 80percent_LOW, 5percent_25,70,100,195 each", ".eps"), plot = ciPlot2, device = "eps")

twoPanelArr[[1]] = ciPlot1
twoPanelArr[[2]] = ciPlot2
multiplotHelper(twoPanelArr, "~/Desktop/plots/1*2 ConfidenceInterval 80percent_LOW, 5percent_25,70,100,195 each.eps", num_col = 2, w = 8, h = 4)

# Zoom in Fe
p1 = plottingHelper(195, TRUE, save_path = save_path, zoom = 1)
p2 = plottingHelper(195, TRUE, save_path = save_path, zoom = 2)
twoPanelArr[[1]] = p1
twoPanelArr[[2]] = p2
multiplotHelper(twoPanelArr, "~/Desktop/plots/1*2 Fe Zoomed Blackwhite.eps", num_col = 2, w = 8, h = 4)

p1 = plottingHelper(195, FALSE, save_path = save_path, zoom = 1)
p2 = plottingHelper(195, FALSE, save_path = save_path, zoom = 2)
twoPanelArr[[1]] = p1
twoPanelArr[[2]] = p2
multiplotHelper(twoPanelArr, "~/Desktop/plots/1*2 Fe Zoomed Color.eps", num_col = 2, w = 8, h = 4)







###############
# Feb 9 plots #
###############

dose1 = 0.01 * 0:100                 # dose range
dose2 = 0.0001 * 0:100
dose3 = 0.00001 * 0:100

df1 = data.frame(d = dose1*100, hg = Calculate.hinC(dose.1 = dose1, L = 185))
df2 = data.frame(d = dose2*100, hg = Calculate.hinC(dose.1 = dose2, L = 185))
df3 = data.frame(d = dose3*100, hg = Calculate.hinC(dose.1 = dose3, L = 185))
p1 = ggplot() + theme_bw() + 
  geom_line(data = df1, aes(x = d, y = hg), colour = "#6AB490", size = 1) + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(panel.border = border_custom()) +
  annotate("text", x = 50, y = 0.1, label = "A", fontface = 2, size = 6) +
  labs(x = "Dose d (cGy)", y = "HG prevalence") 
p2 = ggplot() + theme_bw() + 
  geom_line(data = df2, aes(x = d, y = hg), colour = "#6AB490", size = 1) + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(panel.border = border_custom()) +
  annotate("text", x = 0.5, y = 0.012, label = "B", fontface = 2, size = 6) +
  labs(x = "Dose d (cGy)")
p3 = ggplot() + theme_bw() + 
  geom_line(data = df3, aes(x = d, y = hg), colour = "#6AB490", size = 1) + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  annotate("text", x = 0.05, y = 0.009, label = "C", fontface = 2, size = 6) +
  labs(x = "Dose d (cGy)")

multiplot(p1, p2, p3, cols = 3)

d = c()
nWeight = c()
hgArr = c()

for (row in 1:nrow(dfr)) {
  if ((dfr[row, "L"] == 193 || dfr[row, "L"] == 195) && (dfr[row, "dose.1"] != 1.6)) {
    print(122)
    d = c(d, dfr[row, "dose.1"])
    nWeight = c(nWeight, dfr[row, "NWeight"])
    hgArr = c(hgArr, dfr[row, "HG"])
  }
}

dfData = data.frame(x = d*100, hg = hgArr, nw = nWeight)

p1 = p1 + 
  geom_segment(data = dfData, size = 0.8, aes(x = x, xend = x, y = hg - 1/sqrt(nw), yend = hg + 1/sqrt(nw)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) + 
  geom_segment(data = dfData, size = 0.8, aes(x = x, xend = x, y = hg, yend = hg - 1/sqrt(nw)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) +
  geom_point(data = dfData, aes(x = x, y = hg), size = 2, fill ="black", shape = 21)

multiplot(p1, p2, p3, cols = 3)


###################################################
########### Feb 22 Check lower dose CI ############
###################################################

# Some parameter initializations
sampleN = 200
mod = 1              # 1 if HIN, 0 if HIT

dose <- c(seq(0.0, 1., by = 0.1))

# 5x0.2 Confidence Interval Monte NTE
lowerDoseCI = CIHelper(sampleNum = sampleN, d = dose, r = c(0.2, 0.2, 0.2, 0.2, 0.2), L = c(25, 70, 100, 180), model = 0, seed = 100)

dose <- c(seq(0, .00001, by = 0.000001), 
          seq(.00002, .0001, by=.00001),
          seq(.0002, .001, by=.0001),
          seq(.002, .01, by=.001),
          seq(.02, 1., by=.01))

higherDoseCI = CIHelper(sampleNum = sampleN, d = dose, r = c(0.2, 0.2, 0.2, 0.2, 0.2), L = c(25, 70, 100, 180), model = 0, seed = 99)

err1NTE = abs(higherDoseCI[[2]][1, seq(length = 10, from = 47, by = 10)] - lowerDoseCI[[2]][1, seq(length = 10, from = 2, by = 1)])

# 0.05 0.8 Confidence Interval Monte NTE
dose <- c(seq(0.0, 1., by = 0.1))

lowerDoseCI = CIHelper(sampleNum = sampleN, d = dose, r = c(0.05, 0.05, 0.05, 0.05, 0.8), L = c(25, 70, 100, 180), model = 0, seed = 100)

dose <- c(seq(0, .00001, by = 0.000001), 
          seq(.00002, .0001, by=.00001),
          seq(.0002, .001, by=.0001),
          seq(.002, .01, by=.001),
          seq(.02, 1., by=.01))

higherDoseCI = CIHelper(sampleNum = sampleN, d = dose, r = c(0.05, 0.05, 0.05, 0.05, 0.8), L = c(25, 70, 100, 180), model = 0, seed = 99)

err2NTE = abs(higherDoseCI[[2]][1, seq(length = 10, from = 47, by = 10)] - lowerDoseCI[[2]][1, seq(length = 10, from = 2, by = 1)])

# Some parameter initializations
sampleN = 200
mod = 1              # 1 if HIN, 0 if HIT

dose <- c(seq(0.0, 1., by = 0.1))

# 5x0.2 Confidence Interval Monte TE
lowerDoseCI = CIHelper(sampleNum = sampleN, d = dose, r = c(0.2, 0.2, 0.2, 0.2, 0.2), L = c(25, 70, 100, 180), model = 1, seed = 100)

dose <- c(seq(0, .00001, by = 0.000001), 
          seq(.00002, .0001, by=.00001),
          seq(.0002, .001, by=.0001),
          seq(.002, .01, by=.001),
          seq(.02, 1., by=.01))

higherDoseCI = CIHelper(sampleNum = sampleN, d = dose, r = c(0.2, 0.2, 0.2, 0.2, 0.2), L = c(25, 70, 100, 180), model = 1, seed = 99)

err1TE = abs(higherDoseCI[[2]][1, seq(length = 10, from = 47, by = 10)] - lowerDoseCI[[2]][1, seq(length = 10, from = 2, by = 1)])

# 0.05 0.8 Confidence Interval Monte TE
dose <- c(seq(0.0, 1., by = 0.1))

lowerDoseCI = CIHelper(sampleNum = sampleN, d = dose, r = c(0.05, 0.05, 0.05, 0.05, 0.8), L = c(25, 70, 100, 180), model = 1, seed = 100)

dose <- c(seq(0, .00001, by = 0.000001), 
          seq(.00002, .0001, by=.00001),
          seq(.0002, .001, by=.0001),
          seq(.002, .01, by=.001),
          seq(.02, 1., by=.01))

higherDoseCI = CIHelper(sampleNum = sampleN, d = dose, r = c(0.05, 0.05, 0.05, 0.05, 0.8), L = c(25, 70, 100, 180), model = 1, seed = 99)

err2TE = abs(higherDoseCI[[2]][1, seq(length = 10, from = 47, by = 10)] - lowerDoseCI[[2]][1, seq(length = 10, from = 2, by = 1)])

print("5x0.2 Confidence Interval Monte NTE Error");
print(err1NTE);
print("0.05 0.8 Confidence Interval Monte NTE Error");
print(err2NTE);
print("5x0.2 Confidence Interval Monte TE Error");
print(err1TE);
print("0.05 0.8 Confidence Interval Monte TE Error");
print(err2TE)



########################################################################
#########################   TE VERSION   ###############################
########################################################################


##############################
# Different kinds IDER plots #
##############################

# Some parameter initializations
save_path =  "~/Desktop/plots/"
dose <- c(seq(0, .00001, by = 0.000001), 
          seq(.00002, .0001, by=.00001),
          seq(.0002, .001, by=.0001),
          seq(.002, .01, by=.001),
          seq(.02, 1., by=.01))

# 80% 195, 20% low
r <- c(0.8, 0.2); L <- c(195)
incremental = calculateComplexId.te(r, L, d = dose, lowLET = TRUE)
highVal = Calculate.hitC(dose, 195)
lowVal = CalculateLOW.C(dose, 1)

df1 = data.frame(d = dose, hg = incremental[, 2])
df2 = data.frame(d = dose, hg = highVal)
df3 = data.frame(d = dose, hg = lowVal)
p = ggplot() + theme_bw() + 
  geom_line(data = df2, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df3, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 2) + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm"))
print(p)


# 20% 195, 80% low
r <- c(0.2, 0.8); L <- c(195)
incremental = calculateComplexId.te(r, L, d = dose, lowLET = TRUE)
highVal = Calculate.hitC(dose, 195)
lowVal = CalculateLOW.C(dose, 1)

df1 = data.frame(d = dose, hg = incremental[, 2])
df2 = data.frame(d = dose, hg = highVal)
df3 = data.frame(d = dose, hg = lowVal)
p = ggplot() + theme_bw() + 
  geom_line(data = df2, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df3, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 2) + 
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm"))

ggsave(filename = paste(save_path, "IDER 20percent_195, 80percent_LOW", ".eps"), plot = p, device = "eps")

# 80% low, 5% 25,  5% 70,  5% 100,  5% 195
r <- c(0.05, 0.05, 0.05, 0.05, 0.8); L <- c(25, 70, 100, 195)
incremental = calculateComplexId.te(r, L, d = dose, lowLET = TRUE)
Val195 = Calculate.hitC(dose, 195)
Val25 = Calculate.hitC(dose, 25)
Val70 = Calculate.hitC(dose, 70)
Val100 = Calculate.hitC(dose, 100)
lowVal = CalculateLOW.C(dose, 1)

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df70 = data.frame(d = dose, hg = Val70)
df100 = data.frame(d = dose, hg = Val100)
df195 = data.frame(d = dose, hg = Val195)
dflow = data.frame(d = dose, hg = lowVal)
p = ggplot() + theme_bw() + 
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 2) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm"))

ggsave(filename = paste(save_path, "IDER 80percent_LOW, 5percent_25,70,100,195 each", ".eps"), plot = p, device = "eps")

# 20% low, 20% 25,  20% 70,  20% 100,  20% 195
r <- c(0.2, 0.2, 0.2, 0.2, 0.2); L <- c(25, 70, 100, 195)
incremental = calculateComplexId.te(r, L, d = dose, lowLET = TRUE)
Val195 = Calculate.hitC(dose, 195)
Val25 = Calculate.hitC(dose, 25)
Val70 = Calculate.hitC(dose, 70)
Val100 = Calculate.hitC(dose, 100)
lowVal = CalculateLOW.C(dose, 1)

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df70 = data.frame(d = dose, hg = Val70)
df100 = data.frame(d = dose, hg = Val100)
df195 = data.frame(d = dose, hg = Val195)
dflow = data.frame(d = dose, hg = lowVal)

p = ggplot() + theme_bw() + 
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 2) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm"))

ggsave(filename = paste(save_path, "IDER 20percent_LOW,25,70,100,195 each", ".eps"), plot = p, device = "eps")

# 0.1 each
r <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1); L <- c(25,50,75,100,150,200,300,500,750,1000)
incremental = calculateComplexId.te(r, L, d = dose, lowLET = FALSE)
SEA <- function(dose.1) Calculate.hitC(dose.1/10, 25) + Calculate.hinC(dose.1/10, 50) + Calculate.hinC(dose.1/10, 75) + Calculate.hinC(dose.1/10, 100) + CalculateLOW.C(dose.1/10, 150) + Calculate.hinC(dose.1/10, 200) + CalculateLOW.C(dose.1/10, 300) + Calculate.hinC(dose.1/10, 500) + CalculateLOW.C(dose.1/10, 750)+ CalculateLOW.C(dose.1/10, 1000)
Val25 = Calculate.hitC(dose, 25)
Val50 = Calculate.hitC(dose, 50)
Val75 = Calculate.hitC(dose, 75)
Val100 = Calculate.hitC(dose, 100)
Val150 = Calculate.hitC(dose, 150)
Val200 = Calculate.hitC(dose, 200)
Val300 = Calculate.hitC(dose, 300)
Val500 = Calculate.hitC(dose, 500)
Val750 = Calculate.hitC(dose, 750)
Val1000 = Calculate.hitC(dose, 1000)


dfsea = data.frame(d = dose, hg = SEA(dose))

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df50 = data.frame(d = dose, hg = Val50)
df75 = data.frame(d = dose, hg = Val75)
df100 = data.frame(d = dose, hg = Val100)
df150 = data.frame(d = dose, hg = Val150)
df200 = data.frame(d = dose, hg = Val200)
df300 = data.frame(d = dose, hg = Val300)
df500 = data.frame(d = dose, hg = Val500)
df750 = data.frame(d = dose, hg = Val750)
df1000 = data.frame(d = dose, hg = Val1000)

p = ggplot() + theme_bw() + 
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df50, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df75, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df150, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df200, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df300, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df500, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df750, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df1000, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dfsea, aes(x = d, y = hg), colour = "black", size = 1, linetype = "dashed") +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 1) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave(filename = paste(save_path, "IDER 10percent_LOW,25,50,75,100,150,200,300,500,750,1000 each", ".eps"), plot = p, device = "eps")

#############################################
# Different kinds confidence interval plots #
#############################################

# Some parameter initializations
sampleN = 200
mod = 1              # 1 if HIN, 0 if HIT
dose <- c(seq(0, .00001, by = 0.000001), 
          seq(.00002, .0001, by=.00001),
          seq(.0002, .001, by=.0001),
          seq(.002, .01, by=.001),
          seq(.02, 1., by=.01))

dose <- c(seq(0.0, 1., by = 0.1))

# 5x0.2 Confidence Interval Monte
confidenceIntervals = CIHelper(sampleNum = sampleN, d = dose, r = c(0.2, 0.2, 0.2, 0.2, 0.2), L = c(25, 70, 100, 180), model = 1)
naiveCI = confidenceIntervals[[1]]
monteCarloCI = confidenceIntervals[[2]]
incremental = calculateComplexId.te(r = c(0.2, 0.2, 0.2, 0.2, 0.2), L = c(25, 70, 100, 180), d = dose, lowLET = TRUE)
SEA <- function(dose.1) Calculate.hitC(dose.1/5, 25) + Calculate.hitC(dose.1/5, 70) + Calculate.hitC(dose.1/5, 100) + Calculate.hitC(dose.1/5, 180) + CalculateLOW.C(dose.1/5, 1)
Val195 = Calculate.hitC(dose, 180)
Val25 = Calculate.hitC(dose, 25)
Val70 = Calculate.hitC(dose, 70)
Val100 = Calculate.hitC(dose, 100)
lowVal = CalculateLOW.C(dose, 1)

dfsea = data.frame(d = dose, hg = SEA(dose))

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df70 = data.frame(d = dose, hg = Val70)
df100 = data.frame(d = dose, hg = Val100)
df195 = data.frame(d = dose, hg = Val195)
dflow = data.frame(d = dose, hg = lowVal)
ciPlot1 = ggplot() + theme_bw() + 
  geom_ribbon(aes(x = dose, ymin = monteCarloCI[1, ], ymax = monteCarloCI[2, ]), fill = "yellow") +
  geom_line(data = dfsea, aes(x = d, y = hg), colour = "black", size = 1, linetype = "dashed") +
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 1) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

print(ciPlot1)

ggsave(filename = paste(save_path, "ConfidenceIntervalMonteCarlo 20percent_LOW,25,70,100,195 each", ".eps"), plot = ciPlot1, device = "eps")

# 5x0.2 Confidence Interval Naive
incremental = calculateComplexId.te(r = c(0.2, 0.2, 0.2, 0.2, 0.2), L = c(25, 70, 100, 180), d = dose, lowLET = TRUE)
Val195 = Calculate.hitC(dose, 180)
Val25 = Calculate.hitC(dose, 25)
Val70 = Calculate.hitC(dose, 70)
Val100 = Calculate.hitC(dose, 100)
lowVal = CalculateLOW.C(dose, 1)

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df70 = data.frame(d = dose, hg = Val70)
df100 = data.frame(d = dose, hg = Val100)
df195 = data.frame(d = dose, hg = Val195)
dflow = data.frame(d = dose, hg = lowVal)
ciPlot2 = ggplot() + theme_bw() + 
  geom_ribbon(aes(x = dose, ymin = naiveCI[1, ], ymax = naiveCI[2, ]), fill = "yellow") +
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 1) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave(filename = paste(save_path, "ConfidenceIntervalNaive 20percent_LOW,25,70,100,195 each", ".eps"), plot = ciPlot2, device = "eps")

twoPanelArr = list(0)
twoPanelArr[[1]] = ciPlot1
twoPanelArr[[2]] = ciPlot2
multiplotHelper(twoPanelArr, "~/Desktop/plots/1*2 ConfidenceInterval 20percent_LOW,25,70,100,195 each.eps", num_col = 2, w = 8, h = 4)


# 0.05 0.8 Confidence Interval Monte
r <- c(0.05, 0.05, 0.05, 0.05, 0.8);
confidenceIntervals = CIHelper(sampleNum = sampleN, d = dose, r = c(0.05, 0.05, 0.05, 0.05, 0.8), L = c(25, 70, 100, 180), model = 0)
naiveCI = confidenceIntervals[[1]]
monteCarloCI = confidenceIntervals[[2]]
incremental = calculateComplexId.te(r = c(0.05, 0.05, 0.05, 0.05, 0.8), L = c(25, 70, 100, 180), d = dose, lowLET = TRUE)
SEA <- function(dose.1) Calculate.hitC(dose.1/20, 25) + Calculate.hitC(dose.1/20, 70) + Calculate.hitC(dose.1/20, 100) + Calculate.hitC(dose.1/20, 180) + CalculateLOW.C(4*dose.1/5, 1)
Val195 = Calculate.hitC(dose, 180)
Val25 = Calculate.hitC(dose, 25)
Val70 = Calculate.hitC(dose, 70)
Val100 = Calculate.hitC(dose, 100)
lowVal = CalculateLOW.C(dose, 1)

dfsea = data.frame(d = dose, hg = SEA(dose))

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df70 = data.frame(d = dose, hg = Val70)
df100 = data.frame(d = dose, hg = Val100)
df195 = data.frame(d = dose, hg = Val195)
dflow = data.frame(d = dose, hg = lowVal)
ciPlot1 = ggplot() + theme_bw() + 
  geom_ribbon(aes(x = dose, ymin = monteCarloCI[1, ], ymax = monteCarloCI[2, ]), fill = "yellow") +
  geom_line(data = dfsea, aes(x = d, y = hg), colour = "black", size = 1, linetype = "dashed") +
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 1) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave(filename = paste(save_path, "ConfidenceIntervalMonteCarlo 80percent_LOW, 5percent_25,70,100,195 each", ".eps"), plot = ciPlot1, device = "eps")

# 0.05 0.8 Confidence Interval Naive
incremental = calculateComplexId.te(r = c(0.05, 0.05, 0.05, 0.05, 0.8), L = c(25, 70, 100, 180), d = dose, lowLET = TRUE)
Val195 = Calculate.hitC(dose, 180)
Val25 = Calculate.hitC(dose, 25)
Val70 = Calculate.hitC(dose, 70)
Val100 = Calculate.hitC(dose, 100)
lowVal = CalculateLOW.C(dose, 1)

df1 = data.frame(d = dose, hg = incremental[, 2])
df25 = data.frame(d = dose, hg = Val25)
df70 = data.frame(d = dose, hg = Val70)
df100 = data.frame(d = dose, hg = Val100)
df195 = data.frame(d = dose, hg = Val195)
dflow = data.frame(d = dose, hg = lowVal)
ciPlot2 = ggplot() + theme_bw() + 
  geom_ribbon(aes(x = dose, ymin = naiveCI[1, ], ymax = naiveCI[2, ]), fill = "yellow") +
  geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
  geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
  geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 1) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.ticks.x = element_line(size = 1)) +
  theme(axis.ticks.y = element_line(size = 1)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave(filename = paste(save_path, "ConfidenceIntervalNaive 80percent_LOW, 5percent_25,70,100,195 each", ".eps"), plot = ciPlot2, device = "eps")

twoPanelArr[[1]] = ciPlot1
twoPanelArr[[2]] = ciPlot2
multiplotHelper(twoPanelArr, "~/Desktop/plots/1*2 ConfidenceInterval 80percent_LOW, 5percent_25,70,100,195 each.eps", num_col = 2, w = 8, h = 4)

# Zoom in Fe
p1 = plottingHelper(195, TRUE, save_path = save_path, zoom = 1)
p2 = plottingHelper(195, TRUE, save_path = save_path, zoom = 2)
twoPanelArr[[1]] = p1
twoPanelArr[[2]] = p2
multiplotHelper(twoPanelArr, "~/Desktop/plots/1*2 Fe Zoomed Blackwhite.eps", num_col = 2, w = 8, h = 4)

p1 = plottingHelper(195, FALSE, save_path = save_path, zoom = 1)
p2 = plottingHelper(195, FALSE, save_path = save_path, zoom = 2)
twoPanelArr[[1]] = p1
twoPanelArr[[2]] = p2
multiplotHelper(twoPanelArr, "~/Desktop/plots/1*2 Fe Zoomed Color.eps", num_col = 2, w = 8, h = 4)
















#################################
# Old code, leave for reference #
#################################

# # 0.05 0.8 Confidence Interval Monte
# r <- c(0.05, 0.05, 0.05, 0.05, 0.8); L <- c(25, 70, 100, 195)
# incremental = calculateComplexId(r, L, d = dose, lowLET = TRUE)
# SEA <- function(dose.1) Calculate.hinC(dose.1/20, 25) + Calculate.hinC(dose.1/20, 70) + Calculate.hinC(dose.1/20, 100) + Calculate.hinC(dose.1/20, 180) + CalculateLOW.C(4*dose.1/5, 1)
# Val195 = Calculate.hinC(dose, 180)
# Val25 = Calculate.hinC(dose, 25)
# Val70 = Calculate.hinC(dose, 70)
# Val100 = Calculate.hinC(dose, 100)
# lowVal = CalculateLOW.C(dose, 1)
# 
# dfsea = data.frame(d = dose, hg = SEA(dose))
# 
# df1 = data.frame(d = dose, hg = incremental[, 2])
# df25 = data.frame(d = dose, hg = Val25)
# df70 = data.frame(d = dose, hg = Val70)
# df100 = data.frame(d = dose, hg = Val100)
# df195 = data.frame(d = dose, hg = Val195)
# dflow = data.frame(d = dose, hg = lowVal)
# ciPlot1 = ggplot() + theme_bw() + 
#   geom_ribbon(aes(x = dose, ymin = monteCarloCI[1, ], ymax = monteCarloCI[2, ]), fill = "yellow") +
#   geom_line(data = dfsea, aes(x = d, y = hg), colour = "black", size = 1, linetype = "dashed") +
#   geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
#   geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
#   geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
#   geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
#   geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
#   geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 1) +
#   theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
#   theme(axis.ticks.x = element_line(size = 1), axis.text.x=element_blank()) +
#   theme(axis.ticks.y = element_line(size = 1), axis.text.y=element_blank()) +
#   theme(axis.ticks.length=unit(0.4,"cm")) +
#   theme(plot.margin=unit(c(1,1,1,1),"cm"))
# print(ciPlot1)
# 
# # 0.05 0.8 Confidence Interval Naive
# r <- c(0.05, 0.05, 0.05, 0.05, 0.8); L <- c(25, 70, 100, 195)
# incremental = calculateComplexId(r, L, d = dose, lowLET = TRUE)
# Val195 = Calculate.hinC(dose, 180)
# Val25 = Calculate.hinC(dose, 25)
# Val70 = Calculate.hinC(dose, 70)
# Val100 = Calculate.hinC(dose, 100)
# lowVal = CalculateLOW.C(dose, 1)
# 
# df1 = data.frame(d = dose, hg = incremental[, 2])
# df25 = data.frame(d = dose, hg = Val25)
# df70 = data.frame(d = dose, hg = Val70)
# df100 = data.frame(d = dose, hg = Val100)
# df195 = data.frame(d = dose, hg = Val195)
# dflow = data.frame(d = dose, hg = lowVal)
# ciPlot2 = ggplot() + theme_bw() + 
#   geom_ribbon(aes(x = dose, ymin = naiveCI[1, ], ymax = naiveCI[2, ]), fill = "yellow") +
#   geom_line(data = df25, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
#   geom_line(data = df70, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
#   geom_line(data = df100, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
#   geom_line(data = df195, aes(x = d, y = hg), colour = "#55aaff", size = 1) + 
#   geom_line(data = dflow, aes(x = d, y = hg), colour = "#55aaff", size = 1) +
#   geom_line(data = df1, aes(x = d, y = hg), colour = "red", size = 1) +
#   theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
#   theme(axis.ticks.x = element_line(size = 1), axis.text.x=element_blank()) +
#   theme(axis.ticks.y = element_line(size = 1), axis.text.y=element_blank()) +
#   theme(axis.ticks.length=unit(0.4,"cm")) +
#   theme(plot.margin=unit(c(1,1,1,1),"cm"))
# print(ciPlot2)
# 
# dfChange = data.frame(d = dose, h = (naiveCI[2, ] - naiveCI[1, ])/(monteCarloCI[2, ] - monteCarloCI[1, ]))
# changePlot = ggplot() + theme_bw() +
#   geom_line(data = dfChange, aes(x = d, y = h), colour = "#55aaff", size = 1) +
#   theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
#   theme(axis.ticks.x = element_line(size = 1)) +
#   theme(axis.ticks.y = element_line(size = 1)) +
#   theme(axis.ticks.length=unit(0.4,"cm")) +
#   theme(plot.margin=unit(c(1,1,1,1),"cm"))
# 
# saveRDS(monteCarloCI, file="~/Dropbox/0.2*5_monteCarlo200.Rda")
# monteCarloCI1000 = readRDS(file="~/Dropbox/monteCarlo1000.Rda")
# 
# monteCarloCI1000 <- readRDS(file="~/Desktop/monteCarlo1000.Rda")
# 
# print(changePlot)



# if (Lval == 1) {
#   lowVal = CalculateLOW.C(dose, Lval)
#   df1 = data.frame(d = dose, low = lowVal)
#   df2 = data.frame(x = d, hg = hgArr, nw = nWeight)
#   dfProton = data.frame(x = dProton, hg = protonArr, nw = nWeightProton)
#   
#   p = ggplot() + theme_bw() + 
#     geom_line(data = df1, aes(x = dose, y = low), colour = color1, size = 1) + 
#     geom_segment(data = df2, size = 0.8, aes(x = d, xend = d, y = hg - 1/sqrt(nWeight), yend = hg + 1/sqrt(nWeight)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) + 
#     geom_segment(data = df2, size = 0.8, aes(x = d, xend = d, y = hg, yend = hg - 1/sqrt(nWeight)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) +
#     geom_segment(data = dfProton, size = 0.8, aes(x = dProton, xend = dProton, y = hg - 1/sqrt(nWeightProton), yend = hg + 1/sqrt(nWeightProton)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) + 
#     geom_segment(data = dfProton, size = 0.8, aes(x = dProton, xend = dProton, y = hg, yend = hg - 1/sqrt(nWeightProton)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) +
#     theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
#     theme(axis.ticks.x = element_line(size = 1), axis.text.x=element_blank()) +
#     theme(axis.ticks.y = element_line(size = 1), axis.text.y=element_blank()) +
#     theme(axis.ticks.length=unit(0.4,"cm")) +
#     theme(plot.margin = unit(c(1,1.5,1,1.5), "cm")) +
#     theme(panel.border = border_custom())
#   p = last_plot() + geom_point(data = df2, aes(x = d, y = hg), size = 4, fill ="white", shape = 21)
#   p = last_plot() + geom_point(data = dfProton, aes(x = dProton, y = hg), size = 4, fill ="black", shape = 21)
#   if (blackwhite) {
#     ggsave(filename = paste(save_path, toString(Lval), "Blackwhite", ".eps"), plot = p, device = "eps")
#   } else {
#     ggsave(filename = paste(save_path, toString(Lval), "Color", ".eps"), plot = p, device = "eps")
#   }
#   return(p)
# }
# 
# if (Lval == 193 || Lval == 195) {
#   hinVal = Calculate.hinC(dose, 180)
#   hitVal = Calculate.hitC(dose, 180)
# } else {
#   hinVal = Calculate.hinC(dose, Lval)
#   hitVal = Calculate.hitC(dose, Lval)
# }
# 
# df1 = data.frame(d = dose, hgNTE = hinVal, hgTE = hitVal)
# df2 = data.frame(x = d, hg = hgArr, nw = nWeight)
# p = ggplot() + theme_bw() + 
#        geom_line(data = df1, aes(x = dose, y = hgNTE), colour = color1, size = 1) + 
#        geom_line(data = df1, aes(x = dose, y = hgTE), colour = color2, size = 1, linetype = "dashed") + 
#        geom_segment(data = df2, size = 0.8, aes(x = d, xend = d, y = hg - 1/sqrt(nWeight), yend = hg + 1/sqrt(nWeight)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) + 
#        geom_segment(data = df2, size = 0.8, aes(x = d, xend = d, y = hg, yend = hg - 1/sqrt(nWeight)), arrow = arrow(angle = 90, length = unit(0.2,"cm"))) +
#        theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank()) +
#        theme(axis.ticks.length=unit(0.4,"cm"))
# 
# if (has_axis_text) {
#   p = p + theme(axis.ticks.x = element_line(size = 1), axis.text.x=element_blank()) +
#           theme(axis.ticks.y = element_line(size = 1), axis.text.y=element_blank())
# } else {
#   p = p + theme(axis.ticks.x = element_line(size = 1)) + 
#           theme(axis.ticks.y = element_line(size = 1))
# }
# 
# if (!zoom) {
#   p = p + theme(plot.margin = unit(c(1,1.5,1,1.5), "cm"))
# } else {
#   p = p + theme(plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))
# }
# 
# if (Lval != 100 && Lval != 953) {
#   if (!zoom) {
#     p = last_plot() + theme(panel.border = border_custom())
#   } else if (zoom == 1) {
#     p = last_plot() + theme(panel.border = border_custom())
#   }
# }
# 
# p = last_plot() + geom_point(data = df2, aes(x = d, y = hg), size = 2, fill ="black", shape = 21)
# if (blackwhite) {
#   ggsave(filename = paste(save_path, toString(Lval), "Blackwhite", ".eps"), plot = p, device = "eps")
# } else {
#   ggsave(filename = paste(save_path, toString(Lval), "Color", ".eps"), plot = p, device = "eps")
# }
# return (p)



# ############################### Native plot ###############################
# Lval = 100
# d = c()
# nWeight = c()
# hgArr = c()
# for (row in 1:nrow(dfr)) {
#   if (dfr[row, "L"] == Lval) {
#     d = c(d, dfr[row, "dose.1"])
#     nWeight = c(nWeight, dfr[row, "NWeight"])
#     hgArr = c(hgArr, dfr[row, "HG"])
#   }
# }
# dose <- c(seq(0, .00001, by = 0.000001),
#           seq(.00002, .0001, by=.00001),
#           seq(.0002, .001, by=.0001),
#           seq(.002, .01, by=.001),
#           seq(.02, d[length(d)], by=.01))
# hinVal = Calculate.hinC(dose, Lval)
# hitVal = Calculate.hitC(dose, Lval)
# 
# 
# plot(x = dose, y = hinVal, type='l', col='red', bty='l', ann='F', xaxt = "n", yaxt = "n", lwd = 2)
# axis(side = 1, at = seq(0, d[length(d)], length.out = 10), c(0, rep("", 8), trunc(10*d[length(d)])/10), las = 1)
# axis(side = 2, at = seq(0, hinVal[length(hinVal)] + 1/sqrt(nWeight[1]), length.out = 5), c(0, rep("", 3), trunc(10*(hinVal[length(hinVal)] + 1/sqrt(nWeight[1])))/10), las = 1)
# lines(x = dose, y = hitVal, type='l', col='green', bty='l', ann='F', lwd = 2)
# points(d, hgArr, pch = 19)
# for (i in 1:length(nWeight)) {
#   arrows(d[i], hgArr[i] - 1/sqrt(nWeight[i]), d[i], hgArr[i] + 1/sqrt(nWeight[i]), code = 3, length = 0.04, angle = 90, lwd = 2)
# }