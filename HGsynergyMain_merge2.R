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

####### DATA #######
dfr <- data.frame( # data used in 16Chang; includes data analyzed in .93Alp and .94Alp  
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
  comments = c(".93AlpLooksOK", rep("", 7), ".93AlplooksOK", rep('', 11), ".93Alp.no.iso", "not in 17Cuc (or 16Chang?)", rep("", 3), "16Chang all OK?", rep('', 24), ".94Alp","From graphs",'e.g. in 17Cuc')
) 

# Data for HG induced by photons from Cs-137 or Co-60 beta decay; from 16Chang (and calibration of LQ model)
ddd <- data.frame(dose.1 = c(0, 0.4, 0.8, 1.6, 3.2, 7, 0, .4, .8, .12, 1.6),
                  HG = c(.026, .048, .093, .137, .322, .462, .0497, .054, .067, .128, .202),
                  NWeight = c(6081.2, 4989.5, 1896.8, 981.1, 522.2, 205.2, 7474.1, 2877.6, 1423.7, 689.9, 514.9),
                  Nucleus = c(rep("Cobalt-60", 6), rep("Cesium-137", 5)),
                  Comments = c(rep("TBD", 11))
)
Y <- 0.001 * dfr[, "MeVperu"] #  convert to GeV/u for convenience in a calculation
dfr[, "Katz"] <- round((dfr[, "Z"]) ^2 * (2.57 * Y^2 + 4.781 * Y + 2.233) / (2.57 * Y^2 + 4.781 * Y), 2) #  special relativistic calculation of Z^2/beta^2. The numerics include conversion from GeV to joules and from u to kg.
dfr[, "beta"] <- round(dfr[, "Z"] * sqrt(1 / dfr[, "Katz"]), 3) #  i.e. Z*sqrt(beta^2/Z^2) 
dfr[, "Zeff"] <- round(dfr[, "Z"] * (1 - exp( -125 * dfr[, "Z"] ^ (-2.0 / 3))), 2) #  Barkas formula for Zeff; for us Zeff is almost Z

dfra <- dfr[c(1:19,26:53),] #  removes the zero dose case and the no isograft data

#####  photon and HZE/NTE Models #####
LQ <- lm(HG ~ dose.1 + I(dose.1 ^ 2), data = ddd) # linear model fit on ddd dataset
summary(LQ, correlation = T) 

#HERE OUR NEW HZE/NTE MODEL; works well. Uses 3 adjustable parameters. There is also an HZE/TE model and NTE or TE models for Z<=3, but for a while we now use the HZE/NTE model for Z>3 with data for Z>=8  
dfrHZE <- subset(dfra, Z > 3) # look only at HZE not at much lower LET ions. # In next line phi controls how fast NTE build up from zero; not really needed during calibration since 150*phi*Dose/L>>1 at every observed Dose !=0. phi needed for later synergy calculations.
phi <- 1000 # Even larger phi should give the same final results, but might cause extra problems with R. 
HZEm <- nls(HG ~ .0275 + (1 - exp ( -0.01 * (aa1 * L * dose.1 * exp( -aa2 * L) + (1 - exp( -150 * phi * dose.1 / L)) * kk1))), #  calibrating parameters in a model that modifies the hazard function NTE models in 17Cuc.
            data = dfrHZE, 
            weights = NWeight,
            start = list(aa1 = .9, aa2 = .01, kk1 = 6)) 
summary(HZEm, correlation = T);vcov(HZEm)# parameter values & accuracy; variance-covariance matrix
NTE.HZE.c <- coef(HZEm) #  calibrated central values of the 3 parameters. Next is the IDER, =0 at dose 0
HHC <- function(dose.1,L) 0.01*(NTE.HZE.c[1]*L*dose.1*exp(-NTE.HZE.c[2]*L)+(1-exp(-150*phi*dose.1/L))*NTE.HZE.c[3]) # calibrated hazard function

CalculateHZEC <- function(dose.1, L) {
  1-exp(-HHC(dose.1,L)) #  Calibrated IDER
}

dose <- c(seq(0, .00001, by = 0.000001),seq(.00002,.0001,by=.00001),seq(.0002,.001,by=.0001),seq(.002,.01,by=.001),seq(.02,.5,by=.01))##look carefully near zero, but go out to 0.5 Gy
#dose=dose[1:30] #this can be used to zoom in on the very low dose behavior in the graphs
####### calculate baseline MIXDER I(d) for mixtures of HZE components modeled by NTE IDER #######
IntegrateNTE_HZE_IMIXDER <- function(r, L, d = dose, aa1 = NTE.HZE.c[1], aa2 = NTE.HZE.c[2], kk1 = NTE.HZE.c[3]) {
  dE <- function(yini, State, Pars) {
    aa1 <- aa1; aa2 <- aa2; kk1 <- kk1
    with(as.list(c(State, Pars)), {
      aa = vector(length = length(L))
      u = vector(length = length(L))
      for (i in 1:length(L)) {
        aa[i] = aa1*L[i]*exp(-aa2*L[i])
        u[i] = uniroot(function(d) 1-exp(-0.01*(aa[i]*d+(1-exp(-150*phi*d/L[i]))*kk1)) - I, lower = 0, upper = 20, tol = 10^-10)$root
      }
      dI = vector(length = length(L))
      for (i in 1:length(L)) {
        dI[i] = r[i]*0.01*(aa[i]+exp(-150*phi*u[i]/L[i])*kk1*150*phi/L[i])*exp(-0.01*(aa[i]*u[i]+(1-exp(-150*phi*u[i]/L[i]))*kk1))
        
      }
      dI = sum(dI)
      return(list(c(dI)))
    })
  }
  pars = NULL; yini = c(I= 0); d = d
  out = ode(yini, times = d, dE, pars, method = "radau")
  return(out)
} 

####### Plotting I(d) example

# another example
# r <- .25*1:4; L <- c(25,70,190,250);plot(INCRL, type='l',ylim=c(0,.5),col='red',bty='l',ann='F')
# lines(dose,CalculateHZEC(dose,190), col='green')#component 3
# lines(dose,CalculateHZEC(dose,250),col='green')#component 4
# lines(dose,CalculateHZEC(dose,70),col='green')#component 2
# lines(dose,CalculateHZEC(dose,25),col='green') #component 1
# SEA <- function(dose.1) CalculateHZEC(dose.1/4, 25) + CalculateHZEC(dose.1/4, 70) + CalculateHZEC(dose/4, 190) + CalculateHZEC(dose.1/3, 250)
# lines(dose, SEA(dose), lty=2)

########### Light ion Z <= 3 model 1 ######### 
dfrL <- subset(dfra, Z <= 3) #for Light ions
LOW.m <- nls(HG ~ .0275 + 1-exp(-bet * dose.1),
             data = dfrL,
             weights = NWeight,
             start = list(bet = .5))
summary(LOW.m)
LOW.c <- coef(LOW.m)  # calibrated central values of the parameter
CalculateLOW.C <- function(dose.1, L) { # Calibrated Low LET model. Use L=0, but maybe later will use L >0 but small 
  return(1-exp(-LOW.c[1] * dose.1))
}  

#######Next: visual checks to see if our calibration is consistent with 16Chang, .93Alp, .94Alp and 17Cuc
## Put various values in our calibrated model to check with numbers and graphs in these references
#  L=193; dose.1 = dfrHZE[1:7, "dose.1"]; HGe = dfr[1:7,"HG"] # same for Fe
plot(c(0, 7), c(0, 1), col = 'red', ann = 'F') 
ddose <- 0.01 * 0:700; lines(ddose, CalculateLOW.C(ddose, 0) + .0275)  #calibrated lowLET IDER
points(dfrL[1:8, "dose.1"], dfrL[1:8,"HG"],pch=19) #RKS: Helium data points
points(dfrL[9:12, "dose.1"], dfrL[9:12, "HG"] )  # proton data points 

#####MIXDER# RKS next is just a  fossil I think
dE_1 <- function(d, aa1, aa2, kk1, phi, L) {
   ((150 * kk1 * phi * exp(-150 * phi * d / L) / L + aa1 * L * exp( -aa2 * L)) * 
  exp(-0.01 * (kk1 * (1 - exp( - 150 * phi * d / L)) + aa1 * L * exp(-aa2 * L) * d))) / 100
 }

dE_2 <- function(dose,L) { 
 LOW.c*exp(-LOW.c*dose)  
}

# r1 <- .2; r <- c(r1, 1 - r1) #Proportions. Next plot IDERs and MIXDER
# d <- .01 * 0:300.
# plot(x = d, y = CalculateHZEC(dose.1 = d, L = 173), type = "l", xlab="dose",ylab="HG",bty='l',col='green',lwd=2)
# lines(x = d, y = CalculateLOW.C( d,0), col='green',lwd=2)
# lines(x = d, y = IntegrateNTE_HZE_LOW_IMIXDER(r = r, d = d)[, 2], col = "red", lwd=2) # I(d)

# To be done next. HZE NTE MIXDER 95% CI (Edward)!!! Information criteria and 
# compare with 17Cuc (Mark)! clean up style
# Overall: add TE HZE models; decide on low LET models after seeing Mark's IC 
# results. Decide how to handle the various branches -- probably one ODE for NTE 
# HZE MIXDERs that may have on low LET component and one ODE for TE HZE MIXDERs ditto

#==============================================#
#==========Confidence Interval Part============#
#==============================================#

# Set the pseudorandom seed
set.seed(1)

Generate_CI = function(N = 500, intervalLength = 0.95, d, r, L, HZEmodel = HZEm, method = 0) {
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
      valueArr = c(valueArr, IntegrateNTE_HZE_IMIXDER(r = r, L = L, d = c(0, d), aa1 = monteCarloSamples[, 1][i], aa2 = monteCarloSamples[, 2][i], kk1 = monteCarloSamples[, 3][i])[, 2][2])
    }
    valueArr = sort(valueArr)
    
    # Returning resulting CI
    return (c(valueArr[(1-intervalLength)/2*500], valueArr[(intervalLength + (1-intervalLength)/2)*500]))
  } else {
    #========= Naive =========#
    stdErrArr = summary(HZEmodel)$coefficients[, "Std. Error"] 
    meanArr = summary(HZEmodel)$coefficients[, "Estimate"] 
    upper = IntegrateNTE_HZE_IMIXDER(r = r, L = L, d = c(0, d), aa1 = meanArr["aa1"] + 2*stdErrArr["aa1"], aa2 = meanArr["aa2"] + 2*stdErrArr["aa2"], kk1 = meanArr["kk1"] + 2*stdErrArr["kk1"])[, 2][2]
    lower = IntegrateNTE_HZE_IMIXDER(r = r, L = L, d = c(0, d), aa1 = meanArr["aa1"] - 2*stdErrArr["aa1"], aa2 = meanArr["aa2"] - 2*stdErrArr["aa2"], kk1 = meanArr["kk1"] - 2*stdErrArr["kk1"])[, 2][2]
    return (c(lower, upper))
  }
}

# Parameter initialization
r <- c(1/3, 1/3, 1/3); L <- c(25, 70, 250)
mixderCurve = IntegrateNTE_HZE_IMIXDER(r, L)
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


#egh :: ASSIGNMENT: Calculate I(d) for N > 1 HZE IDERs and one low-LET IDER
calculateComplexId <- function(r, L, d, aa1 = NTE.HZE.c[1], aa2 = NTE.HZE.c[2], kk1 = NTE.HZE.c[3], phi = 3e3, beta = LOW.c, lowLET = FALSE) {
  ## FUNCTION DESCRIPTION
  # Calulates the function I(d) from N > 1 HZE IDERs and one low-LET IDER using 
  # incremental effect additivity
  #
  # new argument: lowLET (FALSE by default, TRUE when one IDER is low-LET)
  dE <- function(yini, State, Pars) { #  Constructing an ode from the IDERS
    aa1 <- aa1; aa2 <- aa2; kk1 <- kk1; beta <- beta; phi <- phi; L <- L
    with(as.list(c(State, Pars)), {
      aa <- vector(length = length(L))  
      u <- vector(length = length(L))  
      for (i in 1:length(L)) {
        aa[i] <- aa1 * L[i] * exp(-aa2 * L[i])
        u[i] <- uniroot(function(d) 1-exp(-0.01*(aa1*L[i]*d*exp(-aa2*L[i])+(1-exp(-150*phi*d/L[i]))*kk1)) - I, lower = 0, upper = 200, extendInt = "yes", tol = 10^-10)$root #egh this is used in the single HZE and lowLET example
      }
      dI <- vector(length = length(L))
      for (i in 1:length(L)) {
        dI[i] <- r[i] * 0.01*(aa[i]+exp(-150*phi*u[i]/L[i])*kk1*150*phi/L[i])*exp(-0.01*(aa[i]*u[i]+(1-exp(-150*phi*u[i]/L[i]))*kk1))
      }
      if (lowLET == TRUE) { # If low-LET IDER is present then include it at the end of the dI vector
        u[length(L) + 1] <- uniroot(function(d) 1-exp(-beta*d) - I, lower = 0, upper = 200, extendInt = "yes", tol = 10^-10)$root
        dI[length(L) + 1] <- r[length(r)] * dE_2(d = u[length(L) + 1], L = 0)
      }
      dI <- sum(dI)
      return(list(c(dI)))
      })
  }
  return(ode(c(I = 0), times = d, dE, parms = NULL, method = "radau")) #  Finds the solution I(d) of the differential equation dE
}

r1 <- .2; r <- c(r1, 1 - r1) #Proportions. Next plot IDERs and MIXDER
d <- .01 * 0:300.
plot(x = d, y = CalculateHZEC(dose.1 = d, L = 173), type = "l", xlab="dose",ylab="HG",bty='l',col='green',lwd=2)
lines(x = d, y = CalculateLOW.C( d,0), col='green', lwd=2)
lines(x = d, y = calculateComplexId(r = r, L = 193, d = d, lowLET = TRUE)[, 2], col = "red", lwd=2) # I(d)

r <- c(1/3, 1/3, 1/3); L <- c(25, 70, 250)
INCRL <- calculateComplexId(r, L, d = dose) # incremental effect additivity 
plot(INCRL, type='l', col='red', bty='l', ann='F') # ,ylim=c(0,.4) ; I(d) plot
lines(dose, CalculateHZEC(dose, 250), col='green') # component 3
lines(dose, CalculateHZEC(dose, 70), col='green') # component 2
lines(dose, CalculateHZEC(dose, 25), col='green') # component 1
SEA <- function(dose.1) CalculateHZEC(dose.1/3, 25) + CalculateHZEC(dose.1/3, 70) + CalculateHZEC(dose.1/3, 250)
lines(dose, SEA(dose), lty=2)

d <- seq(0, .01, .0005)
r <- c(1/20, 1/20, 9/10); L <- c(70, 173)
plot(x = d, y = CalculateHZEC(dose.1 = d, L = 173), type = "l", xlab="dose",ylab="HG",bty='l',col='green',lwd=2)
lines(x = d, y = CalculateHZEC(d, 70), col='green', lwd=2) # component 3
lines(x = d, y = CalculateLOW.C(d, 0), col='green', lwd=2)
lines(x = d, y = calculateComplexId(r, L, d = d, lowLET = TRUE)[, 2], col = 'red', lwd = 2)
