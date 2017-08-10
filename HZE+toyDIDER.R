#Free, open-source software under GNU GPLv3. It comes with no warranty. It illustrates the difference between an IDER and DIDER using actual data for one component of a mixture and simulated data for the component which has a DIDER
#Written by Mark Ebert, Edward Huang, Dae Woong Ham and Ray Sachs, summer 2017. This script will be used in IJRB paper
#Relevant references and their abbreviations in commenting this script are the following.
#".93Alp"=Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-particle radiations." Rad Res 136:382-391 (1993)
#".94Alp"=Alpen et al. "Fluence-based relative biological effectiveness for charged particle carcinogenesis in mouse Harderian gland." Adv Space Res 14(10): 573-581. (1994).  
#"16Chang"=Chang et al. "Harderian Gland Tumorigenesis: Low-Dose and LET Response." Radiat Res 185(5): 449-460. (2016).  
#"16Srn"=Siranart et al."Mixed Beam Murine Harderian Gland Tumorigenesis: Predicted Dose-Effect Relationships if neither Synergism nor Antagonism Occurs." Radiat Res 186(6): 577-591 (2016).  
#"17Cuc"=Cucinotta & Cacao. "Non-Targeted Effects Models Predict Significantly Higher Mars Mission Cancer Risk than Targeted Effects Models." Sci Rep 7(1): 1832. (2017). PMC5431989

library(deSolve) # library for solving differential equations
library(minpack.lm) #for non-linear regression package
rm(list=ls())
dfr=data.frame( #Includes all HZE data used in 16Chang  
  Dose=c(0.2,0.4,0.6,1.2,2.4,3.2,5.1,7,0.05,0.1,0.15,0.2,0.4,0.8,1.6,0.05,0.1,0.2,0.4,0,0.1,0.2,0.4,0.8,1.6,0.4,0.8,1.6,3.2,0.05,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.04,0.08,0.16,0.32,0.033,0.066,0.13,0.26,0.52,.2, .4, .6),
  HG=c(0.091,0.045,0.101,0.169,0.347,0.431,0.667,0.623,0.156,0.215,0.232,0.307,0.325,0.554,0.649,0.123,0.145,0.207,0.31,0.026,0.083,0.25,0.39,0.438,0.424,0.093,0.195,0.302,0.292,0.109,0.054,0.066,0.128,0.286,0.183,0.167,0.396,0.536,0.192,0.234,0.317,0.092,0.131,0.124,0.297,0.082,0.088,0.146,0.236,0.371,.154,.132,.333),#HG prevalence as defined in 16Chang
  NWeight=c(520,2048,1145,584,313,232,293,221,1162,877,455,409,374,223,320,742,661,347,131,6081,1091,251,244,191,131,645,255,199,111,649,378,973,833,201,468,381,197,109,496,257,185,1902,1063,884,350,1767,1408,874,299,261,322,206,67),#nominal weight for weighted least squares regression; see .93Alp. The Lanthanum entries were obtained by measuring the main graph in 17Cuc 
  index=c(rep(1,8),rep(0,17), rep(1,4),  rep(0,24)),#index=0 for Z>3 ions, index=1 for proton p (=1H1) and 2HE4 ions
  L=c(rep(1.6,8),rep(193,7),rep(250,4),rep(195,6),rep(0.4,4),rep(25,5),rep(464,4),rep(193,3),rep(70,4),rep(100,5),rep(953,3)), #L=LET=LET_infinity=stopping power
  Z=c(rep(2,8),rep(26,17),rep(1,4),rep(10,5),rep(43,4),rep(26,3),rep(14,4),rep(22,5),rep(57,3)))#proton #, e.g. 2 for 2He4
dfra=dfr[c(1:19,26:53),] ##removes the zero dose case and the no isograft data
dfrHZE=subset(dfra,Z>3)#for heavy ions
## Our new HZE model; uses 3 adjustable parameters. In next line, phi controls how fast NTE build up from zero; phi is not needed during calibration since 150*phi*Dose/L>>1 at every observed Dose !=0. phi is needed for later synergy calculations and is so large one has to be careful in R
phi=3e3 # even bigger phi would give essentially the same final results but cause extra problems with R  
HZEm=nls(HG~.0275+(1-exp(-0.01*(aa1*L*Dose*exp(-aa2*L) +(1-exp(-150*phi*Dose/L))*kk1))),data=dfrHZE, weights=NWeight, start=list(aa1=.9,aa2=.01,kk1=0.048))#calibrating parameters;
HZEc=coef(HZEm)#calibrated central values of the 3 parameters
HZEC=function(Dose,L) 1-exp(-0.01*(HZEc[1]*L*Dose*exp(-HZEc[2]*L)+(1-exp(-150*phi*Dose/L))*HZEc[3]))#calibrated IDER. The above equations imply: HZEC(dose=0)=0; HZEC increases monotonically and apporaches 1 at large doses.  
summary(HZEm,cor=T)# All 3 adjustable parameters are signficant at the 1e-5 level, an unusually good result

#Dae
low_LET_IDER = function(lambda, beta, d) {
  beta*(1-exp(-lambda*d))
}

d2 = 0.01*0:300
set.seed(15)
simulated_data = low_LET_IDER(lambda = 2, beta = 0.5, d = d2) + 0.05*rweibull(301, 2, 1)
simulated_df = data.frame(d = d2, HG = simulated_data)
simulated_calibration = nls(HG ~ low_LET_IDER(lambda, beta, d),data= simulated_df, start=list(lambda = 2, beta = 0.5))
summary(simulated_calibration)
coefs = coef(simulated_calibration)

#MIXDER
dE_1 = function(d, aa1, aa2, kk1, phi, L) {
  ((150*kk1*phi*exp(-150*phi*d/L)/L + aa1*L*exp(-aa2*L))*exp(-0.01*(kk1*(1-exp(-150*phi*d/L)) + aa1*L*exp(-aa2*L)*d)))/100
}

dE_2 = function(lambda, beta, E) {
  lambda*(beta-E)
}

HZE_MIXDER = function(r, L = 193, aa1 = HZEc[1], aa2 = HZEc[2], kk1 = HZEc[3], phi = 3e3, lambda = coefs[1], beta = coefs[2], d) {
  dE=function(yini,State,Pars){
    aa1 = aa1; aa2 = aa2; kk1 = kk1; lambda = lambda; beta = beta; phi = phi; L = L
    with(as.list(c(State, Pars)), {
      u = vector(length = 1)
      u[1] = uniroot(function(d) 1-exp(-0.01*(aa1*L*d*exp(-aa2*L)+(1-exp(-150*phi*d/L))*kk1)) - I, lower = 0, upper = 10, extendInt = "yes", tol = 10^-10)$root 
      dI = vector(length = 2)
      dI[1] = r[1]*dE_1(d = u[1], aa1 = aa1, aa2 = aa2, kk1 = kk1, phi = phi, L = L)
      dI[2] = r[2]*dE_2(E = I, lambda = lambda, beta = beta)
      dI = sum(dI)
      return(list(c(dI)))
    })
  }
  pars = NULL; yini = c(I= 0); d = d
  out = ode(yini, times = d, dE, pars, method = "radau")
  return(out)
}
d = seq(0, 9.9, 0.01)
HZE_MIXDER(r = c(0.5,0.5), d = d)[, 2]

dE_IDER = function(d, lambda, beta) {
  beta*lambda*exp(-lambda*d)
}

HZE_MIXDER_IDER = function(r, L = 193, aa1 = HZEc[1], aa2 = HZEc[2], kk1 = HZEc[3], phi = 3e3, lambda = coefs[1], beta = coefs[2], d) {
  dE=function(yini,State,Pars){
    aa1 = aa1; aa2 = aa2; kk1 = kk1; lambda = lambda; beta = beta; phi = phi; L = L
    with(as.list(c(State, Pars)), {
      u = vector(length = 2)
      u[1] = uniroot(function(d) 1-exp(-0.01*(aa1*L*d*exp(-aa2*L)+(1-exp(-150*phi*d/L))*kk1)) - I, lower = 0, upper = 10, extendInt = "yes", tol = 10^-10)$root 
      if (I < beta) {
        u[2] = uniroot(function(d) beta*(1-exp(-lambda*d)) - I, lower = 0, upper = 5, extendInt = "yes", tol = 10^-10)$root
      }
      else {
        u[2] = 0
        r[2] = 0
      }
      dI = vector(length = 2)
      dI[1] = r[1]*dE_1(d = u[1], aa1 = aa1, aa2 = aa2, kk1 = kk1, phi = phi, L = L)
      dI[2] = r[2]*dE_IDER(d = u[2], lambda = lambda, beta = beta)
      dI = sum(dI)
      return(list(c(dI)))
    })
  }
  pars = NULL; yini = c(I= 0); d = d
  out = ode(yini, times = d, dE, pars, method = "radau")
  return(out)
}

HZE_MIXDER_IDER(r = c(0.5,0.5), d = d)
plot(x = d, y = HZEC(Dose = d, L = 173), col = "purple", type = "l",lty=2,lwd=2, ylab="HG", xlab="dose")
lines(x = d, y = HZE_MIXDER(r = c(0.5,0.5), d = d)[ ,2], col = "blue")
lines(x = d, y = HZE_MIXDER_IDER(r = c(0.5,0.5), d = d)[, 2], col = "green")
abline(h=coef(simulated_calibration)["beta"])
lines(x = d, y = low_LET_IDER(beta = coefs[2], lambda = coefs[1], d = d), col = "red",lty=2, lwd=2)
