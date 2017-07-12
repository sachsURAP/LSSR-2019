# This is free, open-source software under GNU GPLv3. It comes with no warranty. 
# It concerns radiogenic mouse HG tumorigenesis.
# Written by Mark Ebert, Edward Huang, and Ray Sachs, summer 2017.

# Relevant references and their abbreviations in commenting this script are the following.
# ".93Alp"=Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-particle radiations." Rad Res 136:382-391.(1993)
# ".94Alp"=Alpen et al. "Fluence-based relative biological effectiveness for charged particle carcinogenesis in mouse Harderian gland." Adv Space Res 14(10): 573-581. (1994).  
# "16Chang"=Chang et al. "Harderian Gland Tumorigenesis: Low-Dose and LET Response." Radiat Res 185(5): 449-460. (2016).  
# "16Srn"=1.	Siranart et al."Mixed Beam Murine Harderian Gland Tumorigenesis: Predicted Dose-Effect Relationships if neither Synergism nor Antagonism Occurs." Radiat Res 186(6): 577-591. (2016).  
# "17Cuc"=Cucinotta & Cacao. "Non-Targeted Effects Models Predict Significantly Higher Mars Mission Cancer Risk than Targeted Effects Models." Sci Rep 7(1): 1832. (2017).  PMC5431989

## GENERAL NOTES ON IMPROVING THE SCRIPT - EDWARD##
#
#  Use '<-' to assign variables instead of '='. This is mostly convention,
#  but there is a small technical difference. See https://stackoverflow.com/questions/1741820/assignment-operators-in-r-and 
#  for more detail.
#
#  To maximize readability/organization in our code, we should follow a 
#  set of conventions, e.g. a style guide. One common convention is to keep
#  lines under 80 characters, with a few exceptions, like our dfr() call. 
#  Google has written a set of style guides for many coding languages that
#  accurately reflect the style conventions of most of their users.
#  I recommend reading and mostly adhering to Google's guide on R, 
#  located here: https://google.github.io/styleguide/Rguide.xml
#
#  If other researchers will be poring over our code, it would be helpful to 
#  variable names that are both concise and descriptive. Our variable names 
#  are currently concise but also a bit esoteric to any outsider.
###################################################
 
rm(list=ls()) # EDWARD: this is used to clear the current environment, I recommend
              # removing it from the script and running it as needed in your console.
dfr <- data.frame( # data used in 16Chang  
  Dose = c(0.2,0.4,0.6,1.2,2.4,3.2,5.1,7,0.05,0.1,0.15,0.2,0.4,0.8,1.6,0.05,0.1,0.2,0.4,0,0.1,0.2,0.4,0.8,1.6,0.4,0.8,1.6,3.2,0.05,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.04,0.08,0.16,0.32,0.033,0.066,0.13,0.26,0.52,.2, .4, .6),
  HG = c(0.091,0.045,0.101,0.169,0.347,0.431,0.667,0.623,0.156,0.215,0.232,0.307,0.325,0.554,0.649,0.123,0.145,0.207,0.31,0.026,0.083,0.25,0.39,0.438,0.424,0.093,0.195,0.302,0.292,0.109,0.054,0.066,0.128,0.286,0.183,0.167,0.396,0.536,0.192,0.234,0.317,0.092,0.131,0.124,0.297,0.082,0.088,0.146,0.236,0.371,.154,.132,.333),#HG prevalence as defined in 16Chang
  NWeight = c(520,2048,1145,584,313,232,293,221,1162,877,455,409,374,223,320,742,661,347,131,6081,1091,251,244,191,131,645,255,199,111,649,378,973,833,201,468,381,197,109,496,257,185,1902,1063,884,350,1767,1408,874,299,261,322,206,67),
  #nominal weight for weighted least squaresregression; see .93Alp. The Lanthanum entries were obtained by measuring the main graph in 17Cuc 
  #index=c(rep(1,8),rep(0,17), rep(1,4),  rep(0,24)),#index=0 for Z>3 ions, 1 otherwise. But now no longer needed.
  L = c(rep(1.6,8), rep(193, 7), rep(250, 4), rep(195, 6), rep(0.4, 4), rep(25, 5), rep(464, 4), rep(193, 3),rep(70, 4), rep(100, 5), rep(953, 3)), ##L=LET=LET_infinity=stopping power
  Z = c(rep(2, 8), rep(26, 17), rep(1, 4), rep(10, 5), rep(43, 4), rep(26, 3), rep(14, 4), rep(22, 5), rep(57, 3)),#proton #, e.g. 2 for 2He4
  Zeff = c(rep("TBD", 53)),# effective ion charge according to the formula of W.H Barkas. Zeff <= Z. Calculated below
  beta = c(rep("TBD", 53)),# ion speed, relative to speed of light, calculated below
  MeVperu = c(rep(228, 8), rep(600, 7), rep(300, 4), rep(600, 6), rep(250, 4), rep(670, 5), rep(600, 4), rep(600, 3), rep(260, 4), rep(1000, 5), rep(593, 3)),#Kinetic energy in MeV, divided by atomic mass, e.g. divided by 4u=4x931.5 MeV/c^2 for 2He4
  Katz = c(rep("TBD", 53)),#for fully ionized nuclei, Katz's Z^2/beta^2, Calculated below. It is part of the Bethe Barkas Bloch equation for stopping power. Our calculations don't use Katz, but various similar calculations do.
  ion = c(rep("He4", 8), rep("Fe56", 17), rep("p", 4), rep("Ne20", 5), rep("Nb93", 4), rep("Fe56", 3), rep("Si28", 4), rep("Ti48", 5), rep("La139", 3)),
  comments = c(".93AlpLooksOK", rep("", 7), ".93AlplooksOK", rep('', 11), ".93Alp.no.iso", "not in 17Cuc (or 16Chang?)", rep("", 3), "16Chang all OK?", rep('', 24), ".94Alp","From graphs",'e.g. in 17Cuc')) 

Y <- 0.001 * dfr[, "MeVperu"] # convert to GeV/u for convenience in a calculation
dfr[, "Katz"] <- round((dfr[, "Z"]) ^2 * (2.57 * Y^2 + 4.781 * Y + 2.233) / (2.57 * Y^2 + 4.781 * Y), 2) # special relativistic calculation of Z^2/beta^2. The numerics include conversion from GeV to joules and from u to kg.
dfr[, "beta"] <- round(dfr[, "Z"] * sqrt(1 / dfr[, "Katz"]), 3) # i.e. Z*sqrt(beta^2/Z^2) 
dfr[, "Zeff"] <- round(dfr[, "Z"] * (1 - exp( -125 * dfr[, "Z"] ^ (-2.0 / 3))), 2) # Barkas formula for Zeff; for us Zeff is almost Z

### Data for HG induced by photons from Cs-137 or Co-60 beta decay; from 16Chang
ddd <- data.frame(Dose = c(0, 0.4, 0.8, 1.6, 3.2, 7, 0, .4, .8, .12, 1.6),
                  HG = c(.026, .048, .093, .137, .322, .462, .0497, .054, .067, .128, .202),
                  NWeight = c(6081.2, 4989.5, 1896.8, 981.1, 522.2, 205.2, 7474.1, 2877.6, 1423.7, 689.9, 514.9),
                  Nucleus = c(rep("Cobalt-60", 6), rep("Cesium-137", 5)),
                  Comments = c(rep("TBD", 11))
                  )

# LL <- ddd[,"Dose"]
# HG <- ddd[,"HG"]
# WT <- ddd[,"NWeight"]
# QQ <- LL^2
# LQ <- lm(HG ~ LL + QQ)#,weights=WT)
LQ <- lm(HG ~ Dose + I(Dose ^ 2), data = ddd)  # EDWARD: I summarized lines 54-58 into line 59. Remove 54-58 if desired.
summary(LQ, correlation = T) 


### temporary chunk to enable visual tests of dfr
aa <- 20; bb <- 25# checking dfr for individual ions against graphs in .93Alp, 16Chang, and the new paper 17Cuc
plot(dfr[aa:bb, "Dose"], dfr[aa:bb, "HG"], ann = 'F') #example for checking dfr
####

dfra <- dfr[c(1:19, 26:53), ]  # removes the zero dose case and the no isograft data
############  HERE OUR NEW HZE MODEL; works well. Uses 3 adjustable parameters  ################ 
dfrHZE <- subset(dfra, Z > 3) # look only at HZE not at much lower LET ions. # In next line phi controls how fast NTE build up from zero; not really needed during calibration since 150*phi*Dose/L>>1 at every observed Dose !=0. phi needed for later synergy calculations
phi <- 3e3 # Even larger phi should give the same final results, but might cause extra problems with R.
HZEm <- nls(HG ~ .0275 + (1 - exp ( -0.01 * (aa1 * L * Dose * exp( -aa2 * L) + (1 - exp( -150 * phi * Dose / L)) * kk1))),
            data = dfrHZE, 
            weights = NWeight,
            start = list(aa1 = .9, aa2 = .01, kk1 = 0.048))  # calibrating parameters; 
summary(HZEm, correlation = T)  # ; vcov(HZEm)# parameter values & accuracy; variance-covariance matrix
HZEc <- coef(HZEm)  # calibrated central values of the 3 parameters.  Next is the IDER, =0 at dose 0
HZEC <- function(Dose, L) 1 - exp(-0.01 * (HZEc[1] * L * Dose * exp(-HZEc[2] * L) + (1 - exp(-150 * phi * Dose / L)) * HZEc[3]))#Calibrated IDER

########### Looks good up to here; next need Z < 3 model. Have good one; still looking if there is a better one
dfrL <- subset(dfra, Z <= 3)#for Light ions
L.m <- nls(HG ~ .0275 + bet * Dose^2 * exp(-lam * Dose),
           data = dfrL,
           weights = NWeight,
           start = list(bet = 5, lam = .2))
summary(L.m)

L.c <- coef(L.m)  # calibrated central values of the 2 parameters
L.C <- function(Dose, L) L.c[1] * Dose ^2 * exp(-L.c[2] * Dose)

# HHs=function(Dose,L,index){                 #Calculate slopes fpr later use in calculating I(d)
#   0.01*((1-index)*cc[1]*L*exp(-cc[2]*L)+index*cc[3]*(2*Dose-cc[4]*Dose^2)*exp(-cc[4]*Dose)+
#           exp(-150*phi*Dose/L)*(1-index)*cc[5]*150*phi/L)
#   }
# #Next is just a check that the derivative behaves the way it should.
# L=1.400;index=1; Dose=0.01*1:1000# or, e.g. L=100; index=0; Dose=0.01*1:300 for a heavy ion
# L=100; index=0; Dose=0.01*1:300
# HHv=HH(Dose,L,index)
# HHv2=c(0,HHv[1:(length(Dose)-1)])
# DDer=(HHv-HHv2)
# DDera=HHs(Dose-.005,L,index)
# tail(DDer); tail(DDera)
# ###### end derivative check. Seems to be fine; we can take this whole bit out pretty soon. The next bit is from Dae.

#We will adapt his use of ode() and uniroot() to our model and our data but all the rest of from Dae is irrelevant to us.
# MIXIDER_function = function(r, L, Z.beta, d = seq(0, 0.2, by = 0.001), eta0 = 1.300771e-04, eta1 = 3.164156e-03, sig0 = 2.481817e+00, kap = 2.565276e+02) {
#   dE=function(yini,State,Pars){
#     eta0 = eta0; eta1 = eta1; sig0 = sig0; kap = kap
#     with(as.list(c(State, Pars)), {
#       P = vector(length = length(L))
#       sig = vector(length = length(L))
#       etaa = vector(length = length(L))
#       u = vector(length = length(L))
#       for (i in 1:length(L)) {
#         P[i] = (1-exp(-Z.beta[i]/kap))^2
#         sig[i] = sig0*P[i] + 0.041/6.24*L[i]*P[i]
#         etaa[i] = eta0*L[i]*exp(-eta1*L[i])
#         u[i] = uniroot(function(d) sig[i]*6.24*d/L[i]*(1-exp(-1024*d/L[i])) + etaa[i]*(1-exp(-10^5*d)) - I, lower = 0, upper = 10, tol = 10^-10)$root
#       }
#       dI = vector(length = length(L))
#       for (i in 1:length(L)) {
#         dI[i] = r[i]*(sig[i]*6.24/L[i]*exp(-1024*u[i]/L[i])*(exp(1024*u[i]/L[i]) + 1024*u[i]/L[i] - 1) + etaa[i]*10^5*exp(-10^5*u[i]))
#       }
#       dI = sum(dI)
#       return(list(c(dI)))
#     })
#   }
#   pars = NULL; yini = c(I= 0); d = d
#   out = ode(yini,times = d, dE, pars, method = "radau")
#   return(out)
#}##### End from Dae 
#######Next: visual checks to see if our calibration is consistent with 16Chang, .93Alp, .94Alp and 17Cuc
## Put various values in our calibrated model to check with numbers and graphs in these references
L=1.6; Dose=dfrL[1:8,"Dose"];HGe=dfrL[1:8,"HG"]#He in .93Alp. HGe=experimental HG.
L=.4; Dose=dfrL[9:12,"Dose"];HGe=dfrL[9:12,"HG"]#same for protons. 
#L=193;Dose=dfrHZE[1:7,"Dose"]; HGe=dfr[1:7,"HG"] #same for Fe
plot(c(0,7),c(0,1),col='red', ann='F')# 
lines(Dose,L.C(Dose,L)+.0275)
points(Dose,HGe)#looks great for Helium, OK for protons; very good for iron. Mark: run some checks like these

#Call:
# lm(formula = HG ~ LL + QQ)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.04129 -0.02227 -0.01082  0.02043  0.08009 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  0.036825   0.019012   1.937  0.08877 . 
# LL           0.092911   0.021702   4.281  0.00268 **
#   QQ          -0.004479   0.003085  -1.452  0.18462   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04026 on 8 degrees of freedom
# Multiple R-squared:  0.9299,	Adjusted R-squared:  0.9124 
# F-statistic:  53.1 on 2 and 8 DF,  p-value: 2.408e-05
# 
# Correlation of Coefficients:
#   (Intercept) LL   
# LL -0.72            
# QQ  0.61       -0.96