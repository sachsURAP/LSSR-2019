##This is free, open-source software under GNU GPLv3. It comes with no warranty. Concerns radiogenic mouse HG tumorigenesis
#Written by Mark Ebert and Ray Sachs, summer 2017.
rm(list=ls())
ddfr=data.frame( #data used in Chang16=Chang et al."Harderian Gland..." Rad Res 185(5): 449-460 (2016)  
  Dose=c(0.2,0.4,0.6,1.2,2.4,3.2,5.1,7,0.05,0.1,0.15,0.2,0.4,0.8,1.6,0.05,0.1,0.2,0.4,0,0.1,0.2,0.4,0.8,1.6,0.4,0.8,1.6,3.2,0.05,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.04,0.08,0.16,0.32,0.033,0.066,0.13,0.26,0.52,.2, .4, .6),
  HG=c(0.091,0.045,0.101,0.169,0.347,0.431,0.667,0.623,0.156,0.215,0.232,0.307,0.325,0.554,0.649,0.123,0.145,0.207,0.31,0.026,0.083,0.25,0.39,0.438,0.424,0.093,0.195,0.302,0.292,0.109,0.054,0.066,0.128,0.286,0.183,0.167,0.396,0.536,0.192,0.234,0.317,0.092,0.131,0.124,0.297,0.082,0.088,0.146,0.236,0.371,.154,.132,.333),#HG prevalence as defined in 16Chang
  NWeight=c(520,2048,1145,584,313,232,293,221,1162,877,455,409,374,223,320,742,661,347,131,6081,1091,251,244,191,131,645,255,199,111,649,378,973,833,201,468,381,197,109,496,257,185,1902,1063,884,350,1767,1408,874,299,261,322,206,67),#nominal weight for weighted least squaresregression; see [Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-particle radiations." Rad Res 136(3):382-391.(1993)], "Alp" for short and the earlier Lanthanum data
  index=c(rep(1,8),rep(0,17), rep(1,4),  rep(0,24)),#index=0 for Z>3 ions, index=1 for proton p (=1H1) and 2HE4 ions
  L=c(rep(1.6,8),rep(193,7),rep(250,4),rep(195,6),rep(0.4,4),rep(25,5),rep(464,4),rep(193,3),rep(70,4),rep(100,5),rep(953,3)), ##L=LET=LET_infinity=stopping power
  Z=c(rep(2,8),rep(26,17),rep(1,4),rep(10,5),rep(43,4),rep(26,3),rep(14,4),rep(22,5),rep(57,3)),#proton #, e.g. 2 for 2He4
  Zeff=c(rep("TBD",53)),# effective ion charge according to the formula of W.H Barkas. Zeff <= Z. Calculated below
  beta=c(rep("TBD",53)),# ion speed, relative to speed of light, calculated below
  MeVperu=c(rep(228,8),rep(600,7),rep(300,4),rep(600,6),rep(250,4),rep(670,5),rep(600,4),rep(600,3),rep(260,4),rep(1000,5),rep(593,3)),#Kinetic energy in MeV, divided by atomic mass, e.g. divided by 4u=4x931.5 MeV/c^2 for 2He4
  Katz=c(rep("TBD",53)),#for fully ionized nuclei, Katz's Z^2/beta^2, Calculated below. It is part of the Bethe Barkas Bloch equation for stopping power. Our calculations don't use Katz, but many similar calculations do.
  ion=c(rep("He4",8),rep("Fe56",17),rep("p",4),rep("Ne20",5),rep("Nb93",4),rep("Fe56",3),rep("Si28",4),rep("Ti48",5),rep("La139",3)),
  comments=c("AlpLooksOK",rep("",7),"AlplooksOK",rep('',11),"Alp.no.iso", "not in Cuc (or Chang?)",rep("",3),"Chang all OK?",rep('',24),"Alp","",'')) 
Y=0.001*ddfr[,"MeVperu"]# convert to GeV/u for convenience in a calculation
ddfr[,"Katz"]=round((ddfr[,"Z"])^2*(2.57*Y^2+4.781*Y+2.233)/(2.57*Y^2+4.781*Y),2)#special relativistic calculation of Z^2/beta^2. The numerics include conversion from GeV to joules and from u to kg.
ddfr[,"beta"]=round(ddfr[,"Z"]*sqrt(1/ddfr[,"Katz"]),3)#i.e. Z*sqrt(beta^2/Z^2)
ddfr[,"Zeff"]=round(ddfr[,"Z"]*(1-exp(-125*ddfr[,"Z"]^(-2.0/3))),2)#Barkas formula for Zeff; for us Zeff is almost Z

### temporary chunk to enable visual tests of ddfr
aa=20;bb=25# checking ddfr for individual ions against graphs in 93Alp, 16Chang, and the new paper 17CucNTEhiRisk
plot(ddfr[aa:bb,"Dose"],ddfr[aa:bb,"HG"], ann='F') #example for checking ddfr
####

ddfra=ddfr[c(1:19,26:53),] ##removes the zero dose case and the no isograft data

######Mark: Here our model; works very well. could omit (1-exp(-150*phi*Dose/L) (=1 at every observed dose point !=0)
phi=3e3#controls how fast NTE build up from zero; not really needed since 150*phi*Dose/L>>1 at every observed Dose !=0  
IDERm=nls(HG~1-exp(-0.01*(2.75+(1-index)*aa1*L*Dose*exp(-aa2*L) +index*bb*Dose^2*exp(-ll0*Dose)+
                            (1-exp(-150*phi*Dose/L))*(1-index)*kk1)),data=ddfra, weights=NWeight,
  start=list(aa1=.9,aa2=.01,bb=4.5,ll0=.2, kk1=0.048))#calibrating parameters; only need 5 parameters for very good results
summary(IDERm,correlation=T); vcov(IDERm)# parameter values & accuracy; variance-covariance matrix
cc=coef(IDERm)#calibrated central values of the 5 parameters
HH=function(Dose,L,index){
  0.01*(2.75+(1-index)*cc[1]*L*Dose*exp(-cc[2]*L) +index*cc[3]*Dose^2*exp(-cc[4]*Dose)+
           (1-exp(-150*phi*Dose/L))*(1-index)*cc[5])
}
HGc=function(Dose,L,index) 1-exp(-HH(Dose,L,index))##this is the calibrated model. 

HHs=function(Dose,L,index){                 #Calculate slopes fpr later use in calculating I(d)
  0.01*((1-index)*cc[1]*L*exp(-cc[2]*L)+index*cc[3]*(2*Dose-cc[4]*Dose^2)*exp(-cc[4]*Dose)+
          exp(-150*phi*Dose/L)*(1-index)*cc[5]*150*phi/L)
  }
#Next is just a check that the derivative behaves the way it should.
L=1.400;index=1; Dose=0.01*1:1000# or, e.g. L=100; index=0; Dose=0.01*1:300 for a heavy ion
L=100; index=0; Dose=0.01*1:300
HHv=HH(Dose,L,index)
HHv2=c(0,HHv[1:(length(Dose)-1)])
DDer=(HHv-HHv2)
DDera=HHs(Dose-.005,L,index)
tail(DDer); tail(DDera)
###### end derivative check. Seems to be fine; we can take this whole bit out pretty soon. The next bit is from Dae.
 
#We will adapt his use of ode() and uniroot() to our model and our data but all the rest of from Dae is irrelevant to us.
MIXIDER_function = function(r, L, Z.beta, d = seq(0, 0.2, by = 0.001), eta0 = 1.300771e-04, eta1 = 3.164156e-03, sig0 = 2.481817e+00, kap = 2.565276e+02) {
  dE=function(yini,State,Pars){
    eta0 = eta0; eta1 = eta1; sig0 = sig0; kap = kap
    with(as.list(c(State, Pars)), {
      P = vector(length = length(L))
      sig = vector(length = length(L))
      etaa = vector(length = length(L))
      u = vector(length = length(L))
      for (i in 1:length(L)) {
        P[i] = (1-exp(-Z.beta[i]/kap))^2
        sig[i] = sig0*P[i] + 0.041/6.24*L[i]*P[i]
        etaa[i] = eta0*L[i]*exp(-eta1*L[i])
        u[i] = uniroot(function(d) sig[i]*6.24*d/L[i]*(1-exp(-1024*d/L[i])) + etaa[i]*(1-exp(-10^5*d)) - I, lower = 0, upper = 10, tol = 10^-10)$root
      }
      dI = vector(length = length(L))
      for (i in 1:length(L)) {
        dI[i] = r[i]*(sig[i]*6.24/L[i]*exp(-1024*u[i]/L[i])*(exp(1024*u[i]/L[i]) + 1024*u[i]/L[i] - 1) + etaa[i]*10^5*exp(-10^5*u[i]))
      }
      dI = sum(dI)
      return(list(c(dI)))
    })
  }
  pars = NULL; yini = c(I= 0); d = d
  out = ode(yini,times = d, dE, pars, method = "radau")
  return(out)
}##### End from Dae 

#######Next: visual checks to see if our calibration is consistent with 16Chang, 93Alp, 17CucCacao.
## Put various values in our calibrated simplified model to check with numbers and graphs in these references
L=1.6; Dose=ddfr[1:8,"Dose"];HGe=ddfr[1:8,"HG"];index=1#He in Alpen93. HGe=experimental HG.
#L=.4; Dose=ddfr[26:29,"Dose"];HGe=ddfr[26:29,"HG"];index=1#same for protons. 
L=193;Dose=ddfr[9:15,"Dose"]; HGe=ddfr[9:15,"HG"];index=0#same for Fe
plot(c(0,7),c(0,1),col='red', ann='F')# biggest Dose-effect range ever considered
lines(Dose,HGc(Dose,L,index))
points(Dose,HGe)#looks great for Helium, OK for protons; very good for iron. Mark: run some checks like these
