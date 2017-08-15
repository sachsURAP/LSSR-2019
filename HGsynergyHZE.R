#This is free, open-source software under GNU GPLv3. It comes with no warranty. It concerns radiogenic mouse HG tumorigenesis
#Written by Mark Ebert, Edward Huang, Dae Woong, and Ray Sachs, summer 2017.
#Relevant references and their abbreviations in commenting this script are the following.
#".93Alp"=Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-particle radiations." Rad Res 136:382-391.(1993)
#".94Alp"=Alpen et al. "Fluence-based relative biological effectiveness for charged particle carcinogenesis in mouse Harderian gland." Adv Space Res 14(10): 573-581. (1994).  
#"16Chang"=Chang et al. "Harderian Gland Tumorigenesis: Low-Dose and LET Response." Radiat Res 185(5): 449-460. (2016).  
#"16Srn"=1.	Siranart et al."Mixed Beam Murine Harderian Gland Tumorigenesis: Predicted Dose-Effect Relationships if neither Synergism nor Antagonism Occurs." Radiat Res 186(6): 577-591. (2016).  
#"17Cuc"=Cucinotta & Cacao. "Non-Targeted Effects Models Predict Significantly Higher Mars Mission Cancer Risk than Targeted Effects Models." Sci Rep 7(1): 1832. (2017).  PMC5431989
rm(list=ls())
library(deSolve)
dfr=data.frame( #data used in 16Chang  
  Dose=c(0.2,0.4,0.6,1.2,2.4,3.2,5.1,7,0.05,0.1,0.15,0.2,0.4,0.8,1.6,0.05,0.1,0.2,0.4,0,0.1,0.2,0.4,0.8,1.6,0.4,0.8,1.6,3.2,0.05,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.04,0.08,0.16,0.32,0.033,0.066,0.13,0.26,0.52,.2, .4, .6),
  HG=c(0.091,0.045,0.101,0.169,0.347,0.431,0.667,0.623,0.156,0.215,0.232,0.307,0.325,0.554,0.649,0.123,0.145,0.207,0.31,0.026,0.083,0.25,0.39,0.438,0.424,0.093,0.195,0.302,0.292,0.109,0.054,0.066,0.128,0.286,0.183,0.167,0.396,0.536,0.192,0.234,0.317,0.092,0.131,0.124,0.297,0.082,0.088,0.146,0.236,0.371,.154,.132,.333),#HG prevalence as defined in 16Chang
  NWeight=c(520,2048,1145,584,313,232,293,221,1162,877,455,409,374,223,320,742,661,347,131,6081,1091,251,244,191,131,645,255,199,111,649,378,973,833,201,468,381,197,109,496,257,185,1902,1063,884,350,1767,1408,874,299,261,322,206,67),
  #nominal weight for weighted least squaresregression; see .93Alp. The Lanthanum entries were obtained by measuring the main graph in 17Cuc 
  #index=c(rep(1,8),rep(0,17), rep(1,4),  rep(0,24)),#index=0 for Z>3 ions, 1 otherwise. But now no longer needed.
  L=c(rep(1.6,8),rep(193,7),rep(250,4),rep(195,6),rep(0.4,4),rep(25,5),rep(464,4),rep(193,3),rep(70,4),rep(100,5),rep(953,3)), ##L=LET=LET_infinity=stopping power
  Z=c(rep(2,8),rep(26,17),rep(1,4),rep(10,5),rep(43,4),rep(26,3),rep(14,4),rep(22,5),rep(57,3)),#proton #, e.g. 2 for 2He4
  Zeff=c(rep("TBD",53)),# effective ion charge according to the formula of W.H Barkas. Zeff <= Z. Calculated below
  beta=c(rep("TBD",53)),# ion speed, relative to speed of light, calculated below
  MeVperu=c(rep(228,8),rep(600,7),rep(300,4),rep(600,6),rep(250,4),rep(670,5),rep(600,4),rep(600,3),rep(260,4),rep(1000,5),rep(593,3)),#Kinetic energy in MeV, divided by atomic mass, e.g. divided by 4u=4x931.5 MeV/c^2 for 2He4
  Katz=c(rep("TBD",53)),#for fully ionized nuclei, Katz's Z^2/beta^2, Calculated below. It is part of the Bethe Barkas Bloch equation for stopping power. Our calculations don't use Katz, but various similar calculations do.
  ion=c(rep("He4",8),rep("Fe56",17),rep("p",4),rep("Ne20",5),rep("Nb93",4),rep("Fe56",3),rep("Si28",4),rep("Ti48",5),rep("La139",3)),
  comments=c(".93AlpLooksOK",rep("",7),".93AlplooksOK",rep('',11),".93Alp.no.iso", "not in 17Cuc (or 16Chang?)",rep("",3),"16Chang all OK?",rep('',24),".94Alp","From graphs",'e.g. in 17Cuc')) 
Y=0.001*dfr[,"MeVperu"]# convert to GeV/u for convenience in a calculation
dfr[,"Katz"]=round((dfr[,"Z"])^2*(2.57*Y^2+4.781*Y+2.233)/(2.57*Y^2+4.781*Y),2)#special relativistic calculation of Z^2/beta^2. The numerics include conversion from GeV to joules and from u to kg.
dfr[,"beta"]=round(dfr[,"Z"]*sqrt(1/dfr[,"Katz"]),3)#i.e. Z*sqrt(beta^2/Z^2)
dfr[,"Zeff"]=round(dfr[,"Z"]*(1-exp(-125*dfr[,"Z"]^(-2.0/3))),2)#Barkas formula for Zeff; for us Zeff is almost Z


dfra=dfr[c(1:19,26:53),] ##removes the zero dose case and the no isograft data

##HERE OUR NEW HZE/NTE MODEL; works well. Uses 3 adjustable parameters. There is also an HZE/TE model and NTE or TE models for Z<=3, but this script only has the HZE/NTE model: all models are for Z>3 and all data is for Z>=8  ## 
dfrHZE=subset(dfra,Z>3) #look only at HZE not at much lower LET ions. #In next line phi controls how fast NTE build up from zero; not really needed during calibration since 150*phi*Dose/L>>1 at every observed Dose !=0. phi needed for later synergy calculations
phi=1000# Even larger phi should give the same final results, but might cause extra problems with R.
HZEm=nls(HG~.0275+(1-exp(-0.01*(aa1*L*Dose*exp(-aa2*L)+(1-exp(-150*phi*Dose/L))*kk1))),data=dfrHZE, weights=NWeight,start=list(aa1=.9,aa2=.003, kk1=6))#calibrating parameters; 
summary(HZEm,correlation=T)#; vcov(HZEm)# parameter values & accuracy; variance-covariance matrix
HZEc=coef(HZEm)#calibrated central values of the 3 parameters.  Next is the IDER, =0 at dose 0
HHC=function(Dose,L) 0.01*(HZEc[1]*L*Dose*exp(-HZEc[2]*L)+(1-exp(-150*phi*Dose/L))*HZEc[3])#calibrated hazard function
HZEC=function(Dose,L) 1-exp(-HHC(Dose,L))#Calibrated IDER

HZES=function(Dose,L){                 #Calculate slopes for later use in calculating I(d)
  0.01*(HZEc[1]*L*exp(-HZEc[2]*L)+exp(-150*phi*Dose/L)*HZEc[3]*150*phi/L)*exp(-HHC(Dose,L))
  }# made a numerical check at dose 10 for L=193. Got agreement to six significant figures
#Looks OK up to here
dose=c(seq(0, .00001, by = 0.000001),seq(.00002,.0001,by=.00001),seq(.0002,.001,by=.0001),seq(.002,.01,by=.001),seq(.02,.5,by=.01))##look carefully near zero, but go out to 0.5 Gy
#dose=dose[1:30] #this is used to zoom in on the very low dose behavior in the graphs
#We will adapt Dae's use of ode() and uniroot() to our model and our data but all the rest of from Dae is irrelevant to us.
MIXDER_function = function(r, L,d =dose,aa1=HZEc[1],aa2=HZEc[2],kk1=HZEc[3]) {
   dE=function(yini,State,Pars){
     aa1=aa1;aa2=aa2;kk1=kk1
     with(as.list(c(State, Pars)), {
       aa= vector(length = length(L))
       u = vector(length = length(L))
       for (i in 1:length(L)) {
         aa[i] = aa1*L[i]*exp(-aa2*L[i])
         u[i] = uniroot(function(d) 1-exp(-0.01*(aa[i]*d+(1-exp(-150*phi*d/L[i]))*kk1)) - I, lower = 0, upper = 10, tol = 10^-10)$root
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
   out = ode(yini,times = d, dE, pars, method = "radau")
   return(out)
}##### End from Dae 
r=c(1/3,1/3,1/3)
L=c(25,70,250)
INCRL=MIXDER_function(r,L)#incremental effect additivity
plot(INCRL, type='l',ylim=c(0,.4),col='red',bty='l',ann='F')
lines(dose,HZEC(dose,250),col='green')#component 3
lines(dose,HZEC(dose,70),col='green')#component 2
lines(dose,HZEC(dose,25),col='green') #component 1
SEA=function(Dose) HZEC(Dose/3,25)+HZEC(Dose/3,70)+HZEC(Dose/3,250)
lines(dose,SEA(dose),lty=2)
r=.25*1:4
L=c(25,70,190,250)
plot(INCRL, type='l',ylim=c(0,.5),col='red',bty='l',ann='F')
lines(dose,HZEC(dose,190), col='green')#component 3
lines(dose,HZEC(dose,250),col='green')#component 4
lines(dose,HZEC(dose,70),col='green')#component 2
lines(dose,HZEC(dose,25),col='green') #component 1
SEA=function(Dose) HZEC(Dose/4,25)+HZEC(Dose/4,70)+HZEC(dose/4,190)+HZEC(Dose/3,250)
lines(dose,SEA(dose),lty=2)
#To Be Done: 3. graphs at dose ranges (0,.5) and (0,0.0001)
##Then: add this script to the main script (there is a lot of overlap)
##Then I(d) for mixtures involving light ions
##Then Monte Carlofor 95%CI