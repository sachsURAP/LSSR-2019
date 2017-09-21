##This is free, open-source software under GNU GPLv3. It comes with no warranty. Concerns radiogenic mouse HG gland tumorigenesis. Uses individual 
##dose-response relations (IDER) & parameters in Chang et  al. "Harderian Gland Tumorigenesis: Low-Dose and LET Response" Radiat Res 185 2016  
##User adjustable parameters are those starting near line 147 and also limits on x and y axes for some of the graphs
##Gives main calculations for Siranart, Blakely, Cheng, Handa, Sachs, manuscript submitted to Rad. Res.
#Written by Nopphon Siranart; tested and modified by Cheng, Handa,  and Sachs
library(ggplot2); library(deSolve)
rm(list=ls()) 
###FIXED PARAMETERS###  For individual ion dose-response relations
#next lines have Cucinotta et al. (2013) table 5.4 IDER parameters for ions. TE=targeted effects, NTE = non-TE formerly known as bystander. 
# parameter = list(TE = list(alph0=7.53, alph1=1.261, alph2=0.0037, beta = 0, beta_p = 6.3, lamd0=.25, lamd1=.0051, lamd2=.0034,
#                       e_alph0=3.96, e_alph1=0.213, e_alph2=0.00058, e_beta = 0, e_beta_p = 3.41, e_lamd0=.065, e_lamd1=.0029, e_lamd2=.0027), 
#                  NTE=list(alph0=10.02,alph1=0.679,alph2=0.0033,beta=0,beta_p=5.08,lamd0=.231,lamd1=.0033,lamd2=.005,k1=0.12,k2=0.0053,
#                       e_alph0=2.07,e_alph1=0.187,e_alph2=0.0006,e_beta=0,e_beta_p=3,e_lamd0=.016,e_lamd1=.0042,e_lamd2=.0064,e_k1=0.06,e_k2=0.002))

#next lines: values in Chang et al. 2016. These override the parameters above. 
parameter = list(TE = list(alph0=7.65, alph1=1.25, alph2=0.0038, beta = 0, beta_p = 6.02, lamd0=.243, lamd1=.006, lamd2=.0043,
                           e_alph0=3.94, e_alph1=0.14, e_alph2=0.0004, e_beta = 0, e_beta_p =3.51, e_lamd0=.07, e_lamd1=.0036, e_lamd2=.0027),
                 NTE=list(alph0=10.05,alph1=0.90,alph2=0.0039,beta=0,beta_p=4.61,lamd0=.219,lamd1=.0047,lamd2=.0051,k1=0.048,k2=0.0028,#for NTE1
                          e_alph0=3.56,e_alph1=0.21,e_alph2=0.0009,e_beta=0,e_beta_p=3.33,e_lamd0=.078,e_lamd1=.0059,e_lamd2=.0059,e_k1=0.023,e_k2=0.0019))

alph=function(L, param) {with(param,alph0+alph1*L*exp(-alph2*L))} #slope near zero of the IDER; depends on linear energy transfer LET_inf=L (keV/micron)
lamd=function(L, param) {with(param,lamd0+lamd1*L*exp(-lamd2*L))}
#lamd is controversially interpreted as giving a cell killing modulation of LQ IDER. It leads to IDER which have maxima at doses > 1 Gy.
kappa=function(L,param){with(param, (k1*L)*exp(-k2*L))}#kappa interpreted as NTE, important near d=0; approximated as initial value.

###R FUNCTIONS### Used later (near line 180) after "USER_DEFINED PARAMETERS" have been set 
#Next R function can calculate an IDER "targeted part" TE (Cucinotta et al. 2013); uses proton=0 for ion charge Z>2; proton=1 for H, He ions
T_single = function(dose, L, Eff, proton, param = parameter[[Eff]]){ #if Eff=1=TE=targeted effect model; Eff=2=NTE=non-targeted effect.
  beta = param$beta #initialize to be beta=0 with Z>2, i.e near d=0 TE effect is linear no threshhold (LNT) for Z>2.
  if(proton == 1){ beta= param$beta_p} #ions with Z<3 are approximated as LQ (linear quadratic) near d=0.
  P_E = (alph(L,param)*dose + beta*(dose^2))*exp(-lamd(L,param)*dose) #LQ or LNT modulated by exponentially decreasing factor.
  return(P_E)
}

E_single = function(dose, L, Eff, proton, param = parameter[[Eff]]){#if Eff=1, E_single=T; #if Eff=2, E_single=T+NT, where NT=kappa*exp(-lamd*dose).  
  P_E = T_single(dose,L,Eff,proton, param)
  if (Eff==2){#phi is adjustable parameter for mGy range NTE1; adjust near line 170
    P_E = P_E +(1-exp(-150*phi*dose/L))*kappa(L,param)*exp(-lamd(L,param)*dose)#P_E here is sum of targeted and non-targeted effects
  }
  return(P_E)
} 

#Each E_single has a maximum, often around dose=1 or 2 Gy. Next is a function that can find the maxima and the doses where they occur.
bound_Ed = function(L,Eff,proton, param){
  E = c(); D = c();
  for(ll in 1:length(L)){#when this function is used, L can be a vector of the LETs of the individual components of a mixed-ion beam
    p = proton[ll]; l = L[ll]; 
    b = param[[ll]]$beta; a = alph(l,param[[ll]]); ld = lamd(l,param[[ll]])
    if(p == 1){ b= param[[ll]]$beta_p}
    f = function(x) E_single(x,l,Eff,p,param[[ll]])
    d = optimize(f, c(0, 10), maximum= TRUE,tol=1e-10)[[1]]
    D = c(D,d)
    E = c(E,solve_d(d,l,Eff,p,0, param[[ll]]))
  }
  upper_E = min(E)#smallest of the maxima in a mixed beam
  return(list(upper_E,D))
}

#R function whose zero uniroot can use to get dose from L,Eff,proton, and E_single by treating, as in RBE calculations, 
#E=E_single as independent and d as dependent variable instead of vice versa .
solve_d=function(dose,L,Eff,proton,E, param = parameter[[Eff]]) {
  eq = E_single(dose,L,Eff,proton,param)-E
  return(eq)
}

gamma_param = function(Eff){  #When calculating 95% CI for mixture, IDER parameters are taken to have gamma distributions (thus staying>0)
  param = parameter[[Eff]]    #The gamma distributions have mean M, variance V, shape M^2/V, rate (i.e. 1/scale)=M/V 
  param_error=with(param,list(alph0=rgamma(1,(alph0/e_alph0)^2,alph0/(e_alph0)^2),alph1=rgamma(1,(alph1/e_alph1)^2,alph1/(e_alph1)^2),
   alph2=rgamma(1,(alph2/e_alph2)^2,alph2/(e_alph2)^2),beta=0,beta_p=rgamma(1,(beta_p/e_beta_p)^2,beta_p/(e_beta_p)^2),
   lamd0=rgamma(1,(lamd0/e_lamd0)^2,lamd0/(e_lamd0)^2),lamd1=rgamma(1,(lamd1/e_lamd1)^2,lamd1/(e_lamd1)^2),
   lamd2=rgamma(1,(lamd2/e_lamd2)^2,lamd2/(e_lamd2)^2)))
  if(Eff ==2){
    param_error1 = with(param,list(k1 = rgamma(1,(k1/e_k1)^2,k1/(e_k1)^2), k2 = rgamma(1,(k2/e_k2)^2,k2/(e_k2)^2)))
    param_error = c(param_error,param_error1)
  }
  return(param_error)
}

#R function that generates nL=(no. of mixture components) identical parameter sets during 1 Monte Carlo run in estimates of 95% Confidence Intervals (CI)
gen_param = function(Eff,NL, random){  #NL=nL= number of individual ions in a mixed beam
  param = list()
  p = gamma_param(Eff) #one set of parameters chosen at random
  for(ii in 1:NL){
    param[[ii]] = p #every mixture component gets this set. In next Monte Carlo run, gen_param is used again and makes a new set.
    if (random ==0) param[[ii]] = parameter[[Eff]]
  }
  return(param)
}

#R function for d = D_j(E), where D_j will be the compositional inverse function to the IDER for the jth ion in a mixture (corresponding to jth component of L). 
dose_E=function(E, L, Eff,proton, param){#sometimes used with E a vector of effects, not just a single effect, output a vector of doses
  result0 = bound_Ed(L,Eff,proton, param); upper_E = result0[[1]]; d_max = result0[[2]] #restrict to doses & effects such that all IDER monotonic increasing
  output=c()
  e =min(E, upper_E)
  dose = (uniroot(solve_d,c(0,.001), upper=d_max, extendInt="upX",tol = 1e-6, L = L, 
                    Eff = Eff, proton=proton, E = e, param = param[[1]]))$root
  output= c(output,dose)
  return(output)
}

#Next is R function that can get simple effect additivity default prediction E_A=S(d)
E_A <- function(doses, L, Eff, proton, param = parameter[[Eff]]){
  output = rep(0,cols) #cols is the maximum dose in cGy considered
  for(l in 1:nL){
    output = output + E_single(r[l]*doses,L[l],Eff,proton[l],param[[l]])
  }
  return(output)
}

#R function to get first derivative of a function f evaluated at E, using effect rather than dose as independent variable
first_deriv = function(f,E){
  
  delta = 1e-7
  output = (f(E) - f(E-delta))/delta
  return(output)
}

sum_rE = function(E, param){  #can be used to calculate mixture slope for integrating E_I=I(d); add IDER slopes for given mixture effect E=E_I vector
  output = c()
  for(j in 1:length(E)){
    cumrE = 0
    for(i in 1:nL){
      f = function(x) E_single(x,L[i],Eff,proton[i],param[[i]])
      temp = list(); temp[[1]] = param[[i]]   #note param[[i]], not param[[Eff]]; param[ii] changes from run to run, remains same within one run
      d = dose_E(E[j],L[i],Eff,proton[i],param=temp)
      cumrE = cumrE + r[i]*first_deriv(f,d)
    }
    output = c(output,cumrE)
  }
  return(output)
}

solve_ode <- function(t, state, parameters) {#R function used when integrating the ordinary differential eq. for E_I=I(d)
  with(as.list(c(state, parameters)), {
    dE <- sum_rE(E, parameters)
    #print(dE)
    list(c(dE))
  })
}

E_I <- function(state,doses,L, Eff, proton, param){#R function used to actually integrate
  output = ode(y = state, times = doses, func = solve_ode, parms = param)
  result = output[,2]
   return(result)
}

#####EXAMPLES OF USER-DEFINED PARAMETERS FOR MIXTURE CALCULATIONS#####
#r=mixture proportions (must sum to 1); L=LETs; proton=1 for H, He, and proton =0 for heavier ions (0 means LQ beta is 0). 
#end=largest dose in Gy; N=number of Monte-Carlo runs when calculating 95%CI. Toggle between Eff=1 and Eff=2 as needed. Examples follow
#r=c(.5,.5);L=c(76,175);proton=c(0,0);Eff=2;N=50;end=.5 #2 HZE ions. Use bigger N  (e.g. N=2000) for more accurate CI
#r=c(.70,.30);L=c(170,190);proton=c(0,0);Eff=1;N=50; end=.5 
#r=c(.4,.4,.2); L=c(100,200,300); proton = c(0,0,0); Eff=2; N = 50;end=.5
#Next is a mixture of 10 different heavy ions.
#r=c(.02*c(9,8,7,6,5,5,4,3,2,1)); L=c(25*(1:10)); proton =rep(0,10); Eff=1; N = 50;end=.5 #use N >= 50 to get 95% CI
#r=c(.7,.006*c(9,8,7,6,5,5,4,3,2,1)); L=c(.4,25*(1:10)); proton =c(1,rep(0,10)); Eff=2; N = 50; end=1
#r=c(.7,rep(.003*c(9,8,7,6,5,5,4,3,2,1),2)); L=c(.4,50*(1:5),40*(1:5),30*(1:5),20*(1:5)); proton =c(1,rep(0,20)); Eff=2; N = 500; end=1
#r=1;L=.4; proton=1; Eff=2; N=50;end=.6 #a "mixture" with just one component; 
#r=c(.4,.4,.2); L=c(50,50,200); proton = c(0,0,0); Eff=2; N = 100; end=1  #test consistency with sham mixture principle using next line
#r=c(.8,.2); L=c(50,200); proton = c(0,0); Eff=2; N = 2000; end=1  #should be the same as the previous line
#r=c(.8,.12,.08); L=c(.4,76,174); proton = c(1,0,0); Eff=2; N = 50;end=1
#r=c(.76,.18,.06); L=c(.4,60,190); proton = c(1,0,0); Eff=1; N = 50;end=1
#r=c((1/50)*c(13,11,9,7,5,3,2)); L=c(60,80,35*(3:7)); proton=rep(0,7);Eff=2;N=50; end=1  # 7 HZE for supplement
r=c(.70,.30);L=c(60,190);proton=c(0,0);Eff=2;N=50;end=1# Next are examples for the MS (Table 4).  #N=50 instead of N>=1000 for speed.
#r=c(.6,.4); L=c(.4,76); proton=c(1,0);Eff=2; N=1000; end=1  # panels B & B*
#r=.01*c(50,20,rep(5,6));L=c(.4,1.4,25,21,76,107,174,464);proton=c(1,1,rep(0,6)); Eff=1; N=1000; end=1 #panels D&D* 8-ion expt.
#r=.01*c(28,20,20,12,12,4,4); L=c(60,80,35*(3:7)); proton=rep(0,7);Eff=1;N=1000;end=.5  #F and F* seven heavy ions. absurd S(d)
#r=c(.5,.5); L=c(76,174); proton = c(0,0); Eff=2; N =2000; end=.4 #E and E* 
#r=c(.9,.1);L=c(.4,174);proton=c(1,0);Eff=2;N=1000; end=1  #panels A and A*
#r=c(.60,.40);L=c(.4,174);proton=c(1,0);Eff=2;N=1000;end=.7  #panels C and C* 
#r=c(.6,.2,.2);L=c(.4,76,174);proton=c(1,0,0);Eff=2;N=1000;end=1
by = 0.05#to increase smoothness decrease "by"
monte =T#TRUE/FALSE. If FALSE no mixture 95% CI is calculated and the script runs much faster.
phi=3e3# adjusts how steeply non-targeted effect climbs from 0 to kappa in the mGy range.
#####END OF MAIN USER DEFINED PARAMETERS######

nL = length(L) # The number of ions used
doses <- seq(0, end, by = by);cols = end/by+1 #end=largest dose;cols=total number of dose points if Eff=1
if(Eff==2){#look closely at mGy region
  doses=c(3e-4*(0:100),seq(.031,end,length.out=100));cols=length(doses)
}

# ############E_A with central values of the parameters##############
E_A_mean = E_A(doses, L ,Eff,proton, gen_param(Eff,nL,0))

############E_I with central values of the parameters###############
state <- c(E = 0);
E_I_mean = E_I(state,doses, L, Eff, proton,gen_param(Eff,nL,0))

# ##########95% confidence intervals;  Most CPU intensive part of the calculation
if(monte){
  E_A_mat = matrix(0, N, cols); ##now do Monte Carlo for error structure
  for(i in 1:N){
    E_A_mat[i,] = E_A(doses, L, Eff,proton, gen_param(Eff,nL,1))
  }
  E_A_sort = apply(E_A_mat,2,sort,decreasing=F); 
  E_A_low = E_A_sort[floor(N*0.025),]; E_A_high = E_A_sort[ceiling(N*0.975),]; 
 
  E_I_mat = matrix(0, N, cols);
  for(i in 1:N){
    print(i)
    param = gen_param(Eff,nL,1)
    E_I_mat[i,] = E_I(state,doses,L,Eff,proton, param)
  }
  E_I_sort = apply(E_I_mat,2,sort,decreasing=F); 
  E_I_low = E_I_sort[floor(N*0.025),]; E_I_high = E_I_sort[ceiling(N*0.975),]; 
}
####################GRAPHS####################Comment out unneeded ones
##IDER Plot. Allows adjustment of upper dose limit to less than "end" under user adjsted parameters and of y limit.
# xstart = 0; xend = end; ystart= 0; yend =40#
# plot(c(xstart,xend),c(ystart,yend),type='n',bty="l", ylab = 'E', xlab = 'd')#,xaxp=c(0,end,2))
# lines(c(0,xend),c(yend,yend))     
# for(l in 1:length(L)){
#   lines(doses,E_single(doses,L[l],Eff,proton[l]),type="l", col = colors()[20*l+2],xaxp=c(0,end,4)) #automatic choice of colors 
# }
# 
# #####plot E_I, and E_A together in one graph for the central parameters. Add IDER and/or 95%CI if desired.
xstart = 0; xend =.7; ystart= 0; yend =30
plot(c(xstart,xend),c(ystart,yend))#,type='n',bty="l", ylab = 'E', xlab = 'd', xaxp=c(0,end,2), yaxp=c(0,yend,2))##,xaxs='i',yaxs='i'
abline(h=20,lwd=1,lty=3)
#abline(v=c(.50),col='green', lwd=1.5)
lines(doses,E_I_mean,col='red', lwd=3)
lines(doses,E_A_mean,col='blue', lwd=3)
for(l in 1:length(L)){
  lines(doses,E_single(doses,L[l],Eff,proton[l]),type="l",col="green") 
}
if(monte){
  lines(doses,E_I_low,col='red', lwd=2,lty=2)
  lines(doses,E_I_high,col='red', lwd=2,lty=2)
  lines(doses,E_A_low,col='blue', lwd=2,lty=2)
  lines(doses,E_A_high,col='blue', lwd=2,lty=2)
}
# # ######Elegant ribbon plots######
# ## plot E_I and E_A together in one graph with error bars, allowing the user to adjust the upper dose limit (<2 Gy; if no protons <1)
df= data.frame(md = doses, EI=E_I_mean, EA=E_A_mean)##dev.off()
if(monte){
  p = ggplot(df,aes(md)) + geom_line(aes(y = EI),color ='red') + geom_ribbon(aes(ymin=E_I_low, ymax=E_I_high), alpha=0.2,fill='red')+
  geom_line(aes(y = EA),color ='blue') + geom_ribbon(aes(ymin=E_A_low, ymax=E_A_high), alpha=0.2,fill='blue')+ggtitle("95% CI")
  ggsave("fig8B.eps", device=cairo_ps)
}

###Some numerical values sometimes used### 
# print(c("S(dmax) =", E_A_mean[cols]))# S(end)
# print(c("I(dmax)=", E_I_mean[cols]))# I(end)
# rho<-(E_A_mean-E_I_mean)/(E_I_mean+.00001)# relative difference
# AA=which.max(abs(rho))#index at which |rho| is a maximum
# print(c("maxDELTA=",rho[AA]))
# print(c("at", 50*(doses[AA]+doses[AA-1])), "cGy")
# #plot(doses,rho,type='l')
# I20=min(which(E_I_mean>20))
# d20=50*(doses[I20]+doses[I20-1])
# print(c("d20=",d20, "cGy"))
# a50=min(which(doses>.50))
# print(c("S(50)=", .5*(E_A_mean[a50]+E_A_mean[a50-1])))
# print(c("I(50)=", .5*(E_I_mean[a50]+E_I_mean[a50-1])))
# I30=min(which(E_I_mean>30))
# d30=50*(doses[I30]+doses[I30-1])
# print(c("d30=",d30))

###plot kappa###
# k1=parameter[[2]]$k1
# k2=parameter[[2]]$k2
# L=5*(1:200)
# kkappa=(k1*L)*exp(-k2*L)
# plot(L,kkappa)
