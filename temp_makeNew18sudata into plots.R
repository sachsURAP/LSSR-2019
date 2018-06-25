# Copyright:    (C) 2017-2018 Sachs Undergraduate Research Apprentice Program
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3.
# Filename:     plot_new_results5.31.2018.R
# Purpose:      Run this file instead of running plots.R to get extra Figs.
# Contact:      Rainer K. Sachs 
# Website:      https://github.com/sachsURAP/NASAmouseHG
# Mod history:  13 Jun 2018
# Details:      See hgData.R for further licensing, attribution, references, 
#               and abbreviation information.

source("monteCarlo.R") #  Load Monte Carlo

library(ggplot2) # Ribbon plot functionality
library(grid)  # Plot grids
library(Hmisc) # Error bars

#========== Figure for Chang's new proton data point 5/22/2018=============#
ddose <- 0:82 
plot(ddose, 1-exp(-coef(summary(low_LET_model, correlation = TRUE))[1]*ddose),
     type='l',lwd=2, ann=FALSE, xlim=c(0,80),ylim=c(-.02,.35), bty='u')
errbar(ion_data[1:2, "dose"], ion_data[1:2, "Prev"],yplus=ion_data[1:2, "Prev"]+1.96*ion_data[1:2, "SD"],
       yminus=ion_data[1:2, "Prev"]-1.96*ion_data[1:2, "SD"], pch = 19,cap=0.03, add=TRUE, col='orange',
       errbar.col = 'orange', ann=FALSE,lwd=2) #  RKS: proton data points
errbar(ion_data[5:7, "dose"], ion_data[5:7, "Prev"],yplus=ion_data[5:7, "Prev"]+1.96*ion_data[5:7, "SD"],
       yminus=ion_data[5:7, "Prev"]-1.96*ion_data[5:7, "SD"], pch = 19,cap=0.03, add=TRUE, col='black',
       errbar.col = 'black', ann=FALSE, lwd=2) #alpha particles
errbar(60, 0.081,yplus=.081+.09,yminus=.081-.09, pch = 19,cap=0.05, add=TRUE, col='red',errbar.col = 'red', lwd=3)# use SEM not 1.96 SEM for all error bars perhaps -- both methods are valid.

#=================================================================#
#====== Fig. Si-p mixture July 2018===========#
d2Si=c(0.01*0:9,0.1*1:9,1:101)
corr_ci_2si <- simulate_monte_carlo(n=500, d2Si, c(70,0), c(.4,.6), model = "NTE") 
corrBottom = corr_ci_2si$monte_carlo[1, ]
corrTop = corr_ci_2si$monte_carlo[2, ] #
r =c(1-df_temp["si", "rLow"],df_temp["si", "rLow"])
plot(c(0,101), c(0, 0.4), col="white", ann=FALSE, bty='l')
polygon(x=c(d2Si,rev(d2Si)),y=c(corrTop, rev(corrBottom)),
        xpd=-1,col="yellow",lwd=.4,border="orange") 
lines(d2Si, calibrated_low_LET_der(d2Si, 0), lwd=2, col="blue")
lines(x = d2Si, y = calculate_id(d2Si, c(70, 0), r)[, 2], col = "red", lwd = 2) # I(d)
lines(d2Si, calibrated_HZE_nte_der(dose = d2Si, L = 70), lty=2, lwd=2, col = "blue")
points(100,df_temp["si","prev"], pch = 19)
lines(c(100,100), c(df_temp["si","prev"]-df_temp["si","SD"], df_temp["si","prev"]+df_temp["si","SD"]))

#====== Fig. Fe-p mixture July 2018 TE only===========#
d2FeT=0:70
corr_ci_2fe <- simulate_monte_carlo(n=100, d2FeT, c(193,0), c(3/7,4/7), model = "TE") 
corrBottom = corr_ci_2fe$monte_carlo[1, ]
corrTop = corr_ci_2fe$monte_carlo[2, ] #
r =c(1-df_temp["fe", "rLow"],df_temp["fe", "rLow"])
plot(c(0, 70), c(0, 0.6), col="white", ann=FALSE, bty='l')
polygon(x=c(d2FeT,rev(d2FeT)),y=c(corrTop, rev(corrBottom)),
        xpd=-1,col="yellow",lwd=.4,border="orange") 
lines(d2FeT, calibrated_low_LET_der(d2FeT, 0), lwd=2, col="blue")
lines(x = d2FeT, y = calculate_id(d2FeT, c(193, 0), r, model="TE")[, 2], col = "red", lwd = 3) # I(d)
lines(d2FeT, calibrated_HZE_te_der(dose = d2FeT, L = 193), lty=2, lwd=2, col = "blue")
points(71,df_temp["fe","prev"], pch= 19)
lines(c(71,71), c(df_temp["fe","prev"]-df_temp["fe","95CI"]/2, df_temp["fe","prev"]+df_temp["fe","95CI"]/2))
# In the line above, could use SD instead of 95CI/2; gives alternate, less stringent significance criterion.

#====== Fig. Fe-p mixture July 2018; as above but NTE instead of TE===========#
d2Fe=c(0.01*0:9,0.1*1:9,1:70)
corr_ci_2fe <- simulate_monte_carlo(n=500, d2Fe, c(193,0), c(3/7,4/7), model = "NTE") 
corrBottom = corr_ci_2fe$monte_carlo[1, ]
corrTop = corr_ci_2fe$monte_carlo[2, ] #
r =c(1-df_temp["fe", "rLow"],df_temp["fe", "rLow"])
plot(c(0, 70), c(0, 0.6), col="white", ann=FALSE, bty='l')
polygon(x=c(d2Fe,rev(d2Fe)),y=c(corrTop, rev(corrBottom)),
        xpd=-1,col="yellow",lwd=.4,border="orange") 
lines(d2Fe, calibrated_low_LET_der(d2Fe, 0), lwd=2, col="blue")
lines(x = d2Fe, y = calculate_id(d2Fe, c(193, 0), r)[, 2], col = "red", lwd = 3) # I(d)
lines(d2Fe, calibrated_HZE_nte_der(dose = d2Fe, L = 193), lty=2, lwd=2, col = "blue")
points(71,df_temp["fe","prev"], pch= 19)
lines(c(71,71), c(df_temp["fe","prev"]-df_temp["fe","SD"], df_temp["fe","prev"]+df_temp["fe","SD"]))

#==============================================================================#
#== Fig. 3N. Fe56 (600 MeV/u), Si28 Equal Doses, Confidence Interval new data added==#
#==============================================================================#
# Fe56 (600 MeV/u) and Si28 in equal proportions for a total of 40 cGy
# Declare ratios and LET values for plot
ratios <- c(1/2, 1/2)
LET_vals <- c(193, 70)
d3N <- c(0.01 * 0:9, 0.1 * 1:9, 1:40)
# We use the plot that takes adjustable parameter correlations into account
corr_fig_3N <- simulate_monte_carlo(n = 500, d3N, LET_vals, ratios, model = "NTE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d3N,
                      # Monte Carlo values
                      corrBottom = corr_fig_3N$monte_carlo[1, ],
                      corrTop = corr_fig_3N$monte_carlo[2, ], #
                      
                      # one-ion DERs for comparison
                      fe_six = calibrated_HZE_nte_der(dose = d3N, L = 193),
                      si = calibrated_HZE_nte_der(dose = d3N, L = 70),
                      
                      # IEA baseline mixture DER I(d), denoted by id below
                      i = calculate_id(d3N, LET_vals, ratios, model = "NTE")[, 2])

# We make the ribbon plot for correlated parameters
plot(c(0, 41), c(0, .40), col = "white", bty = 'L', ann = FALSE) # Set plot area
polygon(x = c(d3N, rev(d3N)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon
lines(ci_data[, "dose"], ci_data[, "si"], col = 'blue', lwd = 2) # Si DER
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 3) # I(d)
lines(ci_data[, "dose"], ci_data[, "fe_six"], col = 'blue',lwd=2, lty=2) # Fe
points(40,df_temp["both","prev"], pch= 19)
lines(c(40,40), c(df_temp["both","prev"]-df_temp["both","SD"], df_temp["both","prev"]+df_temp["both","SD"]))

#========= Figure 10B for proposed new experiment June 18, 2018 =========#  
# Fe56 (600 MeV/u), Si and Ti in equal proportions for a total of 60 cGy
# Declare ratios and LET values for plot

ratios <- c(1/3, 1/3, 1/3)
LET_vals <- c(193, 70, 100)
d10B = c(0.01*0:9,0.1*1:9,1:60)
# # We use the plot that takes adjustable parameter correlations into account
corr_ci_10B <- simulate_monte_carlo(n=500, d10B, LET_vals, ratios, model = "NTE") 
# # The first argument, n, is the number of Monte Carlo repeats. Increase for
# #greater accuracy. Decrease to speed up the program.
# # Construct a data.frame for ease of use with ggplot2 if ggplot2 is used
ci_data <- data.frame(dose = d10B,
                      #  Monte Carlo values
                      corrBottom = corr_ci_10B$monte_carlo[1, ],
                      corrTop = corr_ci_10B$monte_carlo[2, ], # 
                      
                      # one-ion DERs for comparison
                      fe_six = calibrated_HZE_nte_der(dose = d10B, L = 193),
                      ti = calibrated_HZE_nte_der(dose = d10B, L = 100),
                      si = calibrated_HZE_nte_der(dose = d10B, L = 70),
                      # #                       
                      # #                       #  IEA basline mixture DER I(d), denoted by id below
                      i = calculate_id(d10B, LET_vals, ratios,
                                       model = "NTE")[, 2]) 
# # 
# #  We make the ribbon plot for correlated parameters
plot(c(0,61), c(0, 0.45), col="white", bty='L', ann=FALSE) #just sets plot area
polygon(x=c(d10B,rev(d10B)),y=c(ci_data[,"corrTop"], rev(ci_data[,"corrBottom"])),
        xpd=-1,col="yellow",lwd=.4,border="orange") ## narrow CI ribbon
lines(ci_data[,"dose"],ci_data[,"si"],col='brown', lwd = 2) # Si DER
lines(ci_data[,"dose"],ci_data[,"fe_six"],col='blue', lwd = 2) #Fe
lines(ci_data[,"dose"],ci_data[,"ti"],col='green', lwd = 3) #Ti
lines(ci_data[,"dose"],ci_data[,"i"],col='red', lwd = 3) # I(d)
abline(h=0.35)
abline(v=c(45,60))

#========Now same idea with equal fluxes rather than doses ==========#
ratios <- c(.2692, .3448, .3859)/(.2692+.3448+.3859)
LET_vals <- c(70, 100,193)
d10C = c(0.01*0:9,0.1*1:9,1:60)
# We use the plot that takes adjustable parameter correlations into account
corr_ci_10C <- simulate_monte_carlo(n=500, d10C, LET_vals, ratios, model = "NTE") 
# # The first argument, n, is the number of Monte Carlo repeats. Increase for
# #greater accuracy. Decrease to speed up the program.
# # Construct a data.frame for ease of use with ggplot2 if ggplot2 is used
ci_data <- data.frame(dose = d10C,
                      #  Monte Carlo values
                      corrBottom = corr_ci_10C$monte_carlo[1, ],
                      corrTop = corr_ci_10C$monte_carlo[2, ], # 
                      
                      # one-ion DERs for comparison
                      fe_six = calibrated_HZE_nte_der(dose = d10C, L = 193),
                      ti = calibrated_HZE_nte_der(dose = d10C, L = 100),
                      si = calibrated_HZE_nte_der(dose = d10C, L = 70),
                      # #                       
                      # #                       #  IEA basline mixture DER I(d), denoted by id below
                      i = calculate_id(d10C, LET_vals, ratios,
                                       model = "NTE")[, 2]) 
# # 
# #  We make the ribbon plot for correlated parameters
plot(c(0,61), c(-.01, 0.65), col="white", bty='L', ann=FALSE) #just sets plot area
polygon(x=c(d10C,rev(d10C)),y=c(ci_data[,"corrTop"], rev(ci_data[,"corrBottom"])),
        xpd=-1,col="yellow",lwd=.4,border="orange") ## narrow CI ribbon
lines(ci_data[,"dose"],ci_data[,"si"],col='blue', lwd = 2,lty=2) # Si DER
lines(ci_data[,"dose"],ci_data[,"fe_six"],col='blue', lwd = 2) #Fe
lines(ci_data[,"dose"],ci_data[,"ti"],col='green', lwd = 3) #Ti
lines(ci_data[,"dose"],ci_data[,"i"],col='red', lwd = 3) # I(d)
abline(h=0.35) # near I(d) at 60 cGy
#abline(v=c(45,60))

#experimental error  bars
mousen=50
prevl=.01*c(14, 35, 55)
eerr=round(1.96*sqrt(prevl*(1-prevl)/mousen),4)
lines(c(45,45), c(prevl[1]+eerr[1], prevl[1]-eerr[1])) #approximates 95% CI
lines(c(60,60), c(prevl[2]+eerr[2], prevl[2]-eerr[2]))
lines(c(45,45), c(prevl[3]+eerr[3], prevl[3]-eerr[3]))

