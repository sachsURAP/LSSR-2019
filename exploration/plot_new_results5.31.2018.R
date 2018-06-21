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
errbar(60, 0.081,yplus=.081+.09,yminus=.081-.09, pch = 19,cap=0.05, add=TRUE, col='red',errbar.col = 'red', lwd=3)

#========= Figure 10A for proposed new experiment June 18, 2018 =========#  
# Fe56 (600 MeV/u) and Ti in equal proportions for a total of 50 cGy
# Declare ratios and LET values for plot

# ratios <- c(1/2, 1/2)
# LET_vals <- c(193, 70)
# d10A = c(0.01*0:9,0.1*1:9,1:60)
# # # We use the plot that takes adjustable parameter correlations into account
# corr_ci_10B <- simulate_monte_carlo(n=500, d10A, LET_vals, ratios, model = "NTE") 
# # # The first argument, n, is the number of Monte Carlo repeats. Increase for
# # #greater accuracy. Decrease to speed up the program.
# # # Construct a data.frame for ease of use with ggplot2 if ggplot2 is used
# ci_data <- data.frame(dose = d10A,
#                        #  Monte Carlo values
#                        corrBottom = corr_ci_10A$monte_carlo[1, ],
#                        corrTop = corr_ci_10A$monte_carlo[2, ], # 
#                        
#                        # one-ion DERs for comparison
#                        fe_six = calibrated_HZE_nte_der(dose = d10A, L = 193),
#                        ti = calibrated_HZE_nte_der(dose = d10A, L = 100),
# # #                       
# # #                       #  IEA basline mixture DER I(d), denoted by id below
#                        i = calculate_id(d10A, LET_vals, ratios,
#                                         model = "NTE")[, 2]) 
# # # 
# # #  We make the ribbon plot for correlated parameters
# plot(c(0,61), c(0, 0.45), col="white", bty='L', ann=FALSE) #just sets plot area
# polygon(x=c(d10A,rev(d10A)),y=c(ci_data[,"corrTop"], rev(ci_data[,"corrBottom"])),
#          xpd=-1,col="yellow",lwd=.4,border="orange") ## narrow CI ribbon
# lines(ci_data[,"dose"],ci_data[,"ti"],col='brown', lwd = 2) # Si DER
# lines(ci_data[,"dose"],ci_data[,"fe_six"],col='blue', lwd = 2) #Fe
# lines(ci_data[,"dose"],ci_data[,"i"],col='red', lwd = 3) # I(d)
# abline(h=0.01*c(25,30,35,40))
# abline(v=c(45,55))

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


