# Copyright:    (C) 2017-2018 Sachs Undergraduate Research Apprentice Program
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3.
# Filename:     plots.R 
# Purpose:      Concerns radiogenic mouse Harderian gland tumorigenesis. 
#               Contains code to generate figures. It is part of the 
#               source code for the NASAmouseHG project.
# Contact:      Rainer K. Sachs 
# Website:      https://github.com/sachsURAP/NASAmouseHG
# Mod history:  17 May 2018
# Details:      See hgData.R for further licensing, attribution, references, 
#               and abbreviation information.

source("monteCarlo.R") #  Load Monte Carlo

library(ggplot2) # Ribbon plot functionality
library(grid)  # Plot grids
library(Hmisc) # Error bars

forty_cGy <- 0:41 
sixty_cGy <- 0:61
seventy_cGy <- 0:71
hundred_cGy <- 0:101
forty_nine_cGy <- 0:50

#===============================================================================
#=== Shape of DER for Fe 600 MeV/u. Paper Fig. 3 (was Fig. 2.2.4.1)  4/23/18 ===
#===============================================================================
d <- 0.01 * 0:16000 # RKD to EGH. dose notation inconsistent with the above forty_cGy etc.
# RKS to EGH. Ideally, the next three plots would align horizontally as panels 
# A, B, and C, about 4" total width with room for ~12 point annotations. 
# However, I will in any case do a lot of fine tuning in Illustrator so even 
# just the right curve shapes panels of almost equal height, and the legend 
# shown are about good enough.
prevalence <- calibrated_HZE_nte_der(d, 193) # RKS to EGH. I suggest changing the name
#of this function to "calibrated_HZE_nte_der and similarly for other functions. # EGH 26 Apr: Done. ADDRESSED
plot(d, prevalence, type = 'l', bty = 'u')
legend(x = "bottomright", 
       legend = "dose in centiGy; Fe 193 zoom in twice",
       cex = 0.6, inset = 0.025)
d <- 2 * 10^-6 * 0:16000
prevalence <- calibrated_HZE_nte_der(d, 193)
plot(d, prevalence, type = 'l', bty='u')
d <- 10^-7 * 0:16000
prevalence <- calibrated_HZE_nte_der(d, 193)
plot(d, prevalence, type = 'l', bty = 'u')

#===============================================================================
#== Low LET data, error bars, and DER. Paper Fig. 5 (was Fig. 3.1.1)  5/17/18 ==
#===============================================================================
errbar(ion_data[1:4, "dose"], ion_data[1:4, "Prev"],yplus=ion_data[1:4, "Prev"]+1.96*ion_data[1:4, "SD"], 
  yminus=ion_data[1:4, "Prev"]-1.96*ion_data[1:4, "SD"], pch=19,cap=0.02,xlim=c(0,700), 
  ylim=c(0,1), bty='l',col='red', errbar.col = 'red', ann=FALSE) #  RKS: proton data points
errbar(ion_data[5:12, "dose"], ion_data[5:12, "Prev"],yplus=ion_data[5:12, "Prev"]+1.96*ion_data[5:12, "SD"],yminus=ion_data[5:12, "Prev"]-1.96*ion_data[5:12, "SD"], pch = 19,cap=0.02, add=TRUE) #  RKS: Helium data points
ddose=0:700 # RKS to EGH: This line and the next need work but do function  # EGH 26 Apr: Ok. ADDRESSED
lines(ddose, 1-exp(-coef(summary(low_LET_model, correlation = TRUE))[1]*ddose)) # RKS to EGH of course 0.00153 is actually from a summary()
legend(x = "topleft", legend = c("protons","4He"), col = c("red", "black"), pch = c(19,19), cex = 1, inset = 0.025)

ddose=0:82 # RKS to EGH: This line and the next need work but do function  # EGH 26 Apr: Ok. ADDRESSED
plot(ddose, 1-exp(-coef(summary(low_LET_model, correlation = TRUE))[1]*ddose), type='l',lwd=2, ann=FALSE, xlim=c(0,80),ylim=c(-.02,.35), bty='u')
errbar(ion_data[1:2, "dose"], ion_data[1:2, "Prev"],yplus=ion_data[1:2, "Prev"]+1.96*ion_data[1:2, "SD"],yminus=ion_data[1:2, "Prev"]-1.96*ion_data[1:2, "SD"], pch = 19,cap=0.05, add=TRUE, col='orange',errbar.col = 'orange', ann=FALSE, xlim=c(0,80),ylim=c(-.02,.35),lwd=2, bty='u') #  RKS: proton data points
errbar(ion_data[5:7, "dose"], ion_data[5:7, "Prev"],yplus=ion_data[5:7, "Prev"]+1.96*ion_data[5:7, "SD"],yminus=ion_data[5:7, "Prev"]-1.96*ion_data[5:7, "SD"], pch = 19,cap=0.05, add=TRUE, col='black',errbar.col = 'black', ann=FALSE, lwd=2) #alpha particles
errbar(60, .081,yplus=.081+.09,yminus=.081-.09, pch = 19,cap=0.05, add=TRUE, col='red',errbar.col = 'red', ann=FALSE, lwd=2)

ddose=0:82 # RKS to EGH: This line and the next need work but do function  # EGH 26 Apr: Ok. ADDRESSED
lines(ddose, 1 - exp(-coef(summary(low_LET_model, correlation = TRUE))[1]*ddose), lwd=2) # RKS to EGH of course 0.00153 is actually from a summary()
#lines(ddose,1 - exp(-coef(summary(low_LET_model, correlation = TRUE))[1] * ddose)) # Replacement for preceeding line # EGH to RKS please verify accuracy
#RKS to EGH 5/19/18 I checked that coef(summary(low_LET_model, correlation = TRUE))[1] is .00153 but did  not yet rerun monte_carlo.R
#======================== PLOTS  ==========================#
# Plot 1 : one HZE one low-LET; RKS->EH: always NTE & TE rather than TE-only in 
# minor paper. For plots 1 and 2 I can and will do the key box specifying which 
# curve is which by hand. When more than 2 ions are involved the key should be 
# part of the script that generates the figure. # EGH 26 Apr: Noted. ADDRESSED
d <- 1 * 0:302
r <- c(.2, .8) #Proportions. Next plot IDERs and MIXDER. 
plot(x = d, y = calibrated_HZE_nte_der(dose = d, L = 193), type = "l", xlab = "dose", ylab = "HG", bty = 'u', col = 'green', lwd = 2) 
lines(x = d, y = calibrated_low_LET_der(d, 0), col = 'orange', lwd = 2) 
lines(x = d, y = calculate_id(d, 193, r, lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
lines(x = d, y = calibrated_HZE_nte_der(dose = .2 * d, L = 193) + calibrated_low_LET_der(.8 * d, 0), lty= 3) # SEA S(d)


# Plot 2 : one HZE one low-LET; NTE & TE rather than TE-only always in minor paper
d <- 1 * 0:302 
r <- c(.8, .2) # Proportions. Next plot IDERs and MIXDER. 
plot(x = d, y = calibrated_HZE_nte_der(dose = d, L = 193), type = "l", xlab = "dose", ylab = "HG", bty = 'u', col = 'green', lwd = 2) 
lines(x = d, y = calibrated_low_LET_der(d, 0), col = 'orange', lwd = 2) 
lines(x = d, y = calculate_id(d, 193, r, lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
lines(x = d, y = calibrated_HZE_nte_der(dose = .8 * d, L = 193)+calibrated_low_LET_der(.2 *d, 0), lty= 3) # SEA S(d)

# Plot 3: four HZE; NTE
dose_vector <- c(0:100)
plot(calculate_id(dose_vector, c(25, 70, 190, 250), rep(0.25, 4)), type = 'l', col = 'red', bty = 'l', ann = 'F') #  I(d) plot
SEA <- function(dose) {
  return(calibrated_HZE_nte_der(dose / 4, 25) + 
           calibrated_HZE_nte_der(dose / 4, 70) + 
           calibrated_HZE_nte_der(dose / 4, 190) + 
           calibrated_HZE_nte_der(dose / 4, 250))
}
lines(dose_vector, SEA(dose_vector), lty = 2)
lines(dose_vector, calibrated_HZE_nte_der(dose_vector, 190), col = 'green') # component 4
lines(dose_vector, calibrated_HZE_nte_der(dose_vector, 250), col = 'green') # component 3
lines(dose_vector, calibrated_HZE_nte_der(dose_vector, 70), col = 'green') # component 2
lines(dose_vector, calibrated_HZE_nte_der(dose_vector, 25), col = 'green') # component 1


# Plot 4: two HZE; NTE; one low-LET
d <- seq(0, 100, 1)
plot(x = d, y = calibrated_HZE_nte_der(dose = d, L = 193), type = "l", xlab = "dose", ylab = "HG", bty = 'l', col = 'green', lwd = 2)
lines(x = d, y = calibrated_HZE_nte_der(d, 70), col = 'green', lwd = 2) # component 3
lines(x = d, y = calibrated_low_LET_der(d, 0), col = 'green', lwd = 2)
lines(x = d, y = calculate_id(d, c(70, 193), c(1/20, 1/20, 9/10), lowLET = TRUE)[, 2], col = 'red', lwd = 2)


#==================================================================#
#=========================== PAPER PLOTS ==========================#
#==================================================================#
# Fig. 3.2.3 - Fe56 (600 MeV/u), Si28, and corresponding IEA and SEA MIXDERS. 
plot(x = forty_cGy, y = calculate_SEA(forty_cGy, c(70, 193), c(1/2, 1/2), n = 2), bty='l', type = "l",  xlab = "Dose (Gy)", ylab = "HG Prevalence", col = "black", lwd = 2, lty = 2) 
lines(x = forty_cGy, y = calibrated_HZE_nte_der(dose = forty_cGy, L = 70), col = "cyan", lwd = 2)
lines(x = forty_cGy, y = calibrated_HZE_nte_der(dose = forty_cGy, L = 193), col = "orange", lwd = 2)
lines(x = forty_cGy, y = calculate_id(forty_cGy, c(70, 193), c(0.5 , 0.5), model = "NTE", lowLET = FALSE)[, 2], col = "red", lwd = 2) # I(d)
legend(x = "topleft", legend = c("Fe56 (600 MeV/u), NTE NTE-TE IDER", "Si28 HZE NTE-TE IDER", "IEA MIXDER (50% Fe56, 50% Si28)", "SEA MIXDER (50% Fe56, 50% Si28)"),
       col = c("orange", "cyan", "red", "black"), lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 2), cex = 0.4, inset = 0.05)


# Fig. for a mixture of 60 cGy H1 (protons) with 40 cGy Si;
plot(x = hundred_cGy, y = calibrated_HZE_nte_der(dose = hundred_cGy, L = 70), type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence", bty = 'l', col = "darkcyan", lwd = 2, xaxs="i")
lines(x = hundred_cGy, y = calibrated_low_LET_der(dose = hundred_cGy, L = 0.4), col = "cyan", lwd = 2)
lines(x = hundred_cGy, y = calculate_SEA(hundred_cGy, c(0.4, 70), c(0.6, 0.4), lowLET = TRUE, n = 2), col = "black", lwd = 2, lty = 2)
lines(x = hundred_cGy, y = calculate_id(hundred_cGy, c(0.4, 70), c(0.6 , 0.4), model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 100, lwd = 2)
legend(x = "topleft", legend = c("Si28 HZE NTE-TE IDER", "H1 Low-LET NTE-TE IDER","IEA MIXDER (60% H1, 40% Si28)", "SEA MIXDER (60% H1, 40% Si28)"),
       col = c("darkcyan", "cyan", "red", "black"), lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 2), cex = 0.5, inset = 0.025)

# Fig. for a mixture of 40 cGy H1 with 30 cGy Fe56 at 600 MeV/u;
plot(x = seventy_cGy, y = calibrated_HZE_nte_der(dose = seventy_cGy, L = 193), col = "darkcyan", type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence", bty = 'l', lwd = 2, lty = 1, xaxs="i")
lines(x = seventy_cGy, y = calibrated_low_LET_der(dose = seventy_cGy, L = 0.4), col = "cyan", lwd = 2)
lines(x = seventy_cGy, y = calculate_SEA(seventy_cGy, c(0.4, 193), c(4/7, 3/7), lowLET = TRUE, n = 2), col = "black", lty = 2, lwd = 2)
lines(x = seventy_cGy, y = calculate_id(seventy_cGy, c(0.4, 193), c(4/7 , 3/7), model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 70, lwd = 2)
legend(x = "topleft", legend = c("Fe56 (600 MeV/u), HZE NTE-TE IDER", "H1 Low-LET NTE-TE IDER","IEA MIXDER (57% H1, 43% Si28)", "SEA MIXDER (57% H1, 43% Si28)"),
       col = c("darkcyan", "cyan", "red", "black"), lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 2), cex = 0.5, inset = 0.025)


# Fig. for a mixture of all 7 HZE ions; total dose 49 cGy, each ion gets 7 cGy;
# assume Hi HZE implies Z > 3
plot(x = forty_nine_cGy, y = calculate_SEA(forty_nine_cGy, c(25, 70, 100, 193, 250, 464, 953), c(1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7)), 
     type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence", bty = 'l', col = "black", lwd = 2, lty = 2, axes = FALSE)

lines(x = forty_nine_cGy, y = calibrated_HZE_nte_der(dose = forty_nine_cGy, L = 25), col = "pink", lwd = 2)
lines(x = forty_nine_cGy, y = calibrated_HZE_nte_der(dose = forty_nine_cGy, L = 70), col = "orange", lwd = 2)
lines(x = forty_nine_cGy, y = calibrated_HZE_nte_der(dose = forty_nine_cGy, L = 100), col = "yellow", lwd = 3.5)
lines(x = forty_nine_cGy, y = calibrated_HZE_nte_der(dose = forty_nine_cGy, L = 193), col = "green", lwd = 3)
lines(x = forty_nine_cGy, y = calibrated_HZE_nte_der(dose = forty_nine_cGy, L = 250), col = "blue", lwd = 2)
lines(x = forty_nine_cGy, y = calibrated_HZE_nte_der(dose = forty_nine_cGy, L = 464), col = "purple", lwd = 2, lty = 2)
lines(x = forty_nine_cGy, y = calibrated_HZE_nte_der(dose = forty_nine_cGy, L = 953), col = "violet", lwd = 2)

lines(x = forty_nine_cGy, y = calculate_id(forty_nine_cGy, c(25, 70, 100, 193, 250, 464, 953),
                                           c(1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7), model = "NTE", lowLET = FALSE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 49, lwd = 1)
axis(2, c(-1, 0, 1, 2, 3, 4, 5, 6, 7))
axis(1, c(0, 10, 20, 30, 40, 49), xaxs = "i")
legend(x = "topleft", legend = c("Ne20 NTE-TE IDER", "Si28 NTE-TE IDER", 
                                 "Ti48 NTE-TE IDER", "Fe56 (600 MeV/u) NTE-TE IDER", 
                                 "Fe56 (300 MeV/u) NTE-TE IDER", "Nb93 NTE-TE IDER",
                                 "La139 NTE-TE IDER",
                                 "IEA MIXDER (Equally Distributed)", "SEA MIXDER (Equally Distributed)"),
       col = c("pink", "orange", "yellow", "green", "blue", "purple", "violet", "red", "black"), 
       lwd = c(2, 2, 2, 2, 2, 2, 2, 2, 2), 
       lty = c(1, 1, 1, 1, 1, 1, 1, 1, 2), cex = 0.3, inset = 0.0125)


# Fig. 3.2.2 
# h1, 0.4, 80%
# L=c(40, 110, 180, 250), proportions (5,5,5,5)%
forty_cGy = c(.001*0:9, .01*1:9, .1*1:9, 1:41) # RKS this plot needed more detail near dose zero
plot(x = forty_cGy, y = calculate_SEA(forty_cGy, c(0.4, 40, 110, 180, 250), c(.8, rep(.05,4)), lowLET = TRUE), 
     type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence", bty = 'l', col = "black", lwd = 2, lty = 2)

lines(x = forty_cGy, y = calibrated_low_LET_der(dose = forty_cGy, L = 0.4), col = "orange", lwd = 2)
lines(x = forty_cGy, y = calibrated_HZE_nte_der(dose = forty_cGy, L = 40), col = "green", lwd = 2)
lines(x = forty_cGy, y = calibrated_HZE_nte_der(dose = forty_cGy, L = 110), col = "purple", lwd = 2)
lines(x = forty_cGy, y = calibrated_HZE_nte_der(dose = forty_cGy, L = 180), col = "blue", lwd = 2)
lines(x = forty_cGy, y = calibrated_HZE_nte_der(dose = forty_cGy, L = 250), col = "aquamarine2", lwd = 2)

lines(x = forty_cGy, y = calculate_id(forty_cGy, c(0.4, 40,110, 180, 250), c(.8, rep(.05,4)), model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
legend(x = "topleft", legend = c("Low-LET","L=40", "L=110", "L=180", "L=250", 
                               "I(d)", "S(d)"),
       col = c("orange", "green","purple","blue", "aquamarine2", "red", "black"), 
       lwd = c(2, 2, 2, 2, 2, 3, 2), 
       lty = c(1, 1, 1, 1,  1, 1, 2), cex = 0.3, inset = 0.025)
forty_cGy=0:41  # RKS: replaced above more detailed forty_cGy with its former valued in case that is needed below

# RKS 5/9/18. I took out an entire incorrect plot here. The problem was having 2 different low_LET components. 
# The script incorrectly interprets the second low LET component as having non-zero NTE effects.

# 8-panel plot: 7 HZE IDERS, one low-LET HZE IDER in bottom right
par(mfrow = c(2, 4))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
LET_values <- c(25, 70, 100, 193, 250, 464, 953)
IDER_names <- c("He4", "Ne20", "Si28", "Ti48", "Fe56 (600 MeV/u)", "Nb93", "La139")
i = 1
while (i < length(LET_values) + 1) {
  plot(x = hundred_cGy, y = calibrated_HZE_nte_der(dose = hundred_cGy, L = LET_values[i]),
       xlim = c(0, 100), ylim = c(0, 1),
       type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l',
       col = "red", lwd = 2, axes = FALSE, ann = FALSE)
  
  if (i %in% c(5, 6, 7))
    axis(1, col = "black", col.axis = "black", at = seq(0, 100, 10))
  if (i %in% c(1, 5))
    axis(2, col = "black", col.axis = "black", at = seq(0, 1, .1))
  grid(col = "gray60")
  box(col = "gray60", lwd = 1.5)
  title(IDER_names[i], line = -2, cex = 0.5)
  i = i + 1
}

plot(x = hundred_cGy, y = calibrated_low_LET_der(dose = hundred_cGy, L = 0.4),
     xlim = c(0, 100), ylim = c(0, 1),
     type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l',
     col = "red", lwd = 2, axes = FALSE, ann = FALSE)
axis(1, col = "black", col.axis = "black", at = seq(0, 100, 10))
box(col = "gray60", lwd = 1.5)
grid(col = "gray60")
title("H1", line = -2, cex = 0.5)

mtext("Dose (cGy)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,
      col = "black")
mtext("HG Prevalence (%)", side = 2, outer = TRUE, cex = 0.7, line = 2.2,
      col = "black")
par(mfrow = c(1, 1))

#==============================================================================#
#====================== Confidence Interval Ribbon Plots ======================#
#==============================================================================#

#  Yimin's original test confidence interval plot
r <- c(0.05, 0.05, 0.05, 0.05, 0.8)
L <- c(25, 70, 100, 193)
calib_ci_test <- simulate_monte_carlo(200, hundred_cGy, L, r, model = "NTE")
ci_data <- data.frame(dose = hundred_cGy,
                      monteCarloBottom = calib_ci_test$monte_carlo[1, ],
                      monteCarloTop = calib_ci_test$monte_carlo[2, ], 
                      i = calculate_id(hundred_cGy, L, r, model = "NTE", 
                                               lowLET = TRUE)[, 2])
ci_plot <- ggplot(data = ci_data, aes = fill) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position="right") +
      labs(title = "Confidence Intervals", x = "Dose (cGy)", y = "HG Prevalence (%)") +
      geom_ribbon(aes(dose, ymin = monteCarloBottom, ymax = monteCarloTop), fill = "aquamarine2", alpha = 0.8) +
      scale_fill_discrete(name="Type",
                          breaks=c("pink"),
                          labels=c("Monte Carlo")) +
      geom_line(aes(dose, y = i), col = "red", size = 1) #  I(d) in red
      
ci_plot

#=============== Correlated vs Uncorrelated CI Overlay Plots ==================#

############### FIGURE 3.2.4 ############# 
# We use the MIXDER shown as Figure 3.2.4 in MS06. This consists of all seven
# HZE ions in our dataset, with equally distributed dosage.

#  Declare ratios and LET values for plot
ratios <- c(1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7)
LET_vals <- c(25, 70, 100, 193, 250, 464, 953)

#  We begin with the calibrated plot
calib_ci_3.2.4 <- simulate_monte_carlo(200, hundred_cGy, LET_vals, ratios, model = "NTE")

#  We now calculate the uncalibrated Monte Carlo
uncorr_ci_3.2.4 <- simulate_monte_carlo(200, hundred_cGy, LET_vals, ratios, model = "NTE", vcov = FALSE)

#  Construct a data.frame for ease of use with ggplot2
ci_data <- data.frame(dose = hundred_cGy,
                      #  Monte Carlo values
                      corrBottom = calib_ci_3.2.4$monte_carlo[1, ],
                      corrTop = calib_ci_3.2.4$monte_carlo[2, ],
                      uncorrTop = uncorr_ci_3.2.4$monte_carlo[1, ],
                      uncorrBottom = uncorr_ci_3.2.4$monte_carlo[2, ],
                      
                      #  DER values
                      ne = calibrated_HZE_nte_der(dose = hundred_cGy, L = 25),
                      si = calibrated_HZE_nte_der(dose = hundred_cGy, L = 70),
                      ti = calibrated_HZE_nte_der(dose = hundred_cGy, L = 100),
                      fe_six = calibrated_HZE_nte_der(dose = hundred_cGy, L = 193),
                      fe_three = calibrated_HZE_nte_der(dose = hundred_cGy, L = 250),
                      nb = calibrated_HZE_nte_der(dose = hundred_cGy, L = 464),
                      la = calibrated_HZE_nte_der(dose = hundred_cGy, L = 953),
                      
                      #  I(d)
                      i = calculate_id(hundred_cGy, LET_vals, ratios,
                                       model = "NTE", lowLET = FALSE)[, 2])
#  Plotting call
ci_plot <- ggplot(data = ci_data, aes = fill) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +
  labs(title = "Confidence Intervals", x = "Dose (cGy)", y = "HG Prevalence (%)") +
  
  #  Ribbon plot of both confidence intervals
  geom_ribbon(aes(dose, ymin = uncorrBottom, ymax = uncorrTop), fill = "aquamarine2", alpha = 1) + #  Uncorrelated in pink
  geom_ribbon(aes(dose, ymin = corrBottom, ymax = corrTop), fill = "yellow", alpha = .7) + #  Correlated in dull blue
  
  # scale_fill_discrete(name="Type",
                      # breaks=c("pink", "blue"),
                      # labels=c("Correlated Monte Carlo", "Uncorrelated Monte Carlo")) +
  
  #  DER plots
  geom_line(aes(dose, y = ne), col = "blue", size = 1) + #  neon in yellow
  geom_line(aes(dose, y = si),  col = "orange", size = 1) + #  silicon in orange
  geom_line(aes(dose, y = ti),  col = "green", size = 2) + # titanium in green
  geom_line(aes(dose, y = fe_six),  col = "purple", size = 1) + #  iron 600 in purple
  geom_line(aes(dose, y = fe_three),  col = "violet", size = 1) + #  iron 300 in violet
  geom_line(aes(dose, y = nb),  col = "darkcyan", size = 1) + #  niobium in dark cyan 
  geom_line(aes(dose, y = la),  col = "darkorange", size = 1) + #  lanthanum in black 
  
  # I(d) plot
  geom_line(aes(dose, y = i), col = "red", size = 1) #  I(d) in red

ci_plot # Print figure


############### FIGURE 3.2.3 #############  
# Fe56 (600 MeV/u) and Si28 in equal proportions

#  Declare ratios and LET values for plot
ratios <- c(1/2, 1/2)
LET_vals <- c(193, 70)
  
#  We begin with the calibrated plot
calib_ci_3.2.3 <- simulate_monte_carlo(200, hundred_cGy, LET_vals, ratios, model = "NTE")

#  We now calculate the uncalibrated Monte Carlo
uncorr_ci_3.2.3 <- simulate_monte_carlo(200, hundred_cGy, LET_vals, ratios, model = "NTE", vcov = FALSE)

ci_data <- data.frame(dose = hundred_cGy,
                      #  Monte Carlo values
                      corrBottom = calib_ci_3.2.3$monte_carlo[1, ],
                      corrTop = calib_ci_3.2.3$monte_carlo[2, ],
                      uncorrTop = uncorr_ci_3.2.3$monte_carlo[1, ],
                      uncorrBottom = uncorr_ci_3.2.3$monte_carlo[2, ],
                      
                      #  DER values
                      si = calibrated_HZE_nte_der(dose = hundred_cGy, L = 70),
                      fe_six = calibrated_HZE_nte_der(dose = hundred_cGy, L = 193),
                      
                      #  I(d)
                      i = calculate_id(hundred_cGy, LET_vals, ratios,
                                       model = "NTE", lowLET = TRUE)[, 2])

#  Plotting call
ci_plot <- ggplot(data = ci_data, aes = fill) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +
  labs(title = "Confidence Intervals", x = "Dose (cGy)", y = "HG Prevalence (%)") +
  
  #  Ribbon plot of both confidence intervals
  geom_ribbon(aes(dose, ymin = uncorrBottom, ymax = uncorrTop), fill = "aquamarine2", alpha = 1) +  #  Uncorrelated in pink
  geom_ribbon(aes(dose, ymin = corrBottom, ymax = corrTop), fill = "yellow", alpha = .7) + #  Correlated in dull blue
  
  # scale_fill_discrete(name="Type",
  # breaks=c("pink", "blue"),
  # labels=c("Correlated Monte Carlo", "Uncorrelated Monte Carlo")) +
  
  # DER plots
  geom_line(aes(dose, y = si),  col = "orange", size = 1) + #  iron in green
  geom_line(aes(dose, y = fe_six),  col = "darkblue", size = 1) + #  silicon in dark green
  
  # I(d) plot
  geom_line(aes(dose, y = i), col = "red", size = 1) #  I(d) in red
  
ci_plot # Print figure
