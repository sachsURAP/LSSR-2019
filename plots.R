# Copyright:    (C) 2017-2018 Sachs Undergraduate Research Apprentice Program
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3.
# Filename:     plots.R 
# Purpose:      Concerns radiogenic mouse Harderian gland tumorigenesis. 
#               Contains code to generate many of the REBP paper figures. It is
#               part of the source code for the NASAmouseHG project.
# Contact:      Rainer K. Sachs 
# Website:      https://github.com/sachsURAP/NASAmouseHG
# Mod history:  18 Jun 2018
# Details:      See dataAndInfo.R for further licensing, attribution, 
#               references, and abbreviation information.

source("monteCarlo.R") # Load Monte Carlo.

library(Hmisc) # Error bars.

#==============================================================================#
#==== Fig. 1. SEA Synergy Theory is Unreliable for Quadratic Component DERs ===#
#==============================================================================#
d1 <- 0.01 * 0:100
d <- 0.5 * d1
E1 <- d1 ^ 2
E2 <- 2 * d1 ^ 2
SEA <- d ^ 2 + 2 * d ^ 2
plot(d1, E2, type = 'l', lwd = 3, bty = 'l', ann = FALSE)
lines(d1, E1, lwd = 3)
lines(d1, SEA, lwd = 3, lty = 2)

#==============================================================================#
#==================== Fig. 2. Convex, Concave, Standard =======================#
#==============================================================================#
d2 <- 0.01 * 0:200
a <- 2; b <- .6; c <- 4
E1 <- a * d2 + b * d2 ^ 2  # Convex
E2 <- a * d2  # Linear no-threshold (LNT) same initial slope
E3 <- c * (1 - (exp(- a * d2 / c))) # Concave, same initial slope
plot(d2, E1, type = 'l', lwd = 3, bty = 'l', ann = FALSE)
lines(d2, E2, lwd = 3)
lines(d2, E3, lwd = 3)

a <- 0.45; b <- 1 / 8; c <- 0.8
E1 <- b * d2 + 0.35 * b * d2 ^ 2  # Convex
E2 <- 0.7 * a * d2  # Linear no-threshold (LNT) 
E3 <- c * (1 - (exp(- 2 * a * d2 / c))) # Concave
plot(d2, E3, type = 'l', lwd = 3, bty = 'l', ann = FALSE)
lines(d2, E1, lwd = 3)
lines(d2, E2, lwd = 3)

#==============================================================================#
#==================== Fig. 3. Shape of DER for Fe 600 MeV/u. ==================#
#==============================================================================#
d3A <- c(0.01 * 0:9, 0.1 * 1:9, 1:150) 
prevalence <- calibrated_HZE_nte_der(d3A, 193) 
plot(d3A, prevalence, type = 'l', bty = 'u', lwd = 3, ann = FALSE)
# legend(x = "bottomright", 
#        legend = "dose in centiGy; Fe 193 zoom in twice",
#        cex = 0.4, inset = 0.025)

d3B <- 2 * 10 ^ -5 * 0:1600 # Zoom in by a factor of 10^4
prevalence <- calibrated_HZE_nte_der(d3B, 193)
plot(d3B, prevalence, type = 'l', bty = 'u', lwd = 3, ann = FALSE)

d3C <- 10 ^ -6 * 0:1600 # Zoom in by another factor of 20
prevalence <- calibrated_HZE_nte_der(d3C, 193)
plot(d3C, prevalence, type = 'l', bty = 'u', lwd = 3)

#==============================================================================#
#=============== Fig. 5. Low LET Data, Error Bars, and DER. ===================#
#==============================================================================#
ddose <- 0:701 
plot(c(0, 701), c(-.02, 1), pch = 19, col = 'white', ann = FALSE, bty = 'u')
lines(ddose, 1 - exp(- coef(summary(low_LET_model, correlation = TRUE))[1] * ddose), lwd = 2)

errbar(ion_data[5:12, "dose"], ion_data[5:12, "Prev"], 
       yplus  = ion_data[5:12, "Prev"] + 1.96 * ion_data[5:12, "SD"],
       yminus = ion_data[5:12, "Prev"] - 1.96 * ion_data[5:12, "SD"],
       pch = 19, cap = 0.02, add = TRUE, col = 'red', errbar.col = 'red', lwd = 2) # Alpha particle data
errbar(ion_data[1:4, "dose"], ion_data[1:4, "Prev"],
       yplus  = ion_data[1:4, "Prev"] + 1.96 * ion_data[1:4, "SD"], 
       yminus = ion_data[1:4, "Prev"] - 1.96 * ion_data[1:4, "SD"], 
       pch = 19, cap = 0.02, add = TRUE, col = 'black', errbar.col = 'black', lwd = 2) # Proton data
# legend(x = "topleft", legend = c("protons","4He"), col = c("red", "black"),
# pch = c(19,19), cex = 1, inset = 0.025)

#==============================================================================#
#======================= Fig. 6. One HZE Ion, One Low-LET =====================#
#==============================================================================#
# We will always use NTE & TE rather than TE-only for HZE in the minor paper. 
d <- 1 * 0:302
r <- c(.2, .8) # Dose proportions.
plot(x = d, y = calibrated_HZE_nte_der(dose = d, L = 193), type = "l", #  Now plot DERs 
     xlab = "dose", ylab = "HG", bty = 'u', col = 'green', lwd = 2) # Fe DER
lines(x = d, y = calibrated_low_LET_der(d, 0), col = 'orange', lwd = 2) # Low LET DER 
lines(x = d, y = calculate_id(d, c(193, 0), r)[, 2], col = "red", lwd = 3) # I(d)
lines(x = d, y = calibrated_HZE_nte_der(dose = .2 * d, L = 193) + 
        calibrated_low_LET_der(.8 * d, 0), lty = 2, lwd = 2) # SEA mixture baseline S(d)

r <- c(.8, .2) # Panel B, proportions reversed with low LET small  
plot(x = d, y = calibrated_HZE_nte_der(dose = d, L = 193), type = "l", 
     xlab = "dose", ylab = "HG", bty = 'u', col = 'green', lwd = 2) 
lines(x = d, y = calibrated_low_LET_der(d, 0), col = 'orange', lwd = 2) 
lines(x = d, y = calculate_id(d, c(193, 0), r)[, 2], col = "red", lwd = 3) # I(d)
lines(x = d, y = calibrated_HZE_nte_der(dose = .8 * d, L = 193) +
        calibrated_low_LET_der(.2 *d, 0), lty= 2, lwd=2) 

#==============================================================================#
#===================== Fig. 7. 80% LowLET and Four HZE Ions ===================#
#==============================================================================#
d7 <- c(0.01*0:9, 0.1*1:9, 1:41)
plot(c(0,41),c(0,0.35), col="white", xlab = "Dose (cGy)", ylab = "HG", bty = 'l')
lines(x = d7, y = calculate_SEA(d7, c(0.4, 40, 110, 180, 250), c(.8, rep(.05,4)),
                                lowLET = TRUE), col = "black", lwd = 2, lty = 2)
lines(x = d7, y = calibrated_low_LET_der(dose = d7, L = 0.4),
      col = "orange", lwd = 2) # low LET DER
lines(x = d7, y = calibrated_HZE_nte_der(dose = d7, L = 40),
      col = "green", lwd = 2) # HZE DER
lines(x = d7, y = calibrated_HZE_nte_der(dose = d7, L = 110),
      col = "purple", lwd = 2)
lines(x = d7, y = calibrated_HZE_nte_der(dose = d7, L = 180),
      col = "blue", lwd = 2)
lines(x = d7, y = calibrated_HZE_nte_der(dose = d7, L = 250),
      col = "aquamarine2", lwd = 2)
lines(x = d7, y = calculate_id(d7, c(0.4, 40, 110, 180, 250),
                               c(.8, rep(.05, 4)), model = "NTE")[, 2], col = "red", lwd = 3) # I(d)
legend(x = "topleft", legend = c("Low-LET","L=40", "L=110", "L=180", 
                                 "L=250", "I(d)", "S(d)"),
       col = c("orange", "green","purple","blue", "aquamarine2", "red", "black"), 
       lwd = c(2, 2, 2, 2, 2, 3, 2), 
       lty = c(1, 1, 1, 1,  1, 1, 2), cex = 0.3, inset = 0.025)

#==============================================================================#
#=============== Fig. 8. DERs: Fe56 (600 MeV/u), Si28, IEA, SEA ===============#
#==============================================================================#
d8 <- c(0.01 * 0:9, 0.1 * 1:9, 1:41)
plot(c(0, 40), c(0, 0.5), col = "white", bty = 'l', xlab = "Dose (cGy)", ylab = "HG")
lines(x = d8, y = calibrated_HZE_nte_der(dose = d8, L = 70), col = "cyan", lwd = 2)
lines(x = d8, y = calibrated_HZE_nte_der(dose = d8, L = 193), col = "orange", lwd = 2)
lines(x = d8, y = calculate_id(d8, c(70, 193), c(0.5 , 0.5), model = "NTE")[, 2],
      col = "red", lwd = 3) # I(d)
lines(x = d8, y = calculate_SEA(d8, c(70, 193), c(1/2, 1/2), n = 2), col = "black",
      lwd = 2, lty = 2)
abline(v = 40) # Unpublished results 6.10.2018
abline(h = .01 * c(17.7, 32.5, 47.3))
legend(x = "topleft", legend = c("Fe56 (600 MeV/u)", "Si28", "IEA", "SEA"),
       col = c("orange", "cyan", "red", "black"), lwd = c(2, 2, 2, 2),
       lty = c(1, 1, 1, 2), cex = 0.6, inset = 0.05)

#==============================================================================#
#============== Fig. 9. All 7 HZE Ions; Each Ion Contributes 7 cGy ============#
#==============================================================================#
d9 <- c(0.01 * 0:9, 0.1 * 1:9, 0.5 * 2:100)
plot(x = d9, y = calculate_SEA(d9, c(25, 70, 100, 193, 250, 464, 953), rep(1/7, 7)), 
     type = "l", xlab = "Dose (cGy)", ylab = "HG", bty = 'u', col = "black", lwd = 2, lty = 2)
lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 25), col = "pink", lwd = 2)
lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 70), col = "orange", lwd = 2)
lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 100), col = "aquamarine2", lwd = 4)
lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 193), col = "green", lwd = 2)
lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 250), col = "blue", lwd = 2)
lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 464), col = "purple", lwd = 2)
lines(x = d9, y = calibrated_HZE_nte_der(dose = d9, L = 953), col = "violet", lwd = 2)
lines(x = d9, y = calculate_id(d9, c(25, 70, 100, 193, 250, 464, 953),
                               rep(1/7, 7), model = "NTE")[, 2], col = "red", lwd = 3) # I(d)
legend(x = "topleft", legend = c("Ne20 NTE-TE IDER", "Si28 NTE-TE IDER", 
                                 "Ti48 NTE-TE IDER", "Fe56 (600 MeV/u) NTE-TE IDER", 
                                 "Fe56 (300 MeV/u) NTE-TE IDER", "Nb93 NTE-TE IDER",
                                 "La139 NTE-TE IDER", "IEA MIXDER (Equally Distributed)",
                                 "SEA MIXDER (Equally Distributed)"),
       col = c("pink", "orange", "aquamarine2", "green", "blue",
               "purple", "violet", "red", "black"), 
       lwd = c(2, 2, 2, 2, 2, 2, 2, 2, 2), 
       lty = c(1, 1, 1, 1, 1, 1, 1, 1, 2), cex = 0.3, inset = 0.0125)


#+++++++++++++++++++++ Confidence Interval Ribbon Plots +++++++++++++++++++++++#

#==============================================================================#
#====== Fig. 10. Fe56 (600 MeV/u), Si28 Equal Doses, Confidence Interval ======#
#==============================================================================#
# Fe56 (600 MeV/u) and Si28 in equal proportions for a total of 40 cGy
# Declare ratios and LET values for plot
ratios <- c(1/2, 1/2)
LET_vals <- c(193, 70)
d10 <- c(0.01 * 0:9, 0.1 * 1:9, 1:40)
# We use the plot that takes adjustable parameter correlations into account
corr_fig_10 <- simulate_monte_carlo(n = 500, d10, LET_vals, ratios, model = "NTE")
# The first argument, n, is the number of Monte Carlo repeats. Increase for
# greater accuracy. Decrease to speed up the program.
ci_data <- data.frame(dose = d10,
                       # Monte Carlo values
                       corrBottom = corr_fig_10$monte_carlo[1, ],
                       corrTop = corr_fig_10$monte_carlo[2, ], #

                       # one-ion DERs for comparison
                       fe_six = calibrated_HZE_nte_der(dose = d10, L = 193),
                       si = calibrated_HZE_nte_der(dose = d10, L = 70),

                       # IEA baseline mixture DER I(d), denoted by id below
                       i = calculate_id(d10, LET_vals, ratios, model = "NTE")[, 2])

# We make the ribbon plot for correlated parameters
plot(c(0, 41), c(0, .30), col = "white", bty = 'L', ann = FALSE) # Set plot area
polygon(x = c(d10, rev(d10)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
         xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon
lines(ci_data[, "dose"], ci_data[, "si"], col = 'brown', lwd = 2) # Si DER
lines(ci_data[, "dose"], ci_data[, "fe_six"], col = 'blue') # Fe
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 3) # I(d)

#==============================================================================#
#========= Fig. 11. Correlated vs Uncorrelated CI Overlay Plot ================#
#==============================================================================#
# Consists of all 7 HZE ions in our 5/20/2018 data set.
# Declare ratios and LET values for plot
ratios <- c(1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7)
LET_vals <- c(25, 70, 100, 193, 250, 464, 953)
d11 <- c(0.1 * 0:9, 1:50)

# We begin with the correlated plot
corr_fig_11 <- simulate_monte_carlo(n = 500, d11, LET_vals, ratios, model = "NTE")
# Comments for Fig. 10 apply with minor changes here and in some other lines
# We now calculate the uncorrelated Monte Carlo
uncorr_fig_11 <- simulate_monte_carlo(n = 500, d11, LET_vals, ratios, model = "NTE", vcov = FALSE)

ci_data <- data.frame(dose = d11,
                       # Monte Carlo values
                       corrBottom = corr_fig_11$monte_carlo[1, ],
                       corrTop = corr_fig_11$monte_carlo[2, ],
                       uncorrBottom = uncorr_fig_11$monte_carlo[1, ],
                       uncorrTop = uncorr_fig_11$monte_carlo[2, ],

                       # DER values
                       ne = calibrated_HZE_nte_der(dose = d11, L = 25),
                       si = calibrated_HZE_nte_der(dose = d11, L = 70),
                       ti = calibrated_HZE_nte_der(dose = d11, L = 100),
                       fe_six = calibrated_HZE_nte_der(dose = d11, L = 193),
                       fe_three = calibrated_HZE_nte_der(dose = d11, L = 250),
                       nb = calibrated_HZE_nte_der(dose = d11, L = 464),
                       la = calibrated_HZE_nte_der(dose = d11, L = 953),
                       i = calculate_id(d11, LET_vals, ratios, model = "NTE")[, 2])

# Plotting call. 
plot(c(0, 50), c(0, .40), col = "white", bty = 'u', ann = FALSE) # Next is broad CI ribbon
polygon(x = c(d11, rev(d11)), y = c(ci_data[, "uncorrTop"], rev(ci_data[, "uncorrBottom"])),
        xpd = -1, col = "orange", lwd = .5, border = "orange") # Wide CI

polygon(x = c(d11, rev(d11)), y = c(ci_data[, "corrTop"], rev(ci_data[, "corrBottom"])),
        xpd = -1, col = "yellow", border = "orange", lwd = .2) # Narrow CI

lines(ci_data[, "dose"], ci_data[, "fe_three"], col = 'black', lwd = 2)
lines(ci_data[, "dose"], ci_data[, "ne"], col = 'black', lwd = 2)
lines(ci_data[, "dose"], ci_data[, "nb"], col = 'black', lwd = 2)
lines(ci_data[, "dose"], ci_data[, "la"], col = 'black', lwd = 2)
lines(ci_data[, "dose"], ci_data[, "si"], col = 'black', lwd = 2)
lines(ci_data[, "dose"], ci_data[, "ti"], col = 'black', lwd = 2)
lines(ci_data[, "dose"], ci_data[, "fe_six"], col = 'black', lwd = 2)
lines(ci_data[, "dose"], ci_data[, "i"], col = 'red', lwd = 3)
