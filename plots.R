#   Filename: HGsynergyMain.R 
#   Purpose: Concerns radiogenic mouse HG tumorigenesis. Contians code to 
#            generate figures.

source("HGsynergyMain.R") # load data
source("synergyTheory.R") # load models
source("monteCarlo.R") # load Monte Carlo

library(ggplot2) #   Plotting


forty_cGy <- .01 * 0:40
sixty_cGy <- .01 * 0:60
seventy_cGy <- .01 * 0:70
hundred_cGy <- .01 * 0:100
forty_nine_cGy <- .01 * 0:49

# Fig. 3.2.1.1. - Fe56 (600 MeV/u), Si28, and corresponding IEA and SEA MIXDERS.
setEPS()
postscript("fe56_si28_nte.eps")
plot(x = forty_cGy, y = calculate_SEA(forty_cGy, 1/2, c(70, 195)), type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2, xaxs="i")
lines(x = forty_cGy, y = calib_HZE_nte_ider(dose = forty_cGy, L = 70), col = "cyan", lwd = 2)
lines(x = forty_cGy, y = calib_HZE_nte_ider(dose = forty_cGy, L = 195), col = "darkcyan", lwd = 2)
lines(x = forty_cGy, y = calculate_complex_id(r = c(0.5 , 0.5), L = c(70, 195), d = forty_cGy, model = "NTE", lowLET = FALSE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 0.4, lwd = 2)
legend(x = "topleft", legend = c("Fe56 (600 MeV/u), NTE NTE-TE IDER", "Si28 HZE NTE-TE IDER", "IEA MIXDER (50% Fe56, 50% Si28)", "SEA MIXDER (50% Fe56, 50% Si28)"),
       col = c("darkcyan", "cyan", "red", "black"), lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 2), inset = 0.05)
dev.off()


# Fig. for a mixture of 60 cGy H1 (protons) with 40 cGy Si;
setEPS()
postscript("h1_si28_nte.eps")
plot(x = hundred_cGy, y = calculate_SEA(hundred_cGy, c(0.6, .4), c(0.4, 70)), type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2, xaxs="i")
lines(x = hundred_cGy, y = calib_low_LET_ider(dose = hundred_cGy, L = 0.4), col = "cyan", lwd = 2)
lines(x = hundred_cGy, y = calib_HZE_nte_ider(dose = hundred_cGy, L = 70), col = "darkcyan", lwd = 2)
lines(x = hundred_cGy, y = calculate_complex_id(r = c(0.6 , 0.4), L = c(0.4, 70), d = hundred_cGy, model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 1, lwd = 2)
legend(x = "topleft", legend = c("Si28 HZE NTE-TE IDER", "H1 Low-LET NTE-TE IDER","IEA MIXDER (60% H1, 40% Si28)", "SEA MIXDER (60% H1, 40% Si28)"),
       col = c("darkcyan", "cyan", "red", "black"), lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 2), inset = 0.025)
dev.off()

# Fig. for a mixture of 40 cGy H1 with 30 cGy Fe56 at 600 MeV/u;
setEPS()
postscript("h1_fe56_nte.eps")
plot(x = seventy_cGy, y = calculate_SEA(seventy_cGy, r = c(4/7, 3/7), c(0.4, 195)), type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2, xaxs="i")
lines(x = seventy_cGy, y = calib_low_LET_ider(dose = seventy_cGy, L = 0.4), col = "cyan", lwd = 2)
lines(x = seventy_cGy, y = calib_HZE_nte_ider(dose = seventy_cGy, L = 195), col = "darkcyan", lwd = 2)
lines(x = seventy_cGy, y = calculate_complex_id(r = c(4/7 , 3/7), L = c(0.4, 195), d = seventy_cGy, model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 0.7, lwd = 2)
legend(x = "topleft", legend = c("Fe56 (600 MeV/u), HZE NTE-TE IDER", "H1 Low-LET NTE-TE IDER","IEA MIXDER (57% H1, 43% Si28)", "SEA MIXDER (57% H1, 43% Si28)"),
       col = c("darkcyan", "cyan", "red", "black"), lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 2), inset = 0.005)
dev.off()

# Fig. for a mixture of all 7 HZE ions; total dose 49 cGy, each ion gets 7 cGy;
setEPS()
postscript("all_hze_nte.eps")

# assume Hi HZE implies Z > 3
plot(x = forty_nine_cGy, y = calculate_SEA(forty_nine_cGy, r =  c(1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7), c(25, 70, 100, 195, 250, 464, 953)), 
     type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2, axes=FALSE)

lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 25), col = "pink", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 70), col = "orange", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 100), col = "yellow", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 195), col = "green", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 250), col = "blue", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 464), col = "purple", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 953), col = "violet", lwd = 2)

lines(x = forty_nine_cGy, y = calculate_complex_id(r = c(1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7), L = c(25, 70, 100, 195, 250, 464, 953),
                                                   d = forty_nine_cGy, model = "NTE", lowLET = FALSE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 0.49, lwd = 1)
axis(2, c(-0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
axis(1, c(-.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.49), xaxs = "i")
legend(x = "topleft", legend = c("Ne20 NTE-TE IDER", "Si28 NTE-TE IDER", 
                                 "Ti48 NTE-TE IDER", "Fe56 (600 MeV/u) NTE-TE IDER", 
                                 "Fe56 (300 MeV/u) NTE-TE IDER", "Nb93 NTE-TE IDER",
                                 "La139 NTE-TE IDER",
                                 "IEA MIXDER (Equally Distributed)", "SEA MIXDER (Equally Distributed)"),
       col = c("pink", "orange", "yellow", "green", "blue", "purple", "violet", "red", "black"), 
       lwd = c(2, 2, 2, 2, 2, 2, 2, 2, 2), 
       lty = c(1, 1, 1, 1, 1, 1, 1, 1, 2), cex = 0.55, inset = 0.0125)
dev.off()



# h1, 0.4, 60
# he4, 1.6, 20
# o16, 25, 10
# si28, 70, 2.5
# ti28, 100, 2.5
# fe56, 195, 5
setEPS()
postscript("big_mix_nte.eps")
plot(x = forty_nine_cGy, y = calculate_SEA(forty_nine_cGy, r = c(.6, .2, .1, .025, .025, .5), L = c(0.4, 1.6, 25, 70, 100, 195)), 
     type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2, xaxs = "i")

lines(x = forty_nine_cGy, y = calib_low_LET_ider(dose = forty_nine_cGy, L = 0.4), col = "orange", lwd = 2)
lines(x = forty_nine_cGy, y = calib_low_LET_ider(dose = forty_nine_cGy, L = 1.6), col = "yellow", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 25), col = "green", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 70), col = "blue", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 100), col = "purple", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 195), col = "violet", lwd = 2)

lines(x = forty_nine_cGy, y = calculate_complex_id(r = c(.6, .2, .1, .025, .025, .5), L =  c(0.4, 1.6, 25, 70, 100, 195),
                                                   d = forty_nine_cGy, model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 1, lwd = 1)
legend(x = "right", legend = c("H1 Low-LET NTE-TE IDER", "He4 Low-LET NTE-TE IDER", "O16 NTE-TE IDER", 
                               "Si28 NTE-TE IDER", "Ti48 NTE-TE IDER", "Fe56 (600 MeV/u) NTE-TE IDER", 
                               "IEA MIXDER (Equally Distributed)", "SEA MIXDER (Equally Distributed)"),
       col = c("orange", "yellow", "green", "blue", "purple", "violet", "red", "black"), 
       lwd = c(2, 2, 2, 2, 2, 2, 2, 2), 
       lty = c(1, 1, 1, 1, 1, 1, 1, 2), cex = 0.7, inset = 0.0125)
dev.off()


## 8-panel plot: 7 HZE IDERS, one low-LET HZE IDER in bottom right
par(mfrow = c(2, 4))
par(cex = 0.6)
par(mar = c(3, 3, 0, 0), oma = c(1, 1, 1, 1))
for (LET_value in c(25, 70, 100, 195, 464, 953)) {
  plot(x = hundred_cGy, y = calib_low_LET_ider(dose = hundred_cGy, L = LET_value), 
       type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l',
       col = "black", lwd = 2, lty = 2, xaxs = "i")
  
}
plot(x = hundred_cGy, y = calib_low_LET_ider(dose = hundred_cGy, L = 0.4), 
     type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l',
     col = "black", lwd = 2, lty = 2, xaxs = "i")
