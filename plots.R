#Filename: plots.R 
#   Purpose: Concerns radiogenic mouse HG tumorigenesis. Contains code to 
#            generate figures.

#   Copyright: (C) 2017 Mark Ebert, Edward Huang, Dae Woong Ham, Yimin Lin, and Ray Sachs 

source("hgData.R") #  load data
source("synergyTheory.R") #  load models
source("monteCarlo.R") #  load Monte Carlo

library(ggplot2) # ribbon plot functionality
library(grid)  # plot grids
library(Hmisc) #  error bars

forty_cGy <- .01 * 0:41 
# RKS->HG: I changed 0:40 to 0:41 better to use 0:41 or 0:42 to make sure the rightmost tick at 40 cGy appears on the graphs. Similarly for all other dose ranges
sixty_cGy <- .01 * 0:60
seventy_cGy <- .01 * 0:70
hundred_cGy <- .01 * 0:100
forty_nine_cGy <- .01 * 0:49

#par(mfrow = c(1, 1)) RKS->EH: I commented this out to get bigger graphs to look at.

#======================== PLOTS (2017) ==========================#
# Plot 1 : one HZE one low-LET; RKS->EH: always NTE & TE rather than TE-only in minor paper. For plots 1 and 2 I can and will do the key box specifying which curve is which by hand. When more than 2 ions are involved the key should be part of the script that generates the figure.
d <- .01 * 0:302. # RKS->EH use 0:302 to be safe in getting rightmost tick
r1 <- .2
r <- c(r1, 1 - r1) #Proportions. Next plot IDERs and MIXDER. RKS->EH changed L=195 and L=193 to L= 195 in Plots 1, 2 & 4
plot(x = d, y = calib_HZE_nte_ider(dose = d, L = 195), type = "l", xlab = "dose", ylab = "HG", bty = 'u', col = 'green', lwd = 2) #RKS->EH: often, as here the higher LET has more effect, so that using it for plot instead of lines is more likely to make a high y-axis tick conveniently appear.
lines(x = d, y = calib_low_LET_ider(d, 0), col = 'orange', lwd = 2) #RKS->EH changed color
lines(x = d, y = calculate_complex_id(r = r, L = 195, d = d, lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
lines(x=d, y=calib_HZE_nte_ider(dose =r1*d, L = 195)+calib_low_LET_ider((1-r1)*d, 0), lty= 3) # RKS->EH added SEA S(d)


# Plot 2 : one HZE one low-LET; NTE & TE rather than TE-only always in minor paper
d <- .01 * 0:302. # RKS->EH use 0:302 to be safe in getting rightmost tick
r1 <- .8 # RKS->EH this way all the rest is exactly like plot 1
r <- c(r1, 1 - r1) #Proportions. Next plot IDERs and MIXDER. RKS->EH changed L=195 and L=193 to L= 195 in Plots 1, 2 & 4
plot(x = d, y = calib_HZE_nte_ider(dose = d, L = 195), type = "l", xlab = "dose", ylab = "HG", bty = 'u', col = 'green', lwd = 2) #RKS->EH: often, as here the higher LET has more effect, so that using it for plot instead of lines is more likely to make a high y-axis tick conveniently appear.
lines(x = d, y = calib_low_LET_ider(d, 0), col = 'orange', lwd = 2) #RKS->EH changed color
lines(x = d, y = calculate_complex_id(r = r, L = 195, d = d, lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
lines(x=d, y=calib_HZE_nte_ider(dose =r1*d, L = 195)+calib_low_LET_ider((1-r1)*d, 0), lty= 3) # RKS->EH added SEA S(d)

# Plot 2 : one HZE one low-LET, TE & NTE. RKS->EH: this version of Plot 2 is now obsolete and I commented it out
# d <- .01 * 0:300
# r1 <- .2
# r <- c(r1, 1 - r1) #Proportions. Next plot IDERs and MIXDER
# plot(x = d, y = calib_HZE_te_ider(dose = d, L = 195), type = "l", xlab = "dose", ylab = "HG", bty = 'l', col = 'green', lwd = 2)
# lines(x = d, y = calib_low_LET_ider(d, 0), col = 'green', lwd = 2)
# lines(x = d, y = calculate_complex_id(r = r, L = 195, d = d, model = "TE", lowLET = TRUE)[, 2], col = "orange", lwd = 2) # I(d)
# lines(x = d, y = calculate_complex_id(r = r, L = 195, d = d, model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2)

# Plot 3: four HZE; NTE
plot(calculate_complex_id(r = rep(0.25, 4), L = c(25, 70, 190, 250), d = dose_vector), type = 'l', col = 'red', bty = 'l', ann = 'F') #  I(d) plot
SEA <- function(dose) {
  return(calib_HZE_nte_ider(dose / 4, 25) + 
           calib_HZE_nte_ider(dose / 4, 70) + 
           calib_HZE_nte_ider(dose / 4, 190) + 
           calib_HZE_nte_ider(dose / 4, 250))
}
lines(dose_vector, SEA(dose_vector), lty = 2)
lines(dose_vector, calib_HZE_nte_ider(dose_vector, 190), col = 'green') # component 4
lines(dose_vector, calib_HZE_nte_ider(dose_vector, 250), col = 'green') # component 3
lines(dose_vector, calib_HZE_nte_ider(dose_vector, 70), col = 'green') # component 2
lines(dose_vector, calib_HZE_nte_ider(dose_vector, 25), col = 'green') # component 1


# Plot 4: two HZE; NTE; one low-LET
d <- seq(0, .01, .0005)
plot(x = d, y = calib_HZE_nte_ider(dose = d, L = 195), type = "l", xlab = "dose", ylab = "HG", bty = 'l', col = 'green', lwd = 2)
lines(x = d, y = calib_HZE_nte_ider(d, 70), col = 'green', lwd = 2) # component 3
lines(x = d, y = calib_low_LET_ider(d, 0), col = 'green', lwd = 2)
lines(x = d, y = calculate_complex_id(r = c(1/20, 1/20, 9/10), L = c(70, 195), d = d, lowLET = TRUE)[, 2], col = 'red', lwd = 2)


#==================================================================#
#=========================== PAPER PLOTS ==========================#
#==================================================================#
# Fig. 3.2.3 - Fe56 (600 MeV/u), Si28, and corresponding IEA and SEA MIXDERS. RKS->EH changed Fig. number. It will change again later
plot(x = forty_cGy, y = calculate_SEA(forty_cGy, c(1/2, 1/2), c(70, 195), n = 2),bty='l', type = "l",  xlab = "Dose (Gy)", ylab = "HG Prevalence (%)", col = "black", lwd = 2, lty = 2) #RKS->EH removed "i", made boundary open at left
lines(x = forty_cGy, y = calib_HZE_nte_ider(dose = forty_cGy, L = 70), col = "cyan", lwd = 2)
lines(x = forty_cGy, y = calib_HZE_nte_ider(dose = forty_cGy, L = 195), col = "orange", lwd = 2) #RKS->EH changed color
lines(x = forty_cGy, y = calculate_complex_id(r = c(0.5 , 0.5), L = c(70, 195), d = forty_cGy, model = "NTE", lowLET = FALSE)[, 2], col = "red", lwd = 2) # I(d)
legend(x = "topleft", legend = c("Fe56 (600 MeV/u), NTE NTE-TE IDER", "Si28 HZE NTE-TE IDER", "IEA MIXDER (50% Fe56, 50% Si28)", "SEA MIXDER (50% Fe56, 50% Si28)"),
       col = c("orange", "cyan", "red", "black"), lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 2), inset = 0.05)


# Fig. for a mixture of 60 cGy H1 (protons) with 40 cGy Si;
plot(x = hundred_cGy, y = calib_HZE_nte_ider(dose = hundred_cGy, L = 70), type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "darkcyan", lwd = 2, xaxs="i")
lines(x = hundred_cGy, y = calib_low_LET_ider(dose = hundred_cGy, L = 0.4), col = "cyan", lwd = 2)
lines(x = hundred_cGy, y = calculate_SEA(hundred_cGy, c(0.6, 0.4), c(0.4, 70), lowLET = TRUE, n = 2), col = "black", lwd = 2, lty = 2)
lines(x = hundred_cGy, y = calculate_complex_id(r = c(0.6 , 0.4), L = c(0.4, 70), d = hundred_cGy, model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 1, lwd = 2)
legend(x = "topleft", legend = c("Si28 HZE NTE-TE IDER", "H1 Low-LET NTE-TE IDER","IEA MIXDER (60% H1, 40% Si28)", "SEA MIXDER (60% H1, 40% Si28)"),
       col = c("darkcyan", "cyan", "red", "black"), lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 2), inset = 0.025)

# Fig. for a mixture of 40 cGy H1 with 30 cGy Fe56 at 600 MeV/u;
plot(x = seventy_cGy, y = calib_HZE_nte_ider(dose = seventy_cGy, L = 195), col = "darkcyan", type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', lwd = 2, lty = 1, xaxs="i")
lines(x = seventy_cGy, y = calib_low_LET_ider(dose = seventy_cGy, L = 0.4), col = "cyan", lwd = 2)
lines(x = seventy_cGy, y = calculate_SEA(seventy_cGy, c(4/7, 3/7), c(0.4, 195), lowLET = TRUE, n = 2), col = "black", lty = 2, lwd = 2)
lines(x = seventy_cGy, y = calculate_complex_id(r = c(4/7 , 3/7), L = c(0.4, 195), d = seventy_cGy, model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
abline(v = 0.7, lwd = 2)
legend(x = "topleft", legend = c("Fe56 (600 MeV/u), HZE NTE-TE IDER", "H1 Low-LET NTE-TE IDER","IEA MIXDER (57% H1, 43% Si28)", "SEA MIXDER (57% H1, 43% Si28)"),
       col = c("darkcyan", "cyan", "red", "black"), lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 2), inset = 0.005)


# Fig. for a mixture of all 7 HZE ions; total dose 49 cGy, each ion gets 7 cGy;
# assume Hi HZE implies Z > 3
plot(x = forty_nine_cGy, y = calculate_SEA(forty_nine_cGy, c(1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7), c(25, 70, 100, 195, 250, 464, 953)), 
     type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2, axes = FALSE)

lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 25), col = "pink", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 70), col = "orange", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 100), col = "yellow", lwd = 3.5)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 195), col = "green", lwd = 3)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 250), col = "blue", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 464), col = "purple", lwd = 2, lty = 2)
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


#Fig. 3.2.2 RKS->EH added a new Fig.
# h1, 0.4, 80%
# L=c(40, 110, 180, 250), proportions (5,5,5,5)%
plot(x = forty_cGy, y = calculate_SEA(forty_cGy, c(.8, rep(.05,4)), c(0.4, 40,110, 180, 250), lowLET = TRUE), 
     type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2)

lines(x = forty_cGy, y = calib_low_LET_ider(dose = forty_cGy, L = 0.4), col = "orange", lwd = 2)
lines(x = forty_cGy, y = calib_HZE_nte_ider(dose = forty_cGy, L = 40), col = "green", lwd = 2)
lines(x = forty_cGy, y = calib_HZE_nte_ider(dose = forty_cGy, L = 110), col = "purple", lwd = 2)
lines(x = forty_cGy, y = calib_HZE_nte_ider(dose = forty_cGy, L = 180), col = "blue", lwd = 2)
lines(x = forty_cGy, y = calib_HZE_nte_ider(dose = forty_cGy, L = 250), col = "aquamarine2", lwd = 2)

lines(x = forty_cGy, y = calculate_complex_id(r = c(.8, rep(.05,4)), L =  c(0.4, 40,110, 180, 250),d = forty_cGy, model = "NTE", lowLET = TRUE)[, 2], col = "red", lwd = 2) # I(d)
legend(x = "right", legend = c("Low-LET","L=40", "L=110", "L=180", "L=250", 
                               "I(d)", "S(d)"),
       col = c("orange", "green","purple","blue", "aquamarine2", "red", "black"), 
       lwd = c(2, 2, 2, 2, 2, 3, 2), 
       lty = c(1, 1, 1, 1,  1, 1, 2))


# h1, 0.4, 60
# he4, 1.6, 20
# o16, 25, 10
# si28, 70, 2.5
# ti28, 100, 2.5
# fe56, 195, 5
plot(x = forty_nine_cGy, y = calculate_SEA(forty_nine_cGy, c(.6, .2, .1, .025, .025, .05), c(0.4, 1.6, 25, 70, 100, 195), lowLET = TRUE), 
     type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l', col = "black", lwd = 2, lty = 2, xaxs = "i")

lines(x = forty_nine_cGy, y = calib_low_LET_ider(dose = forty_nine_cGy, L = 0.4), col = "orange", lwd = 2)
lines(x = forty_nine_cGy, y = calib_HZE_nte_ider(dose = forty_nine_cGy, L = 1.6), col = "yellow", lwd = 2)
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
       col = c("orange", "yellow", "blue", "blue", "purple", "violet", "red", "black"), 
       lwd = c(2, 2, 2, 2, 2, 2, 2, 2), 
       lty = c(1, 1, 1, 1, 1, 1, 1, 2), cex = 0.7, inset = 0.0125)



# 8-panel plot: 7 HZE IDERS, one low-LET HZE IDER in bottom right
par(mfrow = c(2, 4))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
LET_values <- c(25, 70, 100, 195, 250, 464, 953)
IDER_names <- c("He4", "Ne20", "Si28", "Ti48", "Fe56 (600 MeV/u)", "Nb93", "La139")
i = 1
while (i < length(LET_values) + 1) {
  plot(x = hundred_cGy, y = calib_HZE_nte_ider(dose = hundred_cGy, L = LET_values[i]),
       xlim = c(0, 1), ylim = c(0, .7),
       type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l',
       col = "red", lwd = 2, axes = FALSE, ann = FALSE)
  
  if (i %in% c(5, 6, 7))
    axis(1, col = "black", col.axis = "black", at = seq(0, 1, .1))
  if (i %in% c(1, 5))
    axis(2, col = "black", col.axis = "black", at = seq(0, 0.6, .1))
  grid(col = "gray60")
  box(col = "gray60", lwd = 1.5)
  title(IDER_names[i], line = -2, cex = 0.5)
  i = i + 1
}

plot(x = hundred_cGy, y = calib_low_LET_ider(dose = hundred_cGy, L = 0.4),
     xlim = c(0, 1), ylim = c(0, 0.7),
     type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l',
     col = "red", lwd = 2, axes = FALSE, ann = FALSE)
axis(1, col = "black", col.axis = "black", at = seq(0, 1, .1))
box(col = "gray60", lwd = 1.5)
grid(col = "gray60")
title("H1", line = -2, cex = 0.5)

mtext("Dose (cGy)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,
      col = "black")
mtext("HG Prevalence (%)", side = 2, outer = TRUE, cex = 0.7, line = 2.2,
      col = "black")
# RKS->EH Some lines to see how errbar works
# ?errbar()
# plot(1:3,2:4,type='l')
# errbar(1:3,2:4,yminus=.5*1:3,yplus=c(2:4)+.5,cap=.05, add=TRUE)


##==================== Confidence Interval Ribbon Plots ======================#
ci_data <- data.frame(dose = 1:137,
                      monteCarloBottom = monteCarloCI[1, ],
                      monteCarloTop = monteCarloCI[2, ],
                      naiveBottom = naiveCI[1, ],
                      naiveTop = naiveCI[2, ])
ci_plot <- ggplot(data = ci_data, aes = fill) +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5), legend.position="right") + 
      labs(title = "Confidence Intervals", x = "Dose (cGy)", y = "HG Prevalence (%)") + 
      geom_ribbon(aes(dose, ymin = naiveBottom, ymax = naiveTop, fill = "blue")) +
      geom_ribbon(aes(dose, ymin = monteCarloBottom, ymax = monteCarloTop, fill = "pink")) + 
      scale_fill_discrete(name="Type",
                          breaks=c("blue", "pink"),
                          labels=c("Naive", "Monte Carlo")) 
ci_plot
