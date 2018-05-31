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

forty_cGy <- 0:41 #RKS to EGH 5.31.18: dose-vectors are now usually used only once each
#much better to put them with each plot than globally here. I did that for my plots but not for 
#your ggplots.
sixty_cGy <- 0:61
seventy_cGy <- 0:71
hundred_cGy <- 0:101
forty_nine_cGy <- 0:50

#===============================================================================
#=== Shape of DER for Fe 600 MeV/u. Paper Fig. 3 (was Fig. 2.2.4.1)  4/23/18 ===
#===============================================================================
d3A <- 0.01 * 0:16000 # RKS to EGH. d notation inconsistent with forty_cGy etc. above
prevalence <- calibrated_HZE_nte_der(d3A, 193) 
plot(d3A, prevalence, type = 'l', bty = 'u')
legend(x = "bottomright", 
       legend = "dose in centiGy; Fe 193 zoom in twice",
       cex = 0.6, inset = 0.025)

d3b <- 2 * 10^-5 * 0:1600 # zoom in by a factor of 10^4
prevalence <- calibrated_HZE_nte_der(d3b, 193)
plot(d3b, prevalence, type = 'l', bty='u')

d3c <- 10^-6 * 0:1600 # zoom in by another factor of 20
prevalence <- calibrated_HZE_nte_der(d3c, 193)
plot(d3c, prevalence, type = 'l', bty = 'u')

#======================================================================#
#== Low LET data, error bars, and DER. Paper Fig. 5 (was Fig. 3.1.1) ==#
#======================================================================
d5=0:701 
plot(c(0,701), c(-.02,1), pch=19, col='white', ann=FALSE, bty='u')
lines(d5, 1-exp(-coef(summary(low_LET_model, correlation = TRUE))[1]*d5),lwd=2)
errbar(ion_data[5:12, "dose"], ion_data[5:12, "Prev"],
       yplus=ion_data[5:12, "Prev"]+1.96*ion_data[5:12, "SD"],
       yminus=ion_data[5:12, "Prev"]-1.96*ion_data[5:12, "SD"], pch = 19,cap=0.02,
       add=TRUE, col='red',errbar.col = 'red', lwd=2) # alpha particle data
errbar(ion_data[1:4, "dose"], ion_data[1:4, "Prev"],
       yplus=ion_data[1:4, "Prev"]+1.96*ion_data[1:4, "SD"], 
       yminus=ion_data[1:4, "Prev"]-1.96*ion_data[1:4, "SD"], pch = 19,cap=0.02,
       add=TRUE, col='black',errbar.col = 'black',lwd=2) #  proton data
#legend(x = "topleft", legend = c("protons","4He"), col = c("red", "black"),
#pch = c(19,19), cex = 1, inset = 0.025)

#=========== one HZE ion, one low-LET. Minor Paper Fig. 6  ==================#
# We will always use NTE & TE rather than TE-only for HZE in the minor paper. 
d <- 1 * 0:302
r <- c(.2, .8) # dose proportions. RKS to EGH: at most one low LET, always last in r( ).
plot(x = d, y = calibrated_HZE_nte_der(dose = d, L = 193), type = "l", #  Now plot DERs 
     xlab = "dose", ylab = "HG", bty = 'u', col = 'green', lwd = 2) # Fe DER
lines(x = d, y = calibrated_low_LET_der(d, 0), col = 'orange', lwd = 2) #low LET DER 
lines(x = d, y = calculate_id(d, 193, r, lowLET = TRUE)[, 2], col = "red", lwd = 3) # I(d)
lines(x = d, y = calibrated_HZE_nte_der(dose = .2 * d, L = 193) + 
        calibrated_low_LET_der(.8 * d, 0), lty= 2, lwd=2) # SEA mixture baseline S(d)

r <- c(.8, .2) # Panel B, proportions reversed with low LET small  
plot(x = d, y = calibrated_HZE_nte_der(dose = d, L = 193), type = "l", 
     xlab = "dose", ylab = "HG", bty = 'u', col = 'green', lwd = 2) 
lines(x = d, y = calibrated_low_LET_der(d, 0), col = 'orange', lwd = 2) 
lines(x = d, y = calculate_id(d, 193, r, lowLET = TRUE)[, 2], col = "red", lwd = 3) # I(d)
lines(x = d, y = calibrated_HZE_nte_der(dose = .8 * d, L = 193) +
        calibrated_low_LET_der(.2 *d, 0), lty= 2, lwd=2) 

#=========== minor paper Fig. 7, 80% lowLET + 4 HZE ============#
d7 = c(0.01*0:9, 0.1*1:9, 1:41)
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
lines(x = d7, y = calculate_id(d7, c(0.4, 40,110, 180, 250),
                               c(.8, rep(.05,4)), model = "NTE",
                               lowLET = TRUE)[, 2], col = "red", lwd = 2) #I(d)
legend(x = "topleft", legend = c("Low-LET","L=40", "L=110", "L=180", "L=250", 
                                 "I(d)", "S(d)"),
       col = c("orange", "green","purple","blue", "aquamarine2", "red", "black"), 
       lwd = c(2, 2, 2, 2, 2, 3, 2), 
       lty = c(1, 1, 1, 1,  1, 1, 2), cex = 0.3, inset = 0.025)

#=== Fig. 8 was Fig. 3.2.3. DERs: Fe56 (600 MeV/u), Si28, IEA, SEA ======
d8=c(0.01*0:9, 0.1*1:9, 1:41)
plot(c(0,40),c(0, 0.33), col="white", bty='l', xlab = "Dose (Gy)", ylab = "HG")
lines(x = d8, y = calibrated_HZE_nte_der(dose = d8, L = 70), col = "cyan", lwd = 2)
lines(x = d8, y = calibrated_HZE_nte_der(dose = d8, L = 193), col = "orange", lwd = 2)
lines(x = d8, y = calculate_id(d8, c(70, 193), c(0.5 , 0.5), model = "NTE")[, 2],
      col = "red", lwd = 2) # I(d)
lines(x = d8, y = calculate_SEA(d8, c(70, 193), c(1/2, 1/2), n = 2), col = "black",
      lwd = 2, lty = 2)
legend(x = "topleft", legend = c("Fe56 (600 MeV/u)", "Si28", "IEA", "SEA"),
       col = c("orange", "cyan", "red", "black"), lwd = c(2, 2, 2, 2),
       lty = c(1, 1, 1, 2), cex = 0.6, inset = 0.05)

#===== Minor paper Fig. 9. All 7 HZE ions; each ion contributes 7 cGy =======
d9 = c(0.01*0:9, 0.1*1:9, 0.5*2:100)
plot(x = d9, y = calculate_SEA(d9, c(25, 70, 100, 193, 250, 464, 953), rep(1/7,7)), 
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
                                 "La139 NTE-TE IDER",
                                 "IEA MIXDER (Equally Distributed)", "SEA MIXDER (Equally Distributed)"),
       col = c("pink", "orange", "aquamarine2", "green", "blue", "purple", "violet", "red", "black"), 
       lwd = c(2, 2, 2, 2, 2, 2, 2, 2, 2), 
       lty = c(1, 1, 1, 1, 1, 1, 1, 1, 2), cex = 0.3, inset = 0.0125)

#==============================================================================
# RKS to EGH. Always put low LET last in vectors; never put in two different low LET (just add doses);  
# The script incorrectly interprets the second low LET component as having non-zero NTE effects.

#==============================================================================#
#====================== Confidence Interval Ribbon Plots ======================#
#==============================================================================#

#========= FIGURE 10 (was Fig. 3.2.3) =========#  
# Fe56 (600 MeV/u) and Si28 in equal proportions for a total of 40 cGy
#  Declare ratios and LET values for plot
ratios <- c(1/2, 1/2)
LET_vals <- c(193, 70)
d10 = c(0.01*0:9,0.1*1:9,1:40)
# We begin with the correlated plot # RKS correlated not calibrated similar corrections needed often
# #RKS to EGH 5/20/18: Change ggplot2 to agree with my plot version. Different colors are OK though
calib_ci_3.2.3 <- simulate_monte_carlo(200, d10, LET_vals, ratios, model = "NTE") 
#  Construct a data.frame for ease of use with ggplot2 if ggplot2 is used
ci_data <- data.frame(dose = d10,
                      #  Monte Carlo values
                      corrBottom = calib_ci_3.2.3$monte_carlo[1, ],
                      corrTop = calib_ci_3.2.3$monte_carlo[2, ], # 
                      
                      
                      #  DER values
                      fe_six = calibrated_HZE_nte_der(dose = d10, L = 193),
                      si = calibrated_HZE_nte_der(dose = d10, L = 70),
                    
                      
                      #  I(d)
                      i = calculate_id(d10, LET_vals, ratios,
                                       model = "NTE")[, 2]) 

#  We make the ribbon plot for correlated parameters
plot(c(0,41), c(0,.30), col="white", bty='L')
polygon(x=c(d10,rev(d10)),y=c(ci_data[,"corrTop"], rev(ci_data[,"corrBottom"])),
        xpd=-1,col="yellow",lwd=.4,border="orange") ## narrow CI
lines(ci_data[,"dose"],ci_data[,"si"],col='violet', lwd = 2) 
lines(ci_data[,"dose"],ci_data[,"i"],col='red', lwd = 3)
lines(ci_data[,"dose"],ci_data[,"fe_six"],col='green',type = 'l', lwd = 3, lty = 2)


#  RKS to EGH 5/20/18. we do not need the uncorrelated Monte Carlo which would be
#  uncorr_ci_3.2.3 <- simulate_monte_carlo(200, forty_cGy, LET_vals, ratios, model = "NTE", vcov = FALSE)
#  RKS to EGH 5/23/18: change ggplot2 below so it agrees in substance with the above plot; can keep colors as is
#  Plotting call with ggplot2
# ci_plot <- ggplot(data = ci_data, aes = fill) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5), legend.position="right") +
#   labs(title = "Confidence Intervals", x = "Dose (cGy)", y = "HG Prevalence (%)") +
#   
#  # Ribbon plot of both confidence intervals
#  # geom_ribbon(aes(dose, ymin = uncorrBottom, ymax = uncorrTop), fill = "aquamarine2", alpha = 1) +  #  Uncorrelated in pink
#   geom_ribbon(aes(dose, ymin = corrBottom, ymax = corrTop), fill = "yellow", alpha = .7) + #  Correlated in dull blue
#   
#   # scale_fill_discrete(name="Type",
#   # breaks=c("pink", "blue"),
#   # labels=c("Correlated Monte Carlo", "Uncorrelated Monte Carlo")) +
#   
#   # DER plots
#   geom_line(aes(dose, y = si),  col = "orange", size = 1) + #  iron in green
#   geom_line(aes(dose, y = fe_six),  col = "darkblue", size = 1) + #  silicon in dark green
#   
#   # I(d) plot
#   geom_line(aes(dose, y = i), col = "red", size = 1) #  I(d) in red
#   
# ci_plot # Print figure

#========= FIGURE 11 (was 3.2.4) =================#
#=============== Correlated vs Uncorrelated CI Overlay Plot ==================#

#  We use the MIXDER shown as Figure 9 in MS09. Consists of all 7 HZE ions in our 5/20/2018 data set.
#  Declare ratios and LET values for plot
ratios <- c(1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7)
LET_vals <- c(25, 70, 100, 193, 250, 464, 953)
d11 =c(0.1*0:9, 1:50)
# We begin with the correlated plot # RKS correlated not calibrated similar corrections needed often
calib_ci_3.2.4 <- simulate_monte_carlo(200, d11, LET_vals, ratios, model = "NTE") 
# ###RKS rename: calibrated goes to correlated and 3.2.4 goes to 11 and hundred_cGy goes to d11 
# We now calculate the uncorrelated Monte Carlo
uncorr_ci_3.2.4 <- simulate_monte_carlo(200, d11, LET_vals, ratios, model = "NTE", vcov = FALSE)
# 
# Construct a data.frame for ease of use with ggplot2 if ggplot2 is used
ci_data <- data.frame(dose = d11,
# #                       #  Monte Carlo values
                        corrBottom = calib_ci_3.2.4$monte_carlo[1, ],
                        corrTop = calib_ci_3.2.4$monte_carlo[2, ], 
                        uncorrTop = uncorr_ci_3.2.4$monte_carlo[1, ],
                        uncorrBottom = uncorr_ci_3.2.4$monte_carlo[2, ],
# RKS to EGH: This is higher than uncorrTop so somewhere the two got mixed up
# #                       #  DER values
                        ne = calibrated_HZE_nte_der(dose = d11, L = 25),
                        si = calibrated_HZE_nte_der(dose = d11, L = 70),
                        ti = calibrated_HZE_nte_der(dose = d11, L = 100),
                        fe_six = calibrated_HZE_nte_der(dose = d11, L = 193),
                        fe_three = calibrated_HZE_nte_der(dose = d11, L = 250),
                        nb = calibrated_HZE_nte_der(dose = d11, L = 464),
                        la = calibrated_HZE_nte_der(dose = d11, L = 953),
                        i = calculate_id(d11, LET_vals, ratios, model = "NTE")[, 2])  
# # 
# # #  Plotting call. RKS will use the plot below instead of ggplot2
plot(c(0,50),c(0,.40), col="white", bty='u')
polygon(x=c(d11,rev(d11)),y=c(ci_data[,"uncorrTop"],rev(ci_data[,"uncorrBottom"])),xpd=-1,
        col="orange",lwd=.5,border="orange") # wide CI
# # # RKS: "uncorrTop" must actually be uncorrBottom and vice-versa?
polygon(x=c(d11,rev(d11)),y=c(ci_data[,"corrTop"],rev(ci_data[,"corrBottom"])),xpd=-1,
        col="yellow",border="orange", lwd=.2) # narrow CI
# RKS to EGH solid yellow ribbon hides part of the blue ribbon without opacity commands.
# The opacity commands were giving me problems
lines(ci_data[,"dose"],ci_data[,"fe_three"],col='black',type = 'l', lwd = 2, bty='u', ann=FALSE) 
# This IDER is the highest; so plot it first
lines(ci_data[,"dose"],ci_data[,"ne"],col='black', lwd = 2) #the lowest IDER
lines(ci_data[,"dose"],ci_data[,"nb"],col='black', lwd = 2)# 
lines(ci_data[,"dose"],ci_data[,"la"],col='black', lwd = 2)# 
lines(ci_data[,"dose"],ci_data[,"si"],col='black', lwd = 2)
lines(ci_data[,"dose"],ci_data[,"ti"],col='black', lwd = 2)
lines(ci_data[,"dose"],ci_data[,"fe_six"],col='black', lwd = 2)
lines(ci_data[,"dose"],ci_data[,"i"],col='red', lwd=3, lty=2)

#  Plotting call # RKS to EGH 5/20/18. I think the next might give the same curves as plot() above. Does it?
# RKS to EGH 5/20/18. We never need Monte Carlo for a mixture of low and high LET in the minor paper.
# But for later I think a toggle will be needed in monte_carlo.R, with low LET the last component of r
# ci_plot <- ggplot(data = ci_data, aes = fill) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5), legend.position="right") +
#   labs(title = "Confidence Intervals", x = "Dose (cGy)", y = "HG Prevalence (%)") +
#   
#   #  Ribbon plot of both confidence intervals
#   geom_ribbon(aes(dose, ymin = uncorrBottom, ymax = uncorrTop), fill = "aquamarine2", alpha = 1) + #  Uncorrelated in pink
#   geom_ribbon(aes(dose, ymin = corrBottom, ymax = corrTop), fill = "yellow", alpha = .7) + #  Correlated in dull blue
#   
#   # scale_fill_discrete(name="Type",
#                       # breaks=c("pink", "blue"),
#                       # labels=c("Correlated Monte Carlo", "Uncorrelated Monte Carlo")) +
#   
#   #  DER plots
#   geom_line(aes(dose, y = ne), col = "blue", size = 1) + #  neon in yellow
#   geom_line(aes(dose, y = si),  col = "orange", size = 1) + #  silicon in orange
#   geom_line(aes(dose, y = ti),  col = "green", size = 2) + # titanium in green
#   geom_line(aes(dose, y = fe_six),  col = "purple", size = 1) + #  iron 600 in purple
#   geom_line(aes(dose, y = fe_three),  col = "violet", size = 1) + #  iron 300 in violet
#   geom_line(aes(dose, y = nb),  col = "darkcyan", size = 1) + #  niobium in dark cyan 
#   geom_line(aes(dose, y = la),  col = "darkorange", size = 1) + #  lanthanum in black 
#   
#   # I(d) plot
#   geom_line(aes(dose, y = i), col = "red", size = 1) #  I(d) in red
# 
# ci_plot # Print figure#

# Plot 3: four HZE; NTE: RKS to EGH 5/23/18. I deleted this and some other plots not needed in the minor paper.

#==========Figure for Chang's new proton data point 5/22/2018=============#
ddose=0:82 
plot(ddose, 1-exp(-coef(summary(low_LET_model, correlation = TRUE))[1]*ddose),
     type='l',lwd=2, ann=FALSE, xlim=c(0,80),ylim=c(-.02,.35), bty='u')
errbar(ion_data[1:2, "dose"], ion_data[1:2, "Prev"],yplus=ion_data[1:2, "Prev"]+1.96*ion_data[1:2, "SD"],
       yminus=ion_data[1:2, "Prev"]-1.96*ion_data[1:2, "SD"], pch = 19,cap=0.03, add=TRUE, col='orange',
       errbar.col = 'orange', ann=FALSE,lwd=2) #  RKS: proton data points
errbar(ion_data[5:7, "dose"], ion_data[5:7, "Prev"],yplus=ion_data[5:7, "Prev"]+1.96*ion_data[5:7, "SD"],
       yminus=ion_data[5:7, "Prev"]-1.96*ion_data[5:7, "SD"], pch = 19,cap=0.03, add=TRUE, col='black',
       errbar.col = 'black', ann=FALSE, lwd=2) #alpha particles
errbar(60, 0.081,yplus=.081+.09,yminus=.081-.09, pch = 19,cap=0.05, add=TRUE, col='red',errbar.col = 'red', lwd=3)

#========= 8-panel plot for supplement: 7 HZE IDER panels, 1 low-LET IDER ============#
## RKS to EGH: This gives a figure that looks probably correct, but is highly incomplete. Needs data points with errorbars.
## The par commands sometimes leave my terminal in the wrong graphics state so I
##  commented out the whole 8-panel plot until that has also been fixed.
# plot(1,1)
# par(mfrow = c(2, 4)) 
# # Maybe we should put all plot commands not needed for the minor paper, including even this figure
# # needed in the supplement, into a separate script on auxiliary figures?
# par(cex = 0.6)
# par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
# par(tcl = -0.25)
# par(mgp = c(2, 0.6, 0))
# LET_values <- c(25, 70, 100, 193, 250, 464, 953)
# IDER_names <- c("He4", "Ne20", "Si28", "Ti48", "Fe56 (600 MeV/u)", "Nb93", "La139")
# hundred_cGy = 0:101 # RKS to EGH: for some reason this repeat may be needed?
# i = 1
# while (i < length(LET_values) + 1) {
#   plot(x = hundred_cGy, y = calibrated_HZE_nte_der(dose = hundred_cGy, L = LET_values[i]),
#        xlim = c(0, 100), ylim = c(0, 1),
#        type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l',
#        col = "red", lwd = 2, axes = FALSE, ann = FALSE)
#   
#   if (i %in% c(5, 6, 7))
#     axis(1, col = "black", col.axis = "black", at = seq(0, 100, 10))
#   if (i %in% c(1, 5))
#     axis(2, col = "black", col.axis = "black", at = seq(0, 1, .1))
#   grid(col = "gray60")
#   box(col = "gray60", lwd = 1.5)
#   title(IDER_names[i], line = -2, cex = 0.5)
#   i = i + 1
# }
# 
# plot(x = hundred_cGy, y = calibrated_low_LET_der(dose = hundred_cGy, L = 0.4),
#      xlim = c(0, 100), ylim = c(0, 1),
#      type = "l", xlab = "Dose (cGy)", ylab = "HG Prevalence (%)", bty = 'l',
#      col = "red", lwd = 2, axes = FALSE, ann = FALSE)
# axis(1, col = "black", col.axis = "black", at = seq(0, 100, 10))
# box(col = "gray60", lwd = 1.5)
# grid(col = "gray60")
# title("H1", line = -2, cex = 0.5)
# 
# mtext("Dose (cGy)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,
#       col = "black")
# mtext("HG Prevalence (%)", side = 2, outer = TRUE, cex = 0.7, line = 2.2,
#       col = "black")
# par(mfrow = c(1, 1))
