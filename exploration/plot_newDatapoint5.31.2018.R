# Copyright:    (C) 2017-2018 Sachs Undergraduate Research Apprentice Program
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3.
# Filename:     plot_newDatapoint5.31.2018.R
# Purpose:      Run this file instead of running plots.R if need this non-core plot
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


