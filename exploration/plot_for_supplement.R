# Copyright:    (C) 2017-2018 Sachs Undergraduate Research Apprentice Program
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3.
# Filename:     plots.R 
# Purpose:      Concerns radiogenic mouse Harderian gland tumorigenesis. 
#               Contains code to generate figures. It is part of the 
#               source code for the NASAmouseHG project.
# Contact:      Rainer K. Sachs 
# Website:      https://github.com/sachsURAP/NASAmouseHG
# Mod history:  07 Jun 2018
# Details:      See hgData.R for further licensing, attribution, references, 
#               and abbreviation information.

source("monteCarlo.R") #  Load Monte Carlo

library(ggplot2) # Ribbon plot functionality
library(grid)  # Plot grids
library(Hmisc) # Error bars

#========= 9-panel plot for supplement: 8 HZE DER panels, 1 low-LET placeholder ============#
plot(1,1)
par(mfrow = c(3, 3))
par(cex = .5, tcl=-0.3, pch=19)
par(mar = c(2 , 0, 2, 5), oma = c(4, 4, 0.5, 0.5))
par(mgp = c(2, 0.6, 0))

plot(c(0,700),c(0,82),col="white", ann=FALSE) # low-LET
legend(x = "topleft", legend = "low LET swift light ions:
       see Fig. 5
       ", cex= 1, inset = 0.1)

Ne=ion_data[13:17,] # Ne 13:17
dmNe=ion_data[17,"dose"]
dNe=c(0.01*0:9, 0.1*1:9, 1:dmNe)
plot(dNe, calibrated_HZE_nte_der(dose=dNe,Ne[1,"LET"]),type='l',ylim=c(0,.82),  ann=FALSE, col='red', lwd=2 )
errbar(Ne[, "dose"], Ne[, "Prev"],yplus=Ne[, "Prev"]+1.96*Ne[, "SD"], yminus=Ne[, "Prev"]-1.96*Ne[, "SD"],
       cap=0.03, add=TRUE, col='black',  errbar.col = 'black', lwd=1)

Si=ion_data[18:21,] # Si 18:21
dmSi=ion_data[21,"dose"]
dSi=c(0.01*0:9, 0.1*1:9, 1:dmSi)
plot(dSi, calibrated_HZE_nte_der(dose=dSi,Si[1,"LET"]),type='l',ylim=c(0,.82),  ann=FALSE, col='red', lwd=2 )
errbar(Si[, "dose"], Si[, "Prev"],yplus=Si[, "Prev"]+1.96*Si[, "SD"], yminus=Si[, "Prev"]-1.96*Si[, "SD"],
       cap=0.03, add=TRUE, col='black',  errbar.col = 'black', lwd=1)

Ti=ion_data[22:26,] # Ti 22:26
dmTi=ion_data[26,"dose"]
dTi=c(0.01*0:9, 0.1*1:9, 1:dmTi)
plot(dTi, calibrated_HZE_nte_der(dose=dTi,Ti[1,"LET"]),type='l',ylim=c(0,.82),  ann=FALSE, col='red', lwd=2 )
errbar(Ti[, "dose"], Ti[, "Prev"],yplus=Ti[, "Prev"]+1.96*Ti[, "SD"], yminus=Ti[, "Prev"]-1.96*Ti[, "SD"],
       cap=0.03, add=TRUE, col='black',  errbar.col = 'black', lwd=1)

Fe1=ion_data[27:29,]
dmFe1=ion_data[29,"dose"]
dFe1=c(0.01*0:9, 0.1*1:9, 1:dmFe1)
plot(dFe1, calibrated_HZE_nte_der(dose=dFe1,Fe1[1,"LET"]),type='l',ylim=c(0,.82),  ann=FALSE, col='red', lwd=2 )
errbar(Fe1[, "dose"], Fe1[, "Prev"],yplus=Fe1[, "Prev"]+1.96*Fe1[, "SD"], yminus=Fe1[, "Prev"]-1.96*Fe1[, "SD"],
       cap=0.03, add=TRUE, col='black',  errbar.col = 'black', lwd=1)

Fe2=ion_data[30:36,]
dmFe2=ion_data[36,"dose"]
dFe2=c(0.01*0:9, 0.1*1:9, 1:dmFe2)
plot(dFe2, calibrated_HZE_nte_der(dose=dFe2,Fe2[1,"LET"]),type='l',ylim=c(0,.82),  ann=FALSE, col='red', lwd=2 )
errbar(Fe2[, "dose"], Fe2[, "Prev"],yplus=Fe2[, "Prev"]+1.96*Fe2[, "SD"], yminus=Fe2[, "Prev"]-1.96*Fe2[, "SD"],
       cap=0.03, add=TRUE, col='black',  errbar.col = 'black', lwd=1)

Fe3=ion_data[37:40,] # Fe3 37:40
dmFe3=ion_data[40,"dose"]
dFe3=c(0.01*0:9, 0.1*1:9, 1:dmFe3)
plot(dFe3, calibrated_HZE_nte_der(dose=dFe3,Fe3[1,"LET"]),type='l',ylim=c(0,.82),  ann=FALSE, col='red', lwd=2 )
errbar(Fe3[, "dose"], Fe3[, "Prev"],yplus=Fe3[, "Prev"]+1.96*Fe3[, "SD"], yminus=Fe3[, "Prev"]-1.96*Fe3[, "SD"],
       cap=0.03, add=TRUE, col='black',  errbar.col = 'black', lwd=1)

Nb=ion_data[41:44,] # Nb 41:44
dmNb=ion_data[44,"dose"]
dNb=c(0.01*0:9, 0.1*1:9, 1:dmNb)
plot(dNb, calibrated_HZE_nte_der(dose=dNb,Nb[1,"LET"]),type='l',ylim=c(0,.82),  ann=FALSE, col='red', lwd=2 )
errbar(Nb[, "dose"], Nb[, "Prev"],yplus=Nb[, "Prev"]+1.96*Nb[, "SD"], yminus=Nb[, "Prev"]-1.96*Nb[, "SD"],
       cap=0.03, add=TRUE, col='black',  errbar.col = 'black', lwd=1)

La=ion_data[45:47,] # La 45:47
dmLa=ion_data[47,"dose"]
dLa=c(0.01*0:9, 0.1*1:9, 1:dmLa)
plot(dLa, calibrated_HZE_nte_der(dose=dLa,La[1,"LET"]),type='l',ylim=c(0,.82),  ann=FALSE, col='red', lwd=2 )
errbar(La[, "dose"], La[, "Prev"],yplus=La[, "Prev"]+1.96*La[, "SD"], yminus=La[, "Prev"]-1.96*La[, "SD"],
       cap=0.03, add=TRUE, col='black',  errbar.col = 'black', lwd=1)

# dev.off()
