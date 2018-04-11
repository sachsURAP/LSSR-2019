data=read.csv("errorbar_example.csv")
library(Hmisc) #  error bars
ddose=0:700 #in cGy unfortunately
plot(c(0,7),c(0,1))  
errbar(data[9:12, "dose"]/100, data[9:12, "Prev"],yplus=data[9:12, "Prev"]+data[9:12, "SD"],yminus=data[9:12, "Prev"]-data[9:12, "SD"],pch=1, add=TRUE) #  RKS: proton data points
errbar(data[1:8, "dose"]/100, data[1:8, "Prev"],yplus=data[1:8, "Prev"]+data[1:8, "SD"],yminus=data[1:8, "Prev"]-data[1:8, "SD"], pch = 19, add=TRUE) #  RKS: Helium data points