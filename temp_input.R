#Takes chang data June 2018 on 2 ion HZE-lowerLET mixtures and computes prevalence, SD and 95% CI
df_temp=data.frame("mice"=c(40,44,40), "tumors" = c(19,13,13), row.names = c("fe","si","both"), "dose_HZE"= c(30,40,20), "dose_lowerLET"=c(40,60,20))
df_temp[, "prev"] =round(df_temp[,"tumors"]/df_temp[,"mice"],3)
df_temp [,"SD"] = round(sqrt(df_temp[,"prev"]*(1-df_temp[, "prev"])/df_temp[,"mice"]),4) # Ainsworth formula in Alp93
df_temp [,"95CI"] = round(2*1.96*df_temp [,"SD"],3)
df_temp[ , "mixdose"]= round(df_temp[,"dose_lowerLET"]+df_temp[,"dose_HZE"],3)
df_temp [,"rLow"]=round(df_temp [,"dose_lowerLET"]/df_temp [,"mixdose"],3)