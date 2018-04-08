#   Filename: Data.R 
#   Purpose: Concerns radiogenic mouse HG tumorigenesis.

#   Copyright: (C) 2017 Mark Ebert, Edward Huang, Dae Woong Ham, Yimin Lin, and Ray Sachs 

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License version 3 as published 
#   by the Free Software Foundation.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See 
#   GNU GPLv3 for more details. <http://www.gnu.org/licenses/>.

#   Attribution Information: This R script was developed at UC Berkeley.
#   Written by Dae Woong Ham Summer 2017. Additions, corrections, changes, 
#   quality control, reorganization by Edward Huang, Yimin Lin, Mark Ebert. Yunzhi Zhang 
#   and Ray Sachs UCB 2017-2018 academic year

#   Relevant references and abbreviations:
#   ".93Alp" = Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-particle radiations." Rad Res 136:382-391 (1993)
#   ".94Alp" = Alpen et al. "Fluence-based relative biological effectiveness for charged particle carcinogenesis in mouse Harderian gland." Adv Space Res 14(10): 573-581. (1994).  
#   "16Chang" = Chang et al. "Harderian Gland Tumorigenesis: Low-Dose and LET Response." Radiat Res 185(5): 449-460. (2016).  
#   "16Srn" = Siranart et al."Mixed Beam Murine Harderian Gland Tumorigenesis: Predicted Dose-Effect Relationships if neither Synergism nor Antagonism Occurs." Radiat Res 186(6): 577-591 (2016).  
#   "17Cuc" = Cucinotta & Cacao. "Non-Targeted Effects Models Predict Significantly Higher Mars Mission Cancer Risk than Targeted Effects Models." Sci Rep 7(1): 1832. (2017). PMC5431989
rm(list=ls())
ion_data=data.frame(read.csv("hg_raw_data_ordered.csv", sep=",",header=TRUE)) ## Data used in 16Chang; includes data analyzed in .93Alp and .94Alp. Does not include gamma-ray data or ion background (i.e. control) data.
#These two lines of script plus the .csv file (which will need work as additions and perhaps corrections come up) should be all we need.
# However the following, which shows how to compute ion speed and the Katz track structure parameter, may be useful later. The following data.frame itself is obsolete.
# hg_data <- data.frame 
#   HG = c(0.091,0.045,0.101,0.169,0.347,0.431,0.667,0.623,0.156,0.215,0.232,0.307,0.325,0.554,0.649,0.123,0.145,0.207,0.31,0.026,0.083,0.25,0.39,0.438,0.424,0.093,0.195,0.302,0.292,0.109,0.054,0.066,0.128,0.286,0.183,0.167,0.396,0.536,0.192,0.234,0.317,0.092,0.131,0.124,0.297,0.082,0.088,0.146,0.236,0.371,.154,.132,.333), #  HG prevalence as defined in 16Chang
#   NWeight = c(520,2048,1145,584,313,232,293,221,1162,877,455,409,374,223,320,742,661,347,131,6081,1091,251,244,191,131,645,255,199,111,649,378,973,833,201,468,381,197,109,496,257,185,1902,1063,884,350,1767,1408,874,299,261,322,206,67), # weight for weighted least squaresregression; see .93Alp. The Lanthanum entries were obtained by measuring the main graph in 17Cuc 
#   index = c(rep(1, 8),rep(0,17), rep(1, 4),  rep(0, 24)), #  Index = 0 for Z > 3 ions, 1 otherwise. Not needed in some models
#   L = c(rep(1.6,8), rep(193, 7), rep(250, 4), rep(195, 6), rep(0.4, 4), rep(25, 5), rep(464, 4), rep(193, 3),rep(70, 4), rep(100, 5), rep(953, 3)), #  L = LET = LET_infinity = stopping power (keV/micron)
#   Z = c(rep(2, 8), rep(26, 17), rep(1, 4), rep(10, 5), rep(43, 4), rep(26, 3), rep(14, 4), rep(22, 5), rep(57, 3)), #  Atomic number, charge in units of proton charge on fully ionized atomic nucleus, e.g. 2 for 2He4
#   Zeff = c(rep("TBD", 53)), #  Effective ion charge according to the formula of W.H Barkas. Zeff <= Z. Calculated below. For this data, only very slightly less than Z.
#   beta = c(rep("TBD", 53)), #  Ion speed, relative to speed of light, calculated below
#   MeVperu = c(rep(228, 8), rep(600, 7), rep(300, 4), rep(600, 6), rep(250, 4), rep(670, 5), rep(600, 4), rep(600, 3), rep(260, 4), rep(1000, 5), rep(593, 3)), #  Kinetic energy in MeV, divided by atomic mass, e.g. divided by 4u = 4 x 931.5 MeV/c^2 for 2He4
#   Katz = c(rep("TBD", 53)), #  For fully ionized nuclei, Katz's Z^2/beta^2, Calculated below. It is part of the Bethe Barkas Bloch equation for stopping power. Our calculations don't use Katz, but various similar calculations do.
#   ion = c(rep("He4", 8), rep("Fe56", 17), rep("p", 4), rep("Ne20", 5), rep("Nb93", 4), rep("Fe56", 3), rep("Si28", 4), rep("Ti48", 5), rep("La139", 3)),
#   comments = c(".93AlpLooksOK", rep("", 7), ".93AlplooksOK", rep("", 11), ".93Alp.no.iso", "not in 17Cuc (or 16Chang?)", rep("", 3), "16Chang all OK?", rep('', 8),".93Alp",rep('',3),"16Chang",rep('',11), ".94Alp","From graphs",'e.g. in 17Cuc') # RKS to EH: I added two comments. .93 Alpen should now appear on the first Nb row (row 35) and 16 chang on the first 56Fe row below that (row 39)
# ) 
 
## Data for HG induced by photons from Cs-137 or Co-60 beta decay; from 16Chang 
# beta_decay_data <- data.frame(
#   dose = c(0, 0.4, 0.8, 1.6, 3.2, 7, 0, .4, .8, .12, 1.6),
#   HG = c(.026, .048, .093, .137, .322, .462, .0497, .054, .067, .128, .202),
#   NWeight = c(6081.2, 4989.5, 1896.8, 981.1, 522.2, 205.2, 7474.1, 2877.6, 1423.7, 689.9, 514.9),
#   Nucleus = c(rep("Cobalt-60", 6), rep("Cesium-137", 5)),
#   Comments = c(rep("TBD", 11))
# )
# 
# GeVu <- 0.001 * hg_data[, "MeVperu"] #  Convert to GeV/u for convenience in a calculation
# hg_data[, "Katz"] <- round(hg_data[, "Z"] ^ 2 * (2.57 * GeVu ^2 + 4.781 * GeVu + 2.233) / (2.57 * GeVu ^ 2 + 4.781 * GeVu), 2) #  Special relativistic calculation of Z^2/beta^2. The numerics include conversion from GeV to joules and from u to kg.
# hg_data[, "beta"] <- round(hg_data[, "Z"] * sqrt(1 / hg_data[, "Katz"]), 3) #  i.e. Z * sqrt(beta ^ 2 / Z ^ 2) 
# hg_data[, "Zeff"] <- round(hg_data[, "Z"] * (1 - exp( -125 * hg_data[, "Z"] ^ ( - 2.0 / 3))), 2) #  Barkas formula for Zeff; for us Zeff is almost Z
# 
# hg_data <- within(hg_data, L[L < 200 & ion == 'Fe56'] <- 195) # Set all Fe56 with L < 200 to L = 195 
# clean_hg_data <- hg_data[c(1:19, 21:53), ] # Removes the zero dose case 
# clean_HZE_data <- subset(clean_hg_data, Z > 3) #  Look only at HZE not at much lower Z and LET ions. 
# clean_light_ion_data <- subset(clean_hg_data, Z <= 3) 
# #  NOTE: "Light" refers to ionized atomic nuclei lighter than Beryllium. 
# #  The data published to date has such a big gap between alpha particles (Z=2, LET~1.6) 
# #  and Neon (Z=10, LET~25) that we will here use independent models for Z<3 and Z>3 