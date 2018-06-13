# Copyright:    (C) 2017-2018 Sachs Undergraduate Research Apprentice Program
#               This program and its accompanying materials are distributed 
#               under the terms of the GNU General Public License v3.
# Filename:     dataAndInfo.R 
# Purpose:      Concerns radiogenic mouse Harderian gland tumorigenesis. Loads 
#               ion and tumor prevalence data from CSV files. It is part of the 
#               source code for the NASAmouseHG project.
# Contact:      Rainer K. Sachs 
# Website:      https://github.com/sachsURAP/NASAmouseHG
# Mod history:  15 Apr 2018
# Attribution:  This R script was developed at UC Berkeley. Written by Dae Woong 
#               Ham Summer 2017. Additions, corrections, changes, quality 
#               control, reorganization by Edward Huang, Yimin Lin, Mark Ebert,
#               Yunzhi Zhang and Ray Sachs UCB 2017-2018 academic year.

# Relevant references and abbreviations:
#
#   ".93Alp" = Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-
#                           particle radiations." Rad Res 136:382-391 (1993).
#
#   ".94Alp" = Alpen et al. "Fluence-based relative biological effectiveness for
#                           charged particle carcinogenesis in mouse Harderian 
#                           gland." Adv Space Res 14(10): 573-581. (1994).  
#
#   "16Chang" = Chang et al. "Harderian Gland Tumorigenesis: Low-Dose and LET 
#                            Response." Radiat Res 185(5): 449-460. (2016). 
#
#   "16Srn" = Siranart et al. "Mixed Beam Murine Harderian Gland Tumorigenesis: 
#                             Predicted Dose-Effect Relationships if neither 
#                             Synergism nor Antagonism Occurs." 
#                             Radiat Res 186(6): 577-591 (2016).  
#
#   "17Cuc" = Cucinotta & Cacao. "Non-Targeted Effects Models Predict 
#                                Significantly Higher Mars Mission Cancer Risk 
#                                than Targeted Effects Models." 
#                                Sci Rep 7(1): 1832. (2017). PMC5431989
#
#   "HZE"     = High atomic number and energy
#   "LET"     = Linear energy transfer
#   "NTE"     = Non-targeted effects
#   "TE"      = Targeted effects
#   "DER"     = Dose-effect relation(ship)"
#   Obsolescent: "IDER" = one-ion DER; "MIXDER"  = Mixture baseline DER
#   "SEA"     = Simple Effect Additivity
#   "IEA"     = Incremental Effect Additivity
#   "cGy"     = Centigray

rm(list=ls()) # To be removed when script is finalized

# Data used in 16Chang; includes data analyzed in .93Alp and .94Alp. Does not 
# include gamma-ray data. Includes LET=100 keV/micron for Ti, an ad-hoc compromise
# between lower value at beam entry and higher value at mouse cage.

# These two lines of script plus the .csv file (which will need work as 
# additions and perhaps corrections come up) should be all we need.
ion_data <- data.frame(read.csv("data/raw_data_ordered.csv")) 

# However the older versions, which show how to compute ion speed and the Katz 
# track structure parameter, may be useful later.