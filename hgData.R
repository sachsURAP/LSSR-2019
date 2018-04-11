#   Filename: Data.R 
#   Purpose: Concerns radiogenic mouse HG tumorigenesis.

#   Copyright: (C) 2017 Mark Ebert, Edward Huang, Dae Woong Ham, Yimin Lin,
#                       Yunzhi Zhang, and Ray Sachs.

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License version 3 as published 
#   by the Free Software Foundation.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
#   GNU GPLv3 for more details. <http://www.gnu.org/licenses/>.

#   Attribution Information: This R script was developed at UC Berkeley.
#   Written by Dae Woong Ham Summer 2017. Additions, corrections, changes, 
#   quality control, reorganization by Edward Huang, Yimin Lin, Mark Ebert,
#   Yunzhi Zhang and Ray Sachs UCB 2017-2018 academic year.

#   Relevant references and abbreviations:
#   ".93Alp" = Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-particle radiations." Rad Res 136:382-391 (1993)
#   ".94Alp" = Alpen et al. "Fluence-based relative biological effectiveness for charged particle carcinogenesis in mouse Harderian gland." Adv Space Res 14(10): 573-581. (1994).  
#   "16Chang" = Chang et al. "Harderian Gland Tumorigenesis: Low-Dose and LET Response." Radiat Res 185(5): 449-460. (2016).  
#   "16Srn" = Siranart et al."Mixed Beam Murine Harderian Gland Tumorigenesis: Predicted Dose-Effect Relationships if neither Synergism nor Antagonism Occurs." Radiat Res 186(6): 577-591 (2016).  
#   "17Cuc" = Cucinotta & Cacao. "Non-Targeted Effects Models Predict Significantly Higher Mars Mission Cancer Risk than Targeted Effects Models." Sci Rep 7(1): 1832. (2017). PMC5431989

rm(list=ls()) #   To be removed when script is finalized
ion_data <- data.frame(read.csv("data/raw_data_ordered.csv", 
                                sep = ",", 
                                header = TRUE)) 
#  Data used in 16Chang; includes data analyzed in .93Alp and .94Alp. Does not 
#  include gamma-ray data or ion background (i.e. control) data.
#  These two lines of script plus the .csv file (which will need work as 
#  additions and perhaps corrections come up) should be all we need.
#   However the following, which shows how to compute ion speed and the Katz 
#  track structure parameter, may be useful later.
