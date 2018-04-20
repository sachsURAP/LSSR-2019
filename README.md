# NASAmouseHG  
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)  

Data and script for analyzing ionizing radiation tumorigenesis outside low-earth orbit.

## Repository Contents

#### `NASAmouseHG`
The main directory. Houses our primary R scripts for analysis, licensing information, git integration, and important subfolders.
The primary scripts, in order of development, are: `hgData.R`, `synergyTheory.R`, `monteCarlo.R`, and `plots.R`.

#### `data`
Data files, strictly in CSV format. 

#### `paper`
Paper drafts, supplementary materials, and writing materials. 

#### `exploration`
Scripts for indirectly related analyses, ideas, and experiments.

#### `misc_materials`
Meeting minutes, example scripts for various procedures, past agendas, and past papers.  

## Attribution
    # Copyright:    (C) 2017-2018 Sachs Undergraduate Research Apprentice Program
    #               This program and its accompanying materials are distributed 
    #               under the terms of the GNU General Public License v3.
    # Project:      NASAmouseHG
    # Purpose:      Data and script concerning synergy theory for murine Harderian
    #               gland tumorigenesis after irradiation by mixtures of ionized, 
    #               high-energy, atomic nuclei. 
    # Contact:      Rainer K. Sachs 
    # Website:      https://github.com/sachsURAP/NASAmouseHG
    # Mod history:  19 Apr 2018
    # Attribution:  This R script was developed at UC Berkeley. Written by Dae Woong 
    #               Ham Summer 2017. Additions, corrections, changes, quality 
    #               control, reorganization by Edward Huang, Yimin Lin, Mark Ebert,
    #               Yunzhi Zhang and Ray Sachs UCB 2017-2018 academic year.  

## License Notice
    #  NASAmouseHG is free software: you can redistribute it and/or modify
    #  it under the terms of the GNU General Public License as published by
    #  the Free Software Foundation, either version 3 of the License, or
    #  (at your option) any later version.
    # 
    #  NASAmouseHG is distributed in the hope that it will be useful,
    #  but WITHOUT ANY WARRANTY; without even the implied warranty of
    #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #  GNU General Public License for more details.
    # 
    #  You should have received a copy of the GNU General Public License
    #  along with NASAmouseHG. If not, see <http://www.gnu.org/licenses/>.

## Relevant References and Abbreviations:
    
    #   ".93Alp" = Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-
    #                           particle radiations." Rad Res 136:382-391 (1993)
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
    #   "HZE"    = High atomic number and energy
    #   "LET"    = Linear energy transfer
    #   "NTE"    = Non-targeted effects
    #   "TE"     = Targeted effects
    #   "IDER"   = Individual dose-effect relationship
    #   "MIXDER" = Mixed dose-effect relationship
    #   "SEA"    = Simple Effect Additivity
    #   "IEA"    = Incremental Effect Additivity
    #   "cGy"    = Centigray

