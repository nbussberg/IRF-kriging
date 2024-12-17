### Setup code for IRF manuscript
##  Created by: Nicholas Bussberg
##  Date last updated: 2021-12-16
##  Filename: irfkrg0V3-setup.R
##  Purpose: preliminary setup code to run the rest of the IRF kriging code

##  Code establishes directory structure and packages. This file can be skipped
    ## as long as the directory structure is setup and packages are installed 

## Setup

rm(list=ls())
options(width = 80)


## Creating directory structure

# 1. All R scripts should be placed in the same working directory 
# 2. Set working directory PRIOR to running the code
# 3. The data file "tltmonamg.2015_5.6" should be placed in a subfolder named
    # "Data". Subfolder should be at the same directory level as the R scripts.
# 4. Create an empty subfolder to store graphics using the following code:
dir.create("Graphs")



## Packages to install

install.packages(c("fields", "sp", "maptools", "geoR", "akima", "MASS"))



