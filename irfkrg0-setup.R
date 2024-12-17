## Setup code for IRF manuscript ##

# Script metadata --------------------------------------------------------------

# Created by: Nicholas Bussberg and Chunfeng Huang
# Date last updated: 2024-12-17
# Filename: irfkrg0-setup.R
# Purpose: preliminary setup code to run the rest of the IRF kriging code

# Code establishes directory structure and packages. This file can be skipped
  # if the directory structure is setup and packages are installed 


# Basic setup ------------------------------------------------------------------

rm(list=ls())
options(width = 80)


# Directory structure ----------------------------------------------------------

# 1. All R scripts should be placed in the same working directory 
# 2. Set working directory PRIOR to running the code
# 3. The data file "tltmonamg.2015_5.6" should be placed in a subfolder named
    # "Data". Subfolder should be at the same directory level as the R scripts.
# 4. Create an empty subfolder to store graphics using the following code:
dir.create("graphs")


# Packages to install ----------------------------------------------------------

install.packages(c("fields", "sp", "maps", "geoR", "akima", "MASS"))



