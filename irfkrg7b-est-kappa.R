### Code that estimates degree of non-homogeneity for temperature data
##  Created by: Nicholas Bussberg
##  Date last updated: 2024-12-16
##  Filename: irfkrg7b-est-kappa.R
##  Purpose: estimates the degree of non-homogeneity, kappa, for global
##  temperature data


### Setup ###
rm(list=ls())
options(width = 80)

filename <-"irfkrg7b-est-kappa"

sink(paste(filename, ".log", sep = "")) # create log file

# R version
R.version.string


### 00 Load packages and functions

library(fields)   # needed for image.plot()
library(sp)       # needed for coordinates, SpatialPointsDataFrame
library(maps)     # needed for world map

library(geoR) # for as.geodata()
library(akima) # for interp()

source(file = "irf-fcn01-sphere-harmonics.R")
source(file = "irf-fcn02-legendre-polynomials.R")
source(file = "irf-fcn03-distance.R")
source(file = "irf-fcn04-covariances.R")
source(file = "irf-fcn05-MoM-estimator.R")
source(file = "irf-fcn06-kappa-approx.R")


###  01 Load data

Lmat <- readRDS("data/irfkrg7a-parse-data-Lmat.Rds")

data(wrld_simpl)


# Set bin distances
H <- 10                   # Number Distances at Which to Estimate
eps <- 0.1                # Use to Bin Distances in Estimation
Hvec <- as.vector(seq(0,1,by = 1/H))[1:H]

P <- dim(Lmat)[[1]]

Dvec <- readRDS("data/irfkrg7a-parse-data-Dvec.Rds")



# Load data

full.data <- readRDS(paste("data/irfkrg7a-parse-data-2015.Rds", sep = ""))

# Select Jan and Apr from 2015 for this study

JanApr.data <- full.data[,c(3,6)]



### 02 Estimate kappa for each month, plot criterion plots

kappa.vec <- rep(NA, 2)

for(month in 1:2){
    month.data <- JanApr.data[,month]
    if(month==1) month.lab <- "Jan"
    if(month==2) month.lab <- "Apr"
    
    GN <- order.approx(data = month.data, Lmat, Dvec, Hvec, eps, P)

    G.diff <- G.diff.fcn(GN$G, H, Hvec)

    crit.raw <- criterion.v01(G.diff)
    
    value.max <- which.max(crit.raw)
    
    crit.2 <- crit.raw+0.5*(0:6)
    crit.3 <- which.min(crit.2[value.max:7])+value.max-1 
    
    print(paste("Kappa-hat for ", month.lab, ": ", crit.3-1, sep=""))
    
    
    # Criterion plot for raw S(k)
    pdf(file = paste("graphs/", filename, "-2015-",  month,
                     "Criteria-plot.pdf", sep = ""), width = 10, height = 7)
    
    par(las=1, mar=c(5,6,2,1)+0.1)
    plot(0:6, crit.raw, xlab = "k", ylab = "", cex = 1.5,
         cex.axis = 1.5, cex.lab = 1.5, pch = 16)
    if(month==1) mtext("S(k)", 2, line=1.2, cex = 1.5, las=2, at=0.95)
    if(month==2) mtext("S(k)", 2, line=1.2, cex = 1.5, las=2, at=1.3)
    
    dev.off()
    
    
    # Criterion plot for log(S(k))
    pdf(file = paste("graphs/", filename, "-2015-",  month,
        "logCriteria-plot.pdf", sep = ""), width = 10, height = 7)
    
    par(las=1, mar=c(5,6,2,1)+0.1)
    plot(0:6, log(crit.raw), xlab = "k", ylab = "", cex = 1.5,
         cex.axis = 1.5, cex.lab = 1.5, pch = 16)
    mtext("log(S(k))", 2, line=0.5, cex = 1.5, las=2, at=-0.5)
    
    dev.off()
    
    
    loc.data <- cbind((Lmat[,2])*(180/pi), (Lmat[,1] - pi/2)*(180/pi),
                      month.data)

    val.int <- akima::interp(loc.data[,1], loc.data[,2], month.data)
        
    X.lim <- max(max(month.data), abs(min(month.data)))
        
    pdf(file = paste("graphs/", filename, "-2015-", month, ".pdf", 
                     sep = ""), width = 10, height = 7)
    par(las=1)
    par(mar=c(2.5,3,1,1) + 0.1)
        
    fields::image.plot(val.int, cex.axis = 1.2, 
                       xlim = c(-180, 180), 
                       ylim = c(-90, 90), 
                       zlim = c(-X.lim, X.lim))
        
    maps::map("world", add = T)
    dev.off()
}






# close log
sink()






