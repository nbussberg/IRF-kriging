### Code that estimates degree of non-homogeneity for temperature data
##  Created by: Nicholas Bussberg
##  Date last updated: 2022-12-14
##  Filename: irfkrg7bV2-est-kappa.R
##  Purpose: estimates the degree of non-homogeneity, kappa, for global
##  temperature data


### Setup ###
rm(list=ls())
options(width = 80)

filename <-"irfkrg7bV2-est-kappa"

sink(paste(filename, ".log", sep = "")) # create log file

# R version
R.version.string


### 00 Load packages and functions

library(fields)   # needed for image.plot()
library(sp)       # needed for coordinates, SpatialPointsDataFrame
library(maptools) # needed for wrld_simpl map

library(geoR) # for as.geodata()
library(akima) # for interp()

source(file = "irfkrg-fcn01V3-sphere-harmonics.R")
source(file = "irfkrg-fcn02V3-legendre-polynomials.R")
source(file = "irfkrg-fcn03V3-distance.R")
source(file = "irfkrg-fcn04V3-icf.R")
source(file = "irfkrg-fcn05V3-MoM-estimator.R")
source(file = "irfkrg-fcn06V3-order-approx.R")
source(file = "irfkrg-fcn07V4-criteria.R")


###  01 Load data

Lmat <- readRDS("Data/irfkrg7aV2-parse-data-Lmat.Rds")

data(wrld_simpl)


# Set bin distances
H <- 10                   # Number Distances at Which to Estimate
eps <- 0.1                # Use to Bin Distances in Estimation
Hvec <- as.vector(seq(0,1,by = 1/H))[1:H]

P <- dim(Lmat)[[1]]

Dvec <- readRDS("Data/irfkrg7aV2-parse-data-Dvec.Rds")



# Load data

full.data <- readRDS(paste("Data/irfkrg7aV2-parse-data-2015.Rds", sep = ""))

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
    pdf(file = paste("Graphs/", filename, "-2015-",  month,
                     "Criteria-plot.pdf", sep = ""), width = 10, height = 7)
    
    par(las=1, mar=c(5,5,2,1)+0.1)
    plot(0:6, crit.raw, xlab = "k", ylab = "", cex = 1.2,
         cex.axis = 1.2, cex.lab = 1.2)
    if(month==1) mtext("S(k)", 2, line=1.2, cex = 1.2, las=2, at=0.95)
    if(month==2) mtext("S(k)", 2, line=1.2, cex = 1.2, las=2, at=1.3)
    
    dev.off()
    
    
    # Criterion plot for log(S(k))
    pdf(file = paste("Graphs/", filename, "-2015-",  month,
        "logCriteria-plot.pdf", sep = ""), width = 10, height = 7)
    
    par(las=1, mar=c(5,5,2,1)+0.1)
    plot(0:6, log(crit.raw), xlab = "k", ylab = "", cex = 1.2,
         cex.axis = 1.2, cex.lab = 1.2)
    mtext("log(S(k))", 2, line=0.5, cex = 1.2, las=2, at=-0.5)
    
    dev.off()
    
    
    loc.data <- cbind((Lmat[,2])*(180/pi), (Lmat[,1] - pi/2)*(180/pi),
                      month.data)

    val.int <- interp(loc.data[,1], loc.data[,2], month.data)
        
    X.lim <- max(max(month.data), abs(min(month.data)))
        
    pdf(file = paste("Graphs/", filename, "-2015-", month, ".pdf", 
                     sep = ""), width = 10, height = 7)
    par(las=1)
    par(mar=c(2.5,3,1,1) + 0.1)
        
    image.plot(val.int, cex.axis = 1.2, xlim = c(-180, 180), 
               ylim = c(-90, 90), zlim = c(-X.lim, X.lim))
        

    plot(wrld_simpl,border = "black",add = T)
    dev.off()
}






# close log
sink()






