### README for IRF Kriging (irfkrg) Manuscript Code ###

This file provides instructions on how to reproduce the results for the associated manuscript. 



## Manuscript Information ##

Manuscript title: Intrinsic Random Function Kriging on the Sphere

Authors: Nicholas W. Bussberg (nbussberg@elon.edu), Jacob Shields, Chunfeng Huang



## Directory Requirements ##

1. All R scripts, including functions, should be placed in the same working directory
2. Working directory should be set prior to running the code
3. The data file "tltmonamg.2015_5.6" should be placed in a subfolder named "Data". 
	- Subfolder should be at the directory level as the R scripts.
	- Other datasets (simulated and subsets of the tlt file) will be placed in this folder
4. An empty subfolder named "Graphs" should be added as well to store the graphics
	- The file "irfkrg0V3-setup.R" will perform this task. If preferred to setup manually, this line should be commented out in the setup.R script. 



## Software and Package Requirements ##

R version 4.0.5 was used to run the analyses. 

The packages "fields", "sp", "maptools", "geoR", "akima", and "MASS" are required. These packages are installed in the file "irfkrg0V3-setup.R". If preferred to install manually, this line should be commented out in the setup.R script. 

The software XQuartz is required for the spatial graphics packages that are run in our final script file (see below). 



## Data Information ##

For the simulation studies, all datasets are generated from the scripts. 

For the real data study, the data are global temperature anomaly data from the National Space Science and Technology Center (NSSTC, 2017). 
	- Dataset name: tltmonamg.2015_5.6
	- Only the year 2015 was used for this study
	- This data consists of satellite radiance measurements to determine global temperature anomalies in the lower troposphere for each month
	- The satellite readings are represented on a regular latitude-longitude grid, ranging from -88.75 to 88.75 degrees latitude and -178.75 to 178.75 degrees longitude 



## Scripts and Code Execution ##

# Function Files #

The following scripts (irfkrg abbreviation left off the filename for simplicity) contain functions and will be sourced by other scripts

	- fcn01V3-sphere-harmonics.R - relevant spherical harmonic functions
	- fcn02V3-legendre-polynomials.R - Legendre polynomial functions
	- fcn03V3-distance.R - function to calculate Great Circle distance and function to calculate distances for all location pairs
	- fcn04V3-icf.R - functions that simulate IRF processes
	- fcn05V3-MoM-estimator.R - the MoM estimator function
	- fcn06V3-order-approx.R - function that uses the MoM estimator for each order (or degree) for the data
	- fcn07V4-criteria.R - functions that generate the criterion, S(k)
	- fcn08V3-r-estimation.R - function to estimate the parameter for the ICF
	- fcn09V3-kriger.R - kriging functions 
	- fcn10V2-parse-data.R - function to reformat real dataset


# Executable Scripts #

The code will run start to finish in alpha-numeric order. Note that when each file (other than the setup and function files) is Sourced with Echo in R, the file will generate an txt file with the extension .log. These files simply capture the input and output using the sink() function in R. 

Here is the list of scripts (in run order) and what they do. 

	1. irfkrg0V3-setup.R - installs the necessary packages and creates a graph directory structure
		- Running this code is optional if the graph subfolder and packages are installed separately

	2. irfkrg1V3-irf-data.R - simulates the IRF processes on the sphere
		- WARNING: This script can take 3-4 hours to run (on personal laptop)
		- Y00, Y1n1, Y10, etc. are the real spherical harmonic functions (see irfkrg-fcn01). 
			These are written out in Section 2 of our manuscript.
		- The A matrix is defined in Shields (2017, p. 49). 
			It is needed to help define the p_nu(x) functions within H(x, y) (Section 3 in our manuscript). 
		- Cmat matrix is the simulated IRFk covariance matrix. 
			The CF2() and CF3() functions are the reproducing kernel H(x, y) from Section 3 (see irfkrg-fcn04). 
		- G (and thus Gmat) help determine the estimated kappa. 
			G(k, h_j) is defined in Section 2.1 in our manuscript. It is essentially the MoM estimator for the covariance function.
			G is created with the order.approx() (see irfkrg-fcn06) and G_Hat() functions (see irfkrg-fcn05).

	3. irfkrg2V4-criterion-graphs.R - creates the criterion graphs found in the article.

	4. irfkrg3V3-kriging-split.R - simulates IRF2 and IRF3 datasets and splits the data into train/test subsets

	5. irfkrg4V4-kriging-k2.R - compares the IRF kriging procedure vs. ordinary kriging on simulated IRF2 data

	6. irfkrg5V4-kriging-k3.R - compares the IRF kriging procedure vs. ordinary kriging on simulated IRF3 data 

	7. irfkrg6aV2-sim-hist-k2.R - simulates 1000 IRF2 datasets and estimates kappa with numerical procedure for each
		- kappa estimates are plotted in histogram
		- WARNING: This script can take 1-2 hours to run

	8. irfkrg6bV2-sim-hist-k3.R - simulates 1000 IRF3 datasets and estimates kappa with numerical procedure for each
		- kappa estimates are plotted in histogram
		- WARNING: This script can take 1-2 hours to run

	9. irfkrg7aV2-parse-data.R - reformats the "tltmonamg.2015_5.6" dataset from NSSTC and computes the pairwise distances for all locations

	10. irfkrg7bV2-est-kappa.R - estimates kappa for the NSSTC dataset for January and April 2015




## References ##

NSSTC (2017). Temperature anomaly data from 1978 to present. National Space Science and Technology Center (NSSTC), https://www.nsstc.uah.edu/data/msu/t2lt/.

Shields, J. (2017). Intrinsic Random Functions on Spheres: Theory, Methods, and Application. PhD thesis, Indiana University.







