library(adephylo)
library(diagram)
library(HDInterval)
library(lubridate)
library(maptools)
library(rgdal)
library(seraphim)

source("rS_geo_pat_dist.r")

localTreeDirectories = c("WNV_gamma_all",
						 "GTEV_RRW_125",
						 "LASV_L_align_3",
						 "LASV_S_align_3",
						 "PDCV_Dlat_pols",
						 "RABV_US_skunk",
						 "RABV_PE_L1-L3")

# 1. Exploring the patterns of isolation-by-distance (IBD)

for (i in 1:length(localTreeDirectories)) # for (i in c(9:10,32:41))
	{
		localTreesDirectory = localTreeDirectories[i]; extractionFiles = list.files(localTreeDirectories[i])
		nberOfExtractionFiles = length(extractionFiles[which(grepl("TreeExtraction",extractionFiles))])
		rSs = rS_geo_pat_dist(localTreesDirectory, nberOfExtractionFiles)
		median = round(median(rSs),2); HPD = round(HDInterval::hdi(rSs)[1:2],2)
		cat("\t",localTreeDirectories[i])
		cat(": median rS = ",median,", 95% HPD = [",HPD[1],"-",HPD[2],"]","\n",sep="")
	}
		# WNV_gamma_all:  median rS = 0.00, 95% HPD = [-0.02-0.03]
		# GTEV_RRW_125:   median rS = 0.06, 95% HPD = [0.00-0.14]
		# LASV_L_align_3: median rS = 0.70, 95% HPD = [0.70-0.70]
		# LASV_S_align_3: median rS = 0.68, 95% HPD = [0.64-0.72]
		# PDCV_Dlat_pols: median rS = -0.02, 95% HPD = [-0.07-0.05]
		# RABV_US_skunk:  median rS = 0.49, 95% HPD = [0.48-0.60]
		# RABV_PE_L1-L3:  median rS = 0.73, 95% HPD = [0.55-0.83]

# 2. Estimating the weighted diffusion coefficients

for (i in 1:length(localTreeDirectories))
	{
		localTreesDirectory = localTreeDirectories[i]; extractionFiles = list.files(localTreeDirectories[i])
		nberOfExtractionFiles = length(extractionFiles[which(grepl("TreeExtraction",extractionFiles))])
		timeSlices = 100; onlyTipBranches = F; showingPlots = F
		outputName = paste0("Dispersal_stats/",localTreesDirectory); nberOfCores = 10; slidingWindow = 1
		spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
	}

