# 1. Preparing the GIS files and the different colour scales
# 2. Extracting the spatio-temporal information embedded in annotated trees
# 3. Investigating the profile of long-distance dispersal events
# 4. Investigating the patterns of isolation-by-distance
# 5. Estimating some dispersal statistics for each reconstruction
# 6. Visualising selected phylogenetic trees for SC1 and SC2
# 7. Visualising selected continuous phylogeographic reconstructions
# 8. Displaying and analysing the estimated position of the human ancestor
# 9. Visualising all continuous phylogeographic reconstructions
# 10. Mapping the bat species richness map on the topographic background

library(adephylo)
library(diagram)
library(HDInterval)
library(lubridate)
library(maptools)
library(rgdal)
library(seraphim)

source("Tree_data_extraction1.r")
source("Tree_data_extraction2.r")

writingFiles = TRUE; showingPlots = TRUE
writingFiles = FALSE; showingPlots = FALSE

analysis = "Analyses_2023-01-11"; clades = c("SC1","SC2") # analyses performed for the 1° submission
analysis = "Analyses_2023-10-18"; clades = c("SC2") # analyses performed with the clock rate of the 1° submission but including two ghost sequences
analysis = "Analyses_2023-10-19"; clades = c("SC2_noGhost","SC2_ghosts") # analyses based on the new rate priors based on the late 2020 data
analysis = "Analyses_2023-11-27"; clades = c("SC1_prior_stdevDiv5_2ghosts_combined_1-4",
											 "SC1_prior_stdevDiv5_Comb1-10",
											 "SC2_rateEarly2020_relrate_rootConst",
											 "SC2_rateEarly2020_relrate_rootConst_2ghosts",
											 "SC2_rateLate2020_relRate_rootConst",
											 "SC2_rateLate2020_relrate_rootConst_2ghosts")
analysis = "Analyses_2023-12-08"; clades = c("SC1","SC2") # new analyses for the visualisations (n.b.:)
											# for SC1, "timeTrees" correspond to "SC1_prior_stdevDiv5_Comb1-10", and
											# for SC2, "timeTrees" correspond to "SC2_rateLate2020_relRate_rootConst") 
burnIn = 0; selected_NRRs = c("NRR14","NRR3"); mostRecentSamplingDates = c(2021.534, 2021.285)

humanPangolinSeqs = c("AY394995.1", # human sequence (SC1)
					  "MT040333.1","MT040334.1","MT040335.1","MT040336.1","MT072864.1", # "pangolin" clade (SC1)
					  "MN908947.3", # human sequence (SC2)
					  "MT121216.1", # "isolated" pangolin sequence (SC2)
					  "MT040336.1","MT040334.1","MT040335.1","MT072864.1","MT040333.1", # "pangolin" clade (SC2)
					  "Ghost1","Ghost2") # "ghosts" sequences included in some tests

# 4. Investigating the patterns of isolation-by-distance

cutOff = 50; samplingCoordinates = read.csv("SC1_&_2_metadata.csv", head=T)[,c("name","longitude","latitude")]
directories = c("/SC1/noHumPanLoc/timeTrees/","/SC2/noHumPanLoc/timeTrees/"); rS_list_1 = list(); vS_list_1 = list()
for (h in 1:length(directories)) # to compute the correlation between patristic and geographic distances, as well as the ratio between geographic distances and the tMRCAs
	{
		files = list.files(paste0(analysis,directories[h]))
		files = files[which(grepl(".trees",files))]
		rS_list_2 = list(); all_rSs = c(); vS_list_2 = list(); all_vSs = c()
		for (i in 1:length(files))
			{
				trees = read.nexus(paste0(analysis,directories[h],files[i]))
				trees = trees[(length(trees)-99):(length(trees))]
				tempT = read.nexus(paste0(analysis,gsub("timeTrees","transPoW",directories[h]),files[i]))
				tips = trees[[1]]$tip.label; rSs = rep(NA, length(trees))
				for (j in 1:length(tempT))
					{
						# index = which(names(trees)==names(tempT)[j]); tree = trees[[index]]
						vSs = list(); tree = trees[[1]]
						tree = drop.tip(tree, paste0(clades[h],"_",humanPangolinSeqs))
						tab = read.csv(paste0(analysis,directories[h],gsub(".trees","_ext1",files[i]),"/TreeExtractions_",j,".csv"))
						distTree = as.matrix(distTips(tree, method="patristic")) # 2° method
						distsGeo = matrix(nrow=dim(distTree)[1], ncol=dim(distTree)[2])
						for (k in 2:dim(distsGeo)[1])
							{
								for (l in 1:(k-1))
									{
										index1 = which(tab[,"tipLabel"]==tree$tip.label[k])
										index2 = which(tab[,"tipLabel"]==tree$tip.label[l])
										x1 = cbind(tab[index1,"endLon"], tab[index1,"endLat"])
										x2 = cbind(tab[index2,"endLon"], tab[index2,"endLat"])
										index1 = which(samplingCoordinates[,"name"]==gsub("SC1_","",gsub("SC2_","",tree$tip.label[k])))
										index2 = which(samplingCoordinates[,"name"]==gsub("SC1_","",gsub("SC2_","",tree$tip.label[l])))
										x1 = cbind(samplingCoordinates[index1,"longitude"], samplingCoordinates[index1,"latitude"])
										x2 = cbind(samplingCoordinates[index2,"longitude"], samplingCoordinates[index2,"latitude"])
										distsGeo[k,l] = rdist.earth(x1, x2, miles=F, R=NULL)
										distsGeo[l,k] = distsGeo[k,l]
									}
							}
						tMRCAs = matrix(nrow=dim(distTree)[1], ncol=dim(distTree)[2])
						MRCAs = ape::mrca(tree); nodeHeights = nodeHeights(tree)
						nodeAges = max(nodeHeights)-nodeHeights
						for (k in 2:dim(tMRCAs)[1])
							{
								for (l in 1:(k-1))
									{
										index = which(tree$edge[,2]==MRCAs[l,k])
										if (length(index) == 1)
											{
												tMRCAs[k,l] = nodeAges[index,2]
											}	else		{
												index = which(tree$edge[,1]==MRCAs[l,k])[1]
												tMRCAs[k,l] = nodeAges[index,1]
											}
									}
							}
						distTree = distTree[lower.tri(distTree)]; distsGeo = distsGeo[lower.tri(distsGeo)]; tMRCAs = tMRCAs[lower.tri(tMRCAs)]
						distsGeo = distsGeo[which(tMRCAs<cutOff)]; distTree = distTree[which(tMRCAs<cutOff)]; tMRCAs = tMRCAs[which(tMRCAs<cutOff)]						
						rSs[j] = cor(log(distTree[lower.tri(distTree)]),log(distsGeo[lower.tri(distsGeo)]), method="spearman")
						rSs[j] = cor(distTree[lower.tri(distTree)],distsGeo[lower.tri(distsGeo)], method="spearman")
						vSs[[j]] = distsGeo[lower.tri(distsGeo)]/distTree[lower.tri(distTree)]
					}
				rS_list_2[[i]] = rSs; vS_list_2[[i]] = vSs
			}
		for (i in 1:length(rS_list_2)) all_rSs = c(all_rSs, rS_list_2[[i]])
		for (i in 1:length(vS_list_2))
			{
				for (j in 1:length(vS_list_2[[i]])) all_vSs = c(all_vSs, vS_list_2[[i]][[j]])
			}
		HPD_rS = round(HDInterval::hdi(all_rSs)[1:2],2); HPD_vS = round(HDInterval::hdi(all_vSs)[1:2],2)
		cat("\t",clades[h],": median rS = ",round(median(all_rSs),2),", 95% HPD = [",HPD_rS[1],"-",HPD_rS[2],"]","\n",sep="")
		cat("\t",clades[h],": median vS = ",round(median(all_vSs),2),", 95% HPD = [",HPD_vS[1],"-",HPD_vS[2],"]","\n",sep="")
		rS_list_1[[h]] = rS_list_2; vS_list_1[[h]] = vS_list_2
	}

cutOff = 75; samplingCoordinates = read.csv("SC1_&_2_metadata.csv", head=T)[,c("name","longitude","latitude")]
directories = c("/SC1/noHumPanLoc/timeTrees/","/SC2/noHumPanLoc/timeTrees/"); rS_list_1 = list(); vS_list_1 = list()
for (h in 1:length(directories)) # to compute the correlation between patristic and geographic distances, as well as the ratio between geographic distances and the tMRCAs
	{
		files = list.files(paste0(analysis,directories[h]))
		files = files[which(grepl(".trees",files))]
		rS_list_2 = list(); all_rSs = c(); vS_list_2 = list(); all_vSs = c()
		for (i in 1:length(files))
			{
				trees = read.nexus(paste0(analysis,directories[h],files[i]))
				trees = trees[(length(trees)-99):(length(trees))]
				tempT = read.nexus(paste0(analysis,gsub("timeTrees","transPoW",directories[h]),files[i]))
				tips = trees[[1]]$tip.label; rSs = rep(NA, length(trees))
				for (j in 1:length(tempT))
					{
						# index = which(names(trees)==names(tempT)[j]); tree = trees[[index]]
						vSs = list(); tree = trees[[1]]
						tree = drop.tip(tree, paste0(clades[h],"_",humanPangolinSeqs))
						tab = read.csv(paste0(analysis,directories[h],gsub(".trees","_ext1",files[i]),"/TreeExtractions_",j,".csv"))
						distTree = as.matrix(distTips(tree, method="patristic")) # 2° method
						distsGeo = matrix(nrow=dim(distTree)[1], ncol=dim(distTree)[2])
						for (k in 2:dim(distsGeo)[1])
							{
								for (l in 1:(k-1))
									{
										index1 = which(tab[,"tipLabel"]==tree$tip.label[k])
										index2 = which(tab[,"tipLabel"]==tree$tip.label[l])
										x1 = cbind(tab[index1,"endLon"], tab[index1,"endLat"])
										x2 = cbind(tab[index2,"endLon"], tab[index2,"endLat"])
										index1 = which(samplingCoordinates[,"name"]==gsub("SC1_","",gsub("SC2_","",tree$tip.label[k])))
										index2 = which(samplingCoordinates[,"name"]==gsub("SC1_","",gsub("SC2_","",tree$tip.label[l])))
										x1 = cbind(samplingCoordinates[index1,"longitude"], samplingCoordinates[index1,"latitude"])
										x2 = cbind(samplingCoordinates[index2,"longitude"], samplingCoordinates[index2,"latitude"])
										distsGeo[k,l] = rdist.earth(x1, x2, miles=F, R=NULL)
										distsGeo[l,k] = distsGeo[k,l]
									}
							}
						tMRCAs = matrix(nrow=dim(distTree)[1], ncol=dim(distTree)[2])
						MRCAs = ape::mrca(tree); nodeHeights = nodeHeights(tree)
						nodeAges = max(nodeHeights)-nodeHeights
						for (k in 2:dim(tMRCAs)[1])
							{
								for (l in 1:(k-1))
									{
										index = which(tree$edge[,2]==MRCAs[l,k])
										if (length(index) == 1)
											{
												tMRCAs[k,l] = nodeAges[index,2]
											}	else		{
												index = which(tree$edge[,1]==MRCAs[l,k])[1]
												tMRCAs[k,l] = nodeAges[index,1]
											}
									}
							}
						distTree = distTree[lower.tri(distTree)]; distsGeo = distsGeo[lower.tri(distsGeo)]; tMRCAs = tMRCAs[lower.tri(tMRCAs)]
						distsGeo = distsGeo[which(tMRCAs<cutOff)]; distTree = distTree[which(tMRCAs<cutOff)]; tMRCAs = tMRCAs[which(tMRCAs<cutOff)]						
						rSs[j] = cor(log(distTree[lower.tri(distTree)]),log(distsGeo[lower.tri(distsGeo)]), method="spearman")
						rSs[j] = cor(distTree[lower.tri(distTree)],distsGeo[lower.tri(distsGeo)], method="spearman")
						vSs[[j]] = distsGeo[lower.tri(distsGeo)]/distTree[lower.tri(distTree)]
					}
				rS_list_2[[i]] = rSs; vS_list_2[[i]] = vSs
			}
		for (i in 1:length(rS_list_2)) all_rSs = c(all_rSs, rS_list_2[[i]])
		for (i in 1:length(vS_list_2))
			{
				for (j in 1:length(vS_list_2[[i]])) all_vSs = c(all_vSs, vS_list_2[[i]][[j]])
			}
		HPD_rS = round(HDInterval::hdi(all_rSs)[1:2],2); HPD_vS = round(HDInterval::hdi(all_vSs)[1:2],2)
		cat("\t",clades[h],": median rS = ",round(median(all_rSs),2),", 95% HPD = [",HPD_rS[1],"-",HPD_rS[2],"]","\n",sep="")
		cat("\t",clades[h],": median vS = ",round(median(all_vSs),2),", 95% HPD = [",HPD_vS[1],"-",HPD_vS[2],"]","\n",sep="")
		rS_list_1[[h]] = rS_list_2; vS_list_1[[h]] = vS_list_2
	}

cutOff = 9999; samplingCoordinates = read.csv("SC1_&_2_metadata.csv", head=T)[,c("name","longitude","latitude")]
directories = c("/SC1/noHumPanLoc/timeTrees/","/SC2/noHumPanLoc/timeTrees/"); rS_list_1 = list(); vS_list_1 = list()
for (h in 1:length(directories)) # to compute the correlation between patristic and geographic distances, as well as the ratio between geographic distances and the tMRCAs
	{
		files = list.files(paste0(analysis,directories[h]))
		files = files[which(grepl(".trees",files))]
		rS_list_2 = list(); all_rSs = c(); vS_list_2 = list(); all_vSs = c()
		for (i in 1:length(files))
			{
				trees = read.nexus(paste0(analysis,directories[h],files[i]))
				trees = trees[(length(trees)-99):(length(trees))]
				tempT = read.nexus(paste0(analysis,gsub("timeTrees","transPoW",directories[h]),files[i]))
				tips = trees[[1]]$tip.label; rSs = rep(NA, length(trees))
				for (j in 1:length(tempT))
					{
						# index = which(names(trees)==names(tempT)[j]); tree = trees[[index]]
						vSs = list(); tree = trees[[1]]
						tree = drop.tip(tree, paste0(clades[h],"_",humanPangolinSeqs))
						tab = read.csv(paste0(analysis,directories[h],gsub(".trees","_ext1",files[i]),"/TreeExtractions_",j,".csv"))
						distTree = as.matrix(distTips(tree, method="patristic")) # 2° method
						distsGeo = matrix(nrow=dim(distTree)[1], ncol=dim(distTree)[2])
						for (k in 2:dim(distsGeo)[1])
							{
								for (l in 1:(k-1))
									{
										index1 = which(tab[,"tipLabel"]==tree$tip.label[k])
										index2 = which(tab[,"tipLabel"]==tree$tip.label[l])
										x1 = cbind(tab[index1,"endLon"], tab[index1,"endLat"])
										x2 = cbind(tab[index2,"endLon"], tab[index2,"endLat"])
										index1 = which(samplingCoordinates[,"name"]==gsub("SC1_","",gsub("SC2_","",tree$tip.label[k])))
										index2 = which(samplingCoordinates[,"name"]==gsub("SC1_","",gsub("SC2_","",tree$tip.label[l])))
										x1 = cbind(samplingCoordinates[index1,"longitude"], samplingCoordinates[index1,"latitude"])
										x2 = cbind(samplingCoordinates[index2,"longitude"], samplingCoordinates[index2,"latitude"])
										distsGeo[k,l] = rdist.earth(x1, x2, miles=F, R=NULL)
										distsGeo[l,k] = distsGeo[k,l]
									}
							}
						tMRCAs = matrix(nrow=dim(distTree)[1], ncol=dim(distTree)[2])
						MRCAs = ape::mrca(tree); nodeHeights = nodeHeights(tree)
						nodeAges = max(nodeHeights)-nodeHeights
						for (k in 2:dim(tMRCAs)[1])
							{
								for (l in 1:(k-1))
									{
										index = which(tree$edge[,2]==MRCAs[l,k])
										if (length(index) == 1)
											{
												tMRCAs[k,l] = nodeAges[index,2]
											}	else		{
												index = which(tree$edge[,1]==MRCAs[l,k])[1]
												tMRCAs[k,l] = nodeAges[index,1]
											}
									}
							}
						distTree = distTree[lower.tri(distTree)]; distsGeo = distsGeo[lower.tri(distsGeo)]; tMRCAs = tMRCAs[lower.tri(tMRCAs)]
						distsGeo = distsGeo[which(tMRCAs<cutOff)]; distTree = distTree[which(tMRCAs<cutOff)]; tMRCAs = tMRCAs[which(tMRCAs<cutOff)]						
						rSs[j] = cor(log(distTree[lower.tri(distTree)]),log(distsGeo[lower.tri(distsGeo)]), method="spearman")
						rSs[j] = cor(distTree[lower.tri(distTree)],distsGeo[lower.tri(distsGeo)], method="spearman")
						vSs[[j]] = distsGeo[lower.tri(distsGeo)]/distTree[lower.tri(distTree)]
					}
				rS_list_2[[i]] = rSs; vS_list_2[[i]] = vSs
			}
		for (i in 1:length(rS_list_2)) all_rSs = c(all_rSs, rS_list_2[[i]])
		for (i in 1:length(vS_list_2))
			{
				for (j in 1:length(vS_list_2[[i]])) all_vSs = c(all_vSs, vS_list_2[[i]][[j]])
			}
		HPD_rS = round(HDInterval::hdi(all_rSs)[1:2],2); HPD_vS = round(HDInterval::hdi(all_vSs)[1:2],2)
		cat("\t",clades[h],": median rS = ",round(median(all_rSs),2),", 95% HPD = [",HPD_rS[1],"-",HPD_rS[2],"]","\n",sep="")
		cat("\t",clades[h],": median vS = ",round(median(all_vSs),2),", 95% HPD = [",HPD_vS[1],"-",HPD_vS[2],"]","\n",sep="")
		rS_list_1[[h]] = rS_list_2; vS_list_1[[h]] = vS_list_2
	}

