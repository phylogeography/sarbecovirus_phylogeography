spreadStatistics_mod = function(localTreesDirectory="", nberOfExtractionFiles=1, timeSlices=200, onlyTipBranches=F, showingPlots=TRUE, outputName=gsub(" ","_",date()), nberOfCores=1, slidingWindow=NA, simulations=FALSE, discardExtractionTablesWithMoreThanOneAncestorForWavefrontPlot=FALSE) {

	nberOfStatistics = 6; treeIDs = c()
	registerDoMC(cores=nberOfCores)
	meanStatistics = matrix(nrow=(nberOfExtractionFiles), ncol=nberOfStatistics)
	branchVelocities = c() # not used, just to obtain an overall distribution of velocities
	sd_var_velocity = matrix(nrow=(nberOfExtractionFiles), ncol=2) # not used either
	medianMeanStatistics = matrix(nrow=1, ncol=nberOfStatistics)
	ciMeanStatistics = matrix(nrow=2, ncol=nberOfStatistics)
	waveFrontDistances1List = list()
	waveFrontDistances2List = list()
	meanDispersalVelocityList = list()
	weightedDispersalVelocityList = list()
	numberOfBranchesList = list()
	dispersalOrientationList = list()
	cat("Estimation of summary statistics", "\n", sep="")
	onlyOneAncestor = TRUE; extractionsWithMoreThanOneAncestors = c()
	if ((timeSlices==0)|is.na(timeSlices)|is.null(timeSlices)) onlyOneAncestor = FALSE
	dispersalVelocityGraph = FALSE
	if (!is.na(slidingWindow)) dispersalVelocityGraph = TRUE
	for (t in 1:nberOfExtractionFiles)
		{
			if (simulations == FALSE) fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
			if (simulations == TRUE) fileName = paste(localTreesDirectory,"/TreeSimulations_",t,".csv",sep="")
			data = read.csv(fileName, h=T)
			data = data[with(data, order(endYear,startYear)),]
			if (sum(!data[,"node1"]%in%data[,"node2"]) > 2)
				{
					onlyOneAncestor = FALSE; extractionsWithMoreThanOneAncestors = c(extractionsWithMoreThanOneAncestors, t)
				}
			treeIDs = c(treeIDs, data[1,"treeID"])
			if (onlyTipBranches == TRUE)
				{
					indices = c()
					for (i in 1:dim(data)[1])
						{
							n = data[i,"node2"]
							if (length(data[data[,"node1"]==n,"node1"]) == 0)
								{
									indices = c(indices, i)
								}
						}
					data = data[indices,]	
				}
			nberOfConnections = dim(data)[1]
			if (t == 1)
				{
					minLon = min(min(data[,"endLon"]),min(data[,"startLon"]))
					maxLon = max(max(data[,"endLon"]),max(data[,"startLon"]))
					minLat = min(min(data[,"endLat"]),min(data[,"startLat"]))
					maxLat = max(max(data[,"endLat"]),max(data[,"startLat"]))
					minStartYear = min(data[,"startYear"])
		    		minEndYear = min(data[,"startYear"])
		    		maxEndYear = max(data[,"endYear"])
				}	else	{
					if (minLon > min(min(data[,"endLon"]),min(data[,"startLon"]))) minLon = min(min(data[,"endLon"]),min(data[,"startLon"]))
					if (maxLon < max(max(data[,"endLon"]),max(data[,"startLon"]))) maxLon = max(max(data[,"endLon"]),max(data[,"startLon"]))
					if (minLat > min(min(data[,"endLat"]),min(data[,"startLat"]))) minLat = min(min(data[,"endLat"]),min(data[,"startLat"]))
					if (maxLat < max(max(data[,"endLat"]),max(data[,"startLat"]))) maxLat = max(max(data[,"endLat"]),max(data[,"startLat"]))
					if (minStartYear > min(data[,"startYear"])) minStartYear = min(data[,"startYear"])
					if (minEndYear > min(data[,"endYear"])) minEndYear = min(data[,"endYear"])
					if (maxEndYear < max(data[,"endYear"])) maxEndYear = max(data[,"endYear"])
				}
		}
	if ((onlyOneAncestor == FALSE)&(discardExtractionTablesWithMoreThanOneAncestorForWavefrontPlot == TRUE))
		{
			cat("Discarding",length(extractionsWithMoreThanOneAncestors),"extraction tables with more without a single","\n")
			cat("\t","common ancestor for generating wavefront plots","\n")
		}
	wavefrontDistanceSlices = timeSlices
	wavefrontDistanceTimeInterval = (maxEndYear-minStartYear)/wavefrontDistanceSlices
	dispersalVelocitySlices = timeSlices
	dispersalVelocityTimeInterval = (maxEndYear-minStartYear)/dispersalVelocitySlices
	xLim = c(minStartYear, maxEndYear)
	for (t in 1:nberOfExtractionFiles)
		{
			if (simulations == FALSE) fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
			if (simulations == TRUE) fileName = paste(localTreesDirectory,"/TreeSimulations_",t,".csv",sep="")
			data = read.csv(fileName, h=T)
			data = data[with(data, order(endYear,startYear)),]
			nberOfConnections = dim(data)[1]
			distances = matrix(0, nrow=dim(data)[1], ncol=1)
	 		colnames(distances) = c("greatCircleDist_km")
			for (i in 1:nberOfConnections)
				{
					x1 = cbind(data[i,"startLon"],data[i,"startLat"])
	 				x2 = cbind(data[i,"endLon"],data[i,"endLat"])
	 				distances[i] = rdist.earth(x1, x2, miles=F, R=NULL)
				}
			if (length(grep("greatCircleDist_km",colnames(data))) > 0)
				{
					data[,"greatCircleDist_km"] = distances
				}	else	{
					data = cbind(data,distances)
				}
			branchMeasures = matrix(nrow=nberOfConnections, ncol=2)
			weightedDispersalVelocity_numerator = 0; weightedDispersalVelocity_denominator = 0
			weightedDiffusionCoefficient_numerator = 0; weightedDiffusionCoefficient_denominator = 0
			for (i in 1:nberOfConnections)
				{
					dispersalTime = data[i,"endYear"]-data[i,"startYear"]
		    		branchVelocity = data[i,"greatCircleDist_km"]/(dispersalTime)
		    		branchOriginalDiffusionCoefficient = (data[i,"greatCircleDist_km"]^2)/(4*dispersalTime)
		    		branchMeasures[i,1] = branchVelocity
		    		branchMeasures[i,2] = branchOriginalDiffusionCoefficient
		    		weightedDispersalVelocity_numerator = weightedDispersalVelocity_numerator + data[i,"greatCircleDist_km"]
		    		weightedDispersalVelocity_denominator = weightedDispersalVelocity_denominator + dispersalTime
		    		weightedDiffusionCoefficient_numerator = weightedDiffusionCoefficient_numerator + (data[i,"greatCircleDist_km"]^2)
		    		weightedDiffusionCoefficient_denominator = weightedDiffusionCoefficient_denominator + (4*dispersalTime)
		  		}
		  	branchVelocities = c(branchVelocities, branchMeasures[,1])
			sd_var_velocity[t,] = cbind(sd(branchMeasures[,1]), var(branchMeasures[,1]))
			meanStatistics[t,1] = mean(branchMeasures[,1])
			meanStatistics[t,2] = weightedDispersalVelocity_numerator/weightedDispersalVelocity_denominator
			meanStatistics[t,3] = sd(branchMeasures[,1])/mean(branchMeasures[,1])
			meanStatistics[t,4] = mean(branchMeasures[,2])
		    meanStatistics[t,5] = weightedDiffusionCoefficient_numerator/weightedDiffusionCoefficient_denominator
			meanStatistics[t,6] = sd(branchMeasures[,2])/mean(branchMeasures[,2])
		}
	for (i in 1:nberOfStatistics)
		{
			medianMeanStatistics[1,i] = median(meanStatistics[,i], na.rm=T)
			quantiles = quantile(meanStatistics[,i], probs=c(0.025,0.975), na.rm=T)
			ciMeanStatistics[1,i] = as.numeric(quantiles[1]); ciMeanStatistics[2,i] = as.numeric(quantiles[2])
			HPD = HDInterval::hdi(meanStatistics[,i])[1:2]
			ciMeanStatistics[1,i] = as.numeric(HPD[1]); ciMeanStatistics[2,i] = as.numeric(HPD[2])
		}
	cat("Median value of mean branch dispersal velocity = ",medianMeanStatistics[1,1],"\n	95% HPD = [",ciMeanStatistics[1,1],", ",ciMeanStatistics[2,1],"]","\n",sep="")
	cat("Median value of weighted branch dispersal velocity = ",medianMeanStatistics[1,2],"\n	95% HPD = [",ciMeanStatistics[1,2],", ",ciMeanStatistics[2,2],"]","\n",sep="")	
	cat("Median value of original diffusion coefficient = ",medianMeanStatistics[1,4],"\n	95% HPD = [",ciMeanStatistics[1,4],", ",ciMeanStatistics[2,4],"]","\n",sep="")	
	cat("Median value of weighted diffusion coefficient = ",medianMeanStatistics[1,5],"\n	95% HPD = [",ciMeanStatistics[1,5],", ",ciMeanStatistics[2,5],"]","\n",sep="")	
	colnames(meanStatistics) = c("mean_branch_dispersal_velocity", "weighted_branch_dispersal_velocity", "branch_dispersal_velocity_variation_among_branches", 
	"original_diffusion_coefficient", "weighted_diffusion_coefficient", "diffusion_coefficient_variation_among_branches")
	write.table(meanStatistics, file=paste(outputName,"_estimated_dispersal_statistics.txt",sep=""), quote=F, row.names=F, sep="\t")

	LWD = 0.2
	if (nberOfExtractionFiles > 1)
		{
			if (showingPlots) dev.new(width=5, height=5)
			if (showingPlots == FALSE) pdf(paste(outputName, "_mean_branch_dispersal_velocity_variation.pdf",sep=""), width=5, height=5)
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0), lwd=LWD, col="gray30")
			# H = Hpi(cbind(log10(meanStatistics[,1]), meanStatistics[,3]))
			# kde = kde(cbind(log10(meanStatistics[,1]), meanStatistics[,3]), H=H)
			H = Hpi(cbind(meanStatistics[,1], meanStatistics[,3]))
			kde = kde(cbind(meanStatistics[,1], meanStatistics[,3]), H=H)
			text = "Kernel density estimates of mean branch dispersal velocity parameters"
			colours = c("#FFFFFF","#D2D3D3","#9D9FA3","#6A6A6D")
			xLab = "mean branch velocity"; yLab = "mean branch velocity variation among branches"
			plot(kde, display="filled.contour2", cont=c(50,75,95), col=colours, axes=F, ann=F)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.55, col.main="gray30")
			legend("bottomright", c("","95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=15, pt.cex=1.5, col=colours, box.lty=0, cex=0.6, y.intersp=1.2)
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_mean_branch_dispersal_velocity_variation.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()

			if (showingPlots) dev.new(width=5,height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_weighted_branch_dispersal_velocity_variation.pdf",sep=""),width=5,height=5)
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0), lwd=LWD, col="gray30")
			# H = Hpi(cbind(log10(meanStatistics[,2]), meanStatistics[,3]))
			# kde = kde(cbind(log10(meanStatistics[,2]), meanStatistics[,3]), H=H)
			H = Hpi(cbind(meanStatistics[,2], meanStatistics[,3]))
			kde = kde(cbind(meanStatistics[,2], meanStatistics[,3]), H=H)
			text = "Kernel density estimates of weighted branch dispersal velocity parameters"
			colours = c("#FFFFFF","#D2D3D3","#9D9FA3","#6A6A6D")
			xLab = "weighted dispersal velocity"; yLab="weighted dispersal velocity variation among branches"
			plot(kde, display="filled.contour2", cont=c(50,75,95), col=colours, axes=F, ann=F)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.55, col.main="gray30")
			legend("bottomright", c("","95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=15, pt.cex=1.5, col=colours, box.lty=0, cex=0.6, y.intersp=1.2)
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_weighted_branch_dispersal_velocity_variation.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()

			if (showingPlots) dev.new(width=5,height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_original_diffusion_coefficient_variation.pdf",sep=""),width=5,height=5)
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0), lwd=LWD, col="gray30")
			# H = Hpi(cbind(log10(meanStatistics[,4]), meanStatistics[,6]))
			# kde = kde(cbind(log10(meanStatistics[,4]), meanStatistics[,6]), H=H)
			H = Hpi(cbind(meanStatistics[,4], meanStatistics[,6]))
			kde = kde(cbind(meanStatistics[,4], meanStatistics[,6]), H=H)
			text1 = "Kernel density estimates of original diffusion coefficient parameters"; text2 = "(Pybus et al. 2012)"
			colours = c("#FFFFFF","#D2D3D3","#9D9FA3","#6A6A6D")
			xLab = "mean original diffusion coefficient"; yLab="diffusion coefficient variation among branches"
			plot(kde, display="filled.contour2", cont=c(50,75,95), col=colours, axes=F, ann=F)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text1, cex.main=0.6, col.main="gray30")
			legend("bottomright", c("","95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=15, pt.cex=1.5, col=colours, box.lty=0, cex=0.6, y.intersp=1.2)
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_original_diffusion_coefficient_variation.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
			
			if (showingPlots) dev.new(width=5,height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_weighted_diffusion_coefficient_variation.pdf",sep=""), width=5, height=5)
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0), lwd=LWD, col="gray30")
			# H = Hpi(cbind(log10(meanStatistics[,5]), meanStatistics[,6]))
			# kde = kde(cbind(log10(meanStatistics[,5]), meanStatistics[,6]), H=H)
			H = Hpi(cbind(meanStatistics[,5],meanStatistics[,6]))
			kde = kde(cbind(meanStatistics[,5],meanStatistics[,6]),H=H)
			text1 = "Kernel density estimates of weighted diffusion coefficient parameters"; text2 = "(Trovao et al. 2015)"
			colours = c("#FFFFFF","#D2D3D3","#9D9FA3","#6A6A6D")
			xLab = "mean weighted diffusion coefficient"; yLab="diffusion coefficient variation among branches"
			plot(kde, display="filled.contour2", cont=c(50,75,95), col=colours, axes=F, ann=F)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text1, cex.main=0.6, col.main="gray30")
			legend("bottomright", c("","95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=15, pt.cex=1.5, col=colours, box.lty=0, cex=0.6, y.intersp=1.2)
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_weighted_diffusion_coefficient_variation.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
		}
}
