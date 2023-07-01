rS_geo_pat_dist = function(localTreesDirectory, nberOfExtractionFiles)
	{
		rSs = rep(NA, nberOfExtractionFiles)
		for (j in 1:nberOfExtractionFiles)
			{
				tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), head=T)
				tipNodeIndices = which(!tab[,"node2"]%in%tab[,"node1"])
				distTree = matrix(nrow=length(tipNodeIndices), ncol=length(tipNodeIndices))
				for (k in 2:dim(distTree)[1])
					{
						for (l in 1:(k-1))
							{
								index1 = tipNodeIndices[k]
								index2 = tipNodeIndices[l]
								indices1 = index1; root = FALSE
								while (root == FALSE)
									{	
										if (tab[indices1[length(indices1)],"node1"]%in%tab[,"node2"])
											{
												indices1 = c(indices1, which(tab[,"node2"]==tab[indices1[length(indices1)],"node1"]))
											}	else	{
												root = TRUE
											}
									}
								indices2 = index2; root = FALSE
								while (root == FALSE)
									{	
										if (tab[indices2[length(indices2)],"node1"]%in%tab[,"node2"])
											{
												indices2 = c(indices2, which(tab[,"node2"]==tab[indices2[length(indices2)],"node1"]))
											}	else	{
												root = TRUE
											}
									}
								indices3 = indices1[which(indices1%in%indices2)]; patristic_dis = NULL
								if (length(indices3) == 0)
									{
										patristic_dis = sum(tab[c(indices1,indices2),"length"])
									}	else	{
										patristic_dis = sum(tab[c(indices1[which(!indices1%in%indices3)],indices2[which(!indices2%in%indices3)]),"length"])
									}
								distTree[k,l] = patristic_dis; distTree[l,k] = patristic_dis
							}
					}
				distsGeo = matrix(nrow=dim(distTree)[1], ncol=dim(distTree)[2])
				for (k in 2:dim(distsGeo)[1])
					{
						for (l in 1:(k-1))
							{
								index1 = tipNodeIndices[k]
								index2 = tipNodeIndices[l]
								x1 = cbind(tab[index1,"endLon"], tab[index1,"endLat"])
								x2 = cbind(tab[index2,"endLon"], tab[index2,"endLat"])
								distsGeo[k,l] = rdist.earth(x1, x2, miles=F, R=NULL)
								distsGeo[l,k] = distsGeo[k,l]
							}
					}						
				rSs[j] = cor(log(distTree[lower.tri(distTree)]),log(distsGeo[lower.tri(distsGeo)]), method="spearman")
				rSs[j] = cor(distTree[lower.tri(distTree)],distsGeo[lower.tri(distsGeo)], method="spearman")
			}
		return(rSs)
	}

