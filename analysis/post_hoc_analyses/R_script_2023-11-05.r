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

analysis = "Analyses_11-01-2023"; clades = c("SC1","SC2") # analyses performed for the 1° submission
analysis = "Analyses_18-10-2023"; clades = c("SC2") # analyses performed with the clock rate of the 1° submission but including two ghost sequences
analysis = "Analyses_19-10-2023"; clades = c("SC2_noGhost","SC2_ghosts") # analyses based on the new rate priors based on the late 2020 data
burnIn = 0; selected_NRRs = c("NRR14","NRR3")
mostRecentSamplingDates = c(2021.534, 2021.285)

humanPangolinSeqs = c("AY394995.1", # human sequence (SC1)
					  "MT040333.1","MT040334.1","MT040335.1","MT040336.1","MT072864.1", # "pangolin" clade (SC1)
					  "MN908947.3", # human sequence (SC2)
					  "MT121216.1", # "isolated" pangolin sequence (SC2)
					  "MT040336.1","MT040334.1","MT040335.1","MT072864.1","MT040333.1", # "pangolin" clade (SC2)
					  "Ghost1","Ghost2") # "ghosts" sequences included in some tests

renameFiles = FALSE
if (renameFiles) # has to be ran several times until no warnings are logged
	{
		directories1 = list.files(analysis)
		for (i in 1:length(directories1))
			{
				directories2 = list.files(paste0(analysis,"/",directories1[i]))
				for (j in 1:length(directories2))
					{
						directories3 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j]))
						for (k in 1:length(directories3))
							{
								files = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k]))
								for (l in 1:length(files))
									{
										file = paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k],"/",files[l])
										if (grepl("SC1.SC1",file)) file.rename(file, gsub("SC1.SC1","SC1",file))
										if (grepl("SC2.SC2",file)) file.rename(file, gsub("SC2.SC2","SC2",file))
										if (grepl("MCC_",file)) file.rename(file, gsub("MCC_","",file))
										if (grepl(".MCC.tre",file)) file.rename(file, gsub(".MCC.tre","_t.tree",file))
										if (grepl("PoWTransformed_",file)) file.rename(file, gsub("PoWTransformed_","",file))
										if (grepl("\\.rates\\.",file)) file.rename(file, gsub("\\.rates\\.","_rates_",file))
										if (grepl("_rates_",file)) file.rename(file, gsub("_rates_","_WLDV_",file))
										if (grepl("\\.diffusionCoefficient\\.",file)) file.rename(file, gsub("\\.diffusionCoefficient\\.","_WDC_",file))
										if (directories3[k] == "strictClock")
											{
												if ((grepl(".trees",file))&(!grepl("_a.trees",file))) file.rename(file, gsub(".trees","_a.trees",file))
											}
									}
							}
					}
			}
	}

# 1. Preparing the GIS files and the different colour scales

background = raster("All_Natural_Earth_files/Gray_background.tif")
borders = shapefile("All_Natural_Earth_files/International_borders.shp")
borders = subset(borders, (borders@data[,"adm0_a3_l"]=="CHN")&(borders@data[,"adm0_a3_r"]=="PRK")
						 |(borders@data[,"adm0_a3_l"]=="PRK")&(borders@data[,"adm0_a3_r"]=="CHN")
						 |(borders@data[,"adm0_a3_l"]=="KOR")&(borders@data[,"adm0_a3_r"]=="PRK")
						 |(borders@data[,"adm0_a3_l"]=="PRK")&(borders@data[,"adm0_a3_r"]=="KOR")
						 |(borders@data[,"adm0_a3_l"]=="CHN")&(borders@data[,"adm0_a3_r"]=="MMR")
						 |(borders@data[,"adm0_a3_l"]=="MMR")&(borders@data[,"adm0_a3_r"]=="CHN")					 
						 |(borders@data[,"adm0_a3_l"]=="CHN")&(borders@data[,"adm0_a3_r"]=="VNM")
						 |(borders@data[,"adm0_a3_l"]=="VNM")&(borders@data[,"adm0_a3_r"]=="CHN")
						 |(borders@data[,"adm0_a3_l"]=="CHN")&(borders@data[,"adm0_a3_r"]=="LAO")
						 |(borders@data[,"adm0_a3_l"]=="LAO")&(borders@data[,"adm0_a3_r"]=="CHN")
						 |(borders@data[,"adm0_a3_l"]=="THA")&(borders@data[,"adm0_a3_r"]=="MMR")
						 |(borders@data[,"adm0_a3_l"]=="MMR")&(borders@data[,"adm0_a3_r"]=="THA")					 
						 |(borders@data[,"adm0_a3_l"]=="LAO")&(borders@data[,"adm0_a3_r"]=="MMR")
						 |(borders@data[,"adm0_a3_l"]=="MMR")&(borders@data[,"adm0_a3_r"]=="LAO")					 
						 |(borders@data[,"adm0_a3_l"]=="THA")&(borders@data[,"adm0_a3_r"]=="LAO")
						 |(borders@data[,"adm0_a3_l"]=="LAO")&(borders@data[,"adm0_a3_r"]=="THA")					 
						 |(borders@data[,"adm0_a3_l"]=="THA")&(borders@data[,"adm0_a3_r"]=="KHM")
						 |(borders@data[,"adm0_a3_l"]=="KHM")&(borders@data[,"adm0_a3_r"]=="THA")					 
						 |(borders@data[,"adm0_a3_l"]=="VNM")&(borders@data[,"adm0_a3_r"]=="LAO")
						 |(borders@data[,"adm0_a3_l"]=="LAO")&(borders@data[,"adm0_a3_r"]=="VNM")					 
						 |(borders@data[,"adm0_a3_l"]=="KHM")&(borders@data[,"adm0_a3_r"]=="LAO")
						 |(borders@data[,"adm0_a3_l"]=="LAO")&(borders@data[,"adm0_a3_r"]=="KHM")					 
						 |(borders@data[,"adm0_a3_l"]=="VNM")&(borders@data[,"adm0_a3_r"]=="KHM")
						 |(borders@data[,"adm0_a3_l"]=="KHM")&(borders@data[,"adm0_a3_r"]=="VNM"))
provinces = shapefile("All_Chinese_shapefiles/Admin-2_polygons.shp")
provinces = spTransform(provinces, crs(background))
south_korea = getData("GADM", country="KOR", level=0)
north_korea = getData("GADM", country="PRK", level=0)
japan = getData("GADM", country="JPN", level=0)
laos = getData("GADM", country="LAO", level=0)
vietnam = getData("GADM", country="VNM", level=0)
cambodia = getData("GADM", country="KHM", level=0)
thailand = getData("GADM", country="THA", level=0)
myanmar = getData("GADM", country="MMR", level=0)
study_area = bind(provinces, south_korea, north_korea, japan, laos, vietnam, cambodia, thailand, myanmar)
if (writingFiles)
	{
		study_area_internal_nodes = gUnaryUnion(bind(provinces, south_korea, north_korea, laos, vietnam, cambodia, thailand, myanmar))
		areas = rep(NA, length(study_area_internal_nodes@polygons[[1]]@Polygons))
		for (i in 1:length(areas)) areas[i] = study_area_internal_nodes@polygons[[1]]@Polygons[[i]]@area
		study_area_internal_nodes@polygons[[1]]@Polygons = study_area_internal_nodes@polygons[[1]]@Polygons[which(areas==max(areas))]
		coords = study_area_internal_nodes@polygons[[1]]@Polygons[[1]]@coords
		sink(file=paste0("Internal_nodes_area.kml"))
		cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
		cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
		cat(paste("\t<polygon id=\"","internal_nodes_areas","\" samplingProbability=\"",1,"\">",sep=""))
		cat("\n")
		cat("\t\t<coordinates>"); cat("\n")
		for (j in 1:dim(coords)[1])
			{
				cat(paste("\t\t\t",coords[j,2],",",coords[j,1],",0",sep="")); cat("\n")
			}
		cat("\t\t</coordinates>"); cat("\n")
		cat("\t</polygon>"); cat("\n")
		cat("</kml>"); cat("\n")
		sink(NULL)
	}
background = mask(crop(background,study_area),study_area)
background[background[]==106] = NA; r = background; plottingElevation = FALSE
cols1 = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]
cols2 = gsub("FF","",rev(viridis::viridis(161))[1:101])
cols3 = paste0(hcl.colors(10, palette="terrain2")[3:10],"BF")

# 2. Extracting the spatio-temporal information embedded in annotated trees

directories1 = list.files(analysis)
for (i in 1:length(directories1)) # to retrieve and annotate the MCC trees for the strict clock and PoW-transformed trees
	{
		directories2 = list.files(paste0(analysis,"/",directories1[i]))
		for (j in 1:length(directories2))
			{
				directories3 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j]))
				directories3 = directories3[which((!grepl("statistics",directories3))&(!grepl(".txt",directories3)))]
				for (k in 1:length(directories3))
					{
						files = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k]))
						for (l in 1:length(files))
							{
								extension = unlist(strsplit(files[l],"\\."))[length(unlist(strsplit(files[l],"\\.")))]
								file = paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k],"/",files[l])
								if ((!is.null(extension))&&(!is.na(extension))&&(extension == "trees"))
									{
										txt = scan(file, what="", sep="\n", quiet=T, blank.lines.skip=F)
										if ((txt[length(txt)]!="End;")&(txt[length(txt)]!="END;"))
											{
												write(c(txt,"End;"), file)
											}
										if ((directories3[k] == "strictClock")&(grepl("_a.trees",file)))
											{
												txt1 = scan(file, what="", sep="\n", quiet=T, blank.lines.skip=F)
												temp = gsub("strictClock","transPoW",gsub("_a.trees",".trees",file))
												allTrees2 = readAnnotatedNexus(temp); txt2 = c()
												states = gsub("\tTREE \\* ","",names(allTrees2))
												for (m in 1:length(txt1)) # step to retrieve the 100 trees that were subsequentely PoW-transformed
													{
														if (!grepl("STATE_",txt1[m]))
															{
																txt2 = c(txt2, txt1[m])
															}	else	{
																state = unlist(strsplit(txt1[m]," "))[2]
																if (state%in%states) txt2 = c(txt2, txt1[m])
															}
													}
												write(txt2, gsub("_a.trees","_b.trees",file))
												temp = gsub("strictClock","transPoW",gsub("_a.trees",".tree",file))
												mccTreeTransPoW = read.nexus(temp); allTrees1 = read.nexus(gsub("_a.trees","_b.trees",file))
												dimensions = dim(mccTreeTransPoW$edge); index = c()
												for (m in 1:length(allTrees1))
													{
														if (sum(allTrees1[[m]]$edge == mccTreeTransPoW$edge) == (dimensions[1]*dimensions[2]))
															{
																index = c(index,m)
															}
													}
												if (length(index) != 1)
													{
														print(c(i,j,k,l,m))
													}	else	{
														write.tree(allTrees1[[index]], file=gsub("_a.trees","_t.tree",file))
													}
												file1 = gsub("_a.trees","_b.trees",file); file2 = gsub("_a.trees","_b.tree",file); file3 = gsub("_a.trees","_t.tree",file)					
												system(paste0("BEAST_1104_program/bin/treeannotator -burninTrees 0 -heights keep -target ",file3," ",file1," ",file2), 
													   ignore.stdout=F, ignore.stderr=F)
											}
										if ((directories3[k] == "transPoW")&(grepl(".trees",file)))
											{
												file1 = file; file2 = gsub(".trees",".tree",file)
												system(paste0("BEAST_1104_program/bin/treeannotator -burninTrees 0 -heights keep ",file1," ",file2),
													   ignore.stdout=F, ignore.stderr=F)
											}
									}
							}
					}
			}
	}
directories1 = list.files(analysis)
for (i in 1:length(directories1)) # to extract the spatio-temporal information embedded in all MCC and posterior trees
	{
		metadata = read.csv(paste0(gsub("_ghosts","",gsub("_noGhost","",directories1[i])),"_all_metadata.csv"), head=T, sep=";")
		samplingDates = gsub("\\/","-",metadata[,"date"]); maxSamplingDate = 0
		for (j in 1:length(samplingDates))
			{
				samplingDate = gsub("\\/","-",samplingDates[j])
				samplingDate1 = unlist(strsplit(samplingDate,"-"))
				if (length(samplingDate1) > 1)
					{
						if (length(samplingDate1) == 2)
							{
								if (nchar(samplingDate1[1]) == 4) samplingDate = paste0(samplingDate,"-15")
								if (nchar(samplingDate1[2]) == 4) samplingDate = paste0("15-",samplingDate)
								samplingDate1 = unlist(strsplit(samplingDate,"-"))
							}
						samplingDate2 = NA
						if (length(samplingDate1) == 3)
							{
								if (nchar(samplingDate1[1]) == 4) samplingDate2 = decimal_date(ymd(samplingDate))
								if (nchar(samplingDate1[3]) == 4) samplingDate2 = decimal_date(dmy(samplingDate))	
							}
						if (maxSamplingDate < samplingDate2) maxSamplingDate = samplingDate2
					}
			}
		mostRecentSamplingDatum = maxSamplingDate; # print(mostRecentSamplingDatum)
		directories2 = list.files(paste0(analysis,"/",directories1[i]))
		for (j in 1:length(directories2))
			{
				directories3 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j]))
				directories3 = directories3[which((!grepl("statistics",directories3))&(!grepl(".txt",directories3)))]
				# directories3 = directories3[which(grepl("timeTrees",directories3))]
				for (k in 1:length(directories3))
					{
						files = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k]))
						files = files[which((!grepl("_a.trees",files))&(!grepl("_t.tree",files)))]
						for (l in 1:length(files))
							{
								extension = unlist(strsplit(files[l],"\\."))[length(unlist(strsplit(files[l],"\\.")))]
								file = paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k],"/",files[l])
								if (extension == "tree")
									{
										tree1a = readAnnotatedNexus(file); tree1b = read.nexus(file)
										n_tips = length(tree1a$tip.label); n_splits = (2*n_tips)-3-n_tips
										if (directories3[k] == "strictClock")
											{
												temp = gsub("strictClock","transPoW",gsub("_b.tree",".tree",file))
												tree2a = readAnnotatedNexus(temp); tree2b = read.nexus(temp)
												comparison = comparePhylo(tree1a, tree2a, plot=F, force.rooted=F)
												if (as.numeric(unlist(strsplit(comparison$messages[8]," "))[1]) != n_splits)
													{
														print(c(i,j,k,l,m))
													}
												dimensions = dim(tree2b$edge)
												if (sum(tree1b$edge == tree2b$edge) != (dimensions[1]*dimensions[2]))
													{
														# print(c(i,j,k,l,m))
													}
												tree1a$edge = tree2a$edge # NEW !!
											}
										tab = Tree_data_extraction1(tree1a, mostRecentSamplingDatum)
										hosts = matrix(nrow=dim(tab)[1], ncol=1); colnames(hosts) = "host"
										missingSequences = c()
										for (m in 1:dim(tab)[1])
											{
												if (!is.na(tab[m,"tipLabel"]))
													{
														index = which(metadata[,"name"]==gsub("SC2_","",gsub("SC1_","",tab[m,"tipLabel"])))
														if (length(index) == 0)
															{
																missingSequences = c(missingSequences, tab[m,"tipLabel"])
															}	else	{
																hosts[m,"host"] = metadata[index,"host"]
															}
													}
											}
										if (!"host"%in%colnames(tab))
											{
												tab = cbind(tab, hosts)
											}	else	{
												tab[,"host"] = hosts
											}
										write.csv(tab, gsub("\\.tree","\\_1.csv",file), row.names=F, quote=F)
									}
								if (extension == "trees")
									{
										localTreesDirectory = gsub("\\.trees","_ext1",file)
										dir.create(file.path(localTreesDirectory), showWarnings=F)
										txt = scan(file, what="", sep="\n", quiet=T, blank.lines.skip=F)
										if ((txt[length(txt)]!="End;")&(txt[length(txt)]!="END;"))
											{
												write(c(txt,"End;"), file)
											}
										allTrees1a = readAnnotatedNexus(file)
										allTrees1b = read.nexus(file)
										if (directories3[k] == "strictClock")
											{
												temp = gsub("strictClock","transPoW",gsub("_b.trees",".trees",file))
												allTrees2a = readAnnotatedNexus(temp)
												allTrees2b = read.nexus(temp); buffer = allTrees2a
												names(allTrees2a) = gsub("\tTREE \\* ","",names(allTrees2a))
												for (m in 1:length(buffer))
													{
														index = which(names(allTrees1a)==names(allTrees2a)[m])
														dimensions = dim(allTrees1b[[index]]$edge)
														if (sum(allTrees1b[[index]]$edge == allTrees2b[[m]]$edge) != (dimensions[1]*dimensions[2]))
															{
																print(c(i,j,k,l,m))
															}
														buffer[[m]] = allTrees1a[[index]]
														buffer[[m]]$edge = allTrees2a[[m]]$edge # NEW !!
														names(buffer)[m] = names(allTrees2a)[m]
													}
												allTrees1a = buffer
											}
										if (directories3[k] == "timeTrees")
											{
												temp = gsub("timeTrees","transPoW",gsub("_b.trees",".trees",file))
												if (file.exists(temp))
													{
														allTrees2a = readAnnotatedNexus(temp); allTrees2b = read.nexus(temp)
														buffer = allTrees2a # not used anymore because messing with the tree IDs !!
														names(allTrees2a) = gsub("\tTREE \\* ","",names(allTrees2a)); indices = length(allTrees2a)
														for (m in 1:length(allTrees2a)) # trees (and tree IDs) actually do not correspond between both posterior distributions
															{			 	  # (but this is then simply a way to randomly sample 100 trees from the "timeTrees" distribution)
																indices[m] = which(names(allTrees1a)==names(allTrees2a)[m])
															}
														allTrees1a = allTrees1a[indices]
													}	else	{
														allTrees1a = allTrees1a[(length(allTrees1a)-99):(length(allTrees1a))] # to take the last 100 trees (assuring the remove the burn-in)
													}
											}
										for (m in 1:length(allTrees1a))
											{
												tab = Tree_data_extraction2(allTrees1a[[m]], mostRecentSamplingDatum)
												write.csv(tab, paste0(localTreesDirectory,"/TreeExtractions_",m,".csv"), row.names=F, quote=F)
											}
									}
							}
					}
			}
	}
directories1 = list.files(analysis)
for (i in 1:length(directories1)) # to import the continuous phylogeographic annotations within the PoW-transformed trees
	{
		directories2 = list.files(paste0(analysis,"/",directories1[i]))
		for (j in 1:length(directories2))
			{
				directories3 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j]))
				directories3 = directories3[which(grepl("transPoW",directories3))]
				for (k in 1:length(directories3))
					{
						files1 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k]))
						for (l in 1:length(files1))
							{
								extension = unlist(strsplit(files1[l],"\\."))[length(unlist(strsplit(files1[l],"\\.")))]
								file = paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k],"/",files1[l])
								if (extension == "csv")
									{
										tab1 = read.csv(file, head=T)
										temp = gsub("transPoW","strictClock",gsub(".csv","_b.csv",file))
										tab2 = read.csv(temp, head=T)
										for (n in 1:dim(tab1)[1])
											{
												index = which((tab2[,"node1"]==tab1[n,"node1"])&(tab2[,"node2"]==tab1[n,"node2"]))
												tab1[n,"startLon"] = tab2[index,"startLon"]; tab1[n,"startLat"] = tab2[index,"startLat"]
												tab1[n,"endLon"] = tab2[index,"endLon"]; tab1[n,"endLat"] = tab2[index,"endLat"]
											}
										write.csv(tab1, file, row.names=F, quote=F)
									}
								if (grepl("_ext1",file))
									{
										files2 = list.files(file)
										for (m in 1:length(files2))
											{
												tab1 = read.csv(paste0(file,"/",files2[m]), head=T)
												temp = gsub("transPoW","strictClock",gsub("_ext1","_b_ext1",file))
												tab2 = read.csv(paste0(temp,"/",files2[m]), head=T)
												for (n in 1:dim(tab1)[1])
													{
														index = which((tab2[,"node1"]==tab1[n,"node1"])&(tab2[,"node2"]==tab1[n,"node2"]))
														tab1[n,"startLon"] = tab2[index,"startLon"]; tab1[n,"startLat"] = tab2[index,"startLat"]
														tab1[n,"endLon"] = tab2[index,"endLon"]; tab1[n,"endLat"] = tab2[index,"endLat"]
													}
												write.csv(tab1, paste0(file,"/",files2[m]), row.names=F, quote=F)
											}
									}
							}
					}
			}
	}
directories1 = list.files(analysis)
for (i in 1:length(directories1)) # to generate extraction files that do not include the clades corresponding to human and pangolin sequences (for the strict clock and PoW-transformed trees)
	{
		directories2 = list.files(paste0(analysis,"/",directories1[i]))
		for (j in 1:length(directories2))
			{
				directories3 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j]))
				directories3 = directories3[which((!grepl(".txt",directories3))&(!grepl("timeTrees",directories3))&(!grepl("statistics",directories3)))]
				for (k in 1:length(directories3))
					{
						files1 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k]))
						for (l in 1:length(files1))
							{
								extension = unlist(strsplit(files1[l],"\\."))[length(unlist(strsplit(files1[l],"\\.")))]
								file = paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k],"/",files1[l])
								if (grepl("_ext1",file))
									{
										localTreesDirectory1 = file; localTreesDirectory2 = gsub("_ext1","_ext2",localTreesDirectory1) # negative branch lengths set to 0
													# and human/pangolin tip branches (and clades) discarded in the context of the "noHumPanLoc" phylogeographic analyses
										dir.create(file.path(localTreesDirectory2), showWarnings=F); files2 = list.files(localTreesDirectory1)
										for (m in 0:length(files2))
											{
												if (m == 0) tab1 = read.csv(gsub("_ext1","_1.csv",localTreesDirectory1), head=T)
												if (m != 0) tab1 = read.csv(paste0(localTreesDirectory1,"/",files2[m]), head=T)
												if (directories2[j]=="noHumPanLoc") # to discard human and pangolin tips, as well as the "pangolin" clade for SC2
													{
														nodeClades = matrix(nrow=dim(tab1), ncol=2); colnames(nodeClades) = c("startClade","endClade")
														nodeClades[which(!tab1[,"node2"]%in%tab1[,"node1"]),"endClade"] = "main"
														nodeClades[which(gsub("SC1_","",gsub("SC2_","",tab1[,"tipLabel"]))%in%humanPangolinSeqs),"endClade"] = "humanPangolin"
														tab1 = cbind(tab1, nodeClades); noRemainingCladeToIdentify = FALSE
														while (noRemainingCladeToIdentify == FALSE)
															{
																for (n in 1:dim(tab1)[1])
																	{
																		if (is.na(tab1[n,"endClade"]))
																			{
																				indices = which(tab1[,"node1"]==tab1[n,"node2"])
																				if (sum(!is.na(tab1[indices,"endClade"])) == 2)
																					{
																						if (tab1[indices[1],"endClade"] == tab1[indices[2],"endClade"])
																							{
																								tab1[indices[1],"startClade"] = tab1[indices[1],"endClade"]
																								tab1[indices[2],"startClade"] = tab1[indices[2],"endClade"]
																								tab1[n,"endClade"] = tab1[indices[1],"endClade"]
																							}	else	{
																								tab1[indices[1],"startClade"] = "main"
																								tab1[indices[2],"startClade"] = "main"
																								tab1[n,"endClade"] = "main"
																							}
																					}
																			}
																	}
																if (sum(is.na(tab1[,"startClade"])) == 2)
																	{
																		tab1[which(!tab1[,"node1"]%in%tab1[,"node2"]),"startClade"] = "main"
																	}				
																if (sum(is.na(tab1[,"startClade"])) == 0) noRemainingCladeToIdentify = TRUE
															}
														tab1 = tab1[-which(tab1[,"endClade"]=="humanPangolin"),]
													}
												tab2 = tab1; correctedBranches = c()
												for (n in 1:dim(tab2)[1])
													{
														if (tab2[n,"length"] < 0)
															{
																buffer = tab2[n,"length"]; tab2[n,"length"] = 0
																tab2[n,"startYear"] = tab2[n,"startYear"]+buffer
																branchToCorrect = which(tab2[,"node2"]==tab2[n,"node1"])
																if (!branchToCorrect%in%correctedBranches)
																	{
																		tab2[branchToCorrect,"length"] = tab2[branchToCorrect,"length"]+buffer
																		tab2[branchToCorrect,"endYear"] = tab2[branchToCorrect,"endYear"]+buffer
																		correctedBranches = c(correctedBranches, branchToCorrect)
																	}
															}				
													}
												if (m == 0) write.csv(tab2, gsub("_ext2","_2.csv",localTreesDirectory2), row.names=F, quote=F)
												if (m != 0) write.csv(tab2, paste0(localTreesDirectory2,"/",files2[m]), row.names=F, quote=F)
											}
										localTreesDirectory3 = gsub("_ext1","_ext3",localTreesDirectory1)
										dir.create(file.path(localTreesDirectory3), showWarnings=F); files2 = list.files(localTreesDirectory2)
										for (m in 0:length(files2))
											{
												if (m == 0) tab2 = read.csv(gsub("_ext2","_2.csv",localTreesDirectory2), head=T)
												if (m != 0) tab2 = read.csv(paste0(localTreesDirectory2,"/",files2[m]), head=T)
												tab3 = tab2; tab3 = tab3[which(tab2[,"length"]>0),]
												if (m == 0) write.csv(tab3, gsub("_ext3","_3.csv",localTreesDirectory3), row.names=F, quote=F)
												if (m != 0) write.csv(tab3, paste0(localTreesDirectory3,"/",files2[m]), row.names=F, quote=F)
											}
									}
							}
					}
			}
	}
directories1 = list.files(analysis)
for (i in 1:length(directories1)) # to generate extraction files that do not include the clades corresponding to human, pangolin, and ghost sequences (for the "timeTrees" trees)
	{
		directories2 = list.files(paste0(analysis,"/",directories1[i]))
		for (j in 1:length(directories2))
			{
				directories3 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j]))
				directories3 = directories3[which(grepl("timeTrees",directories3))]
				for (k in 1:length(directories3))
					{
						files1 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k]))
						for (l in 1:length(files1))
							{
								extension = unlist(strsplit(files1[l],"\\."))[length(unlist(strsplit(files1[l],"\\.")))]
								file = paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k],"/",files1[l])
								if (grepl("_ext1",file))
									{
										localTreesDirectory1 = file; localTreesDirectory2 = gsub("_ext1","_ext2",localTreesDirectory1) # negative branch lengths set to 0
													# and human/pangolin tip branches (and clades) discarded in the context of the "noHumPanLoc" phylogeographic analyses
										dir.create(file.path(localTreesDirectory2), showWarnings=F); files2 = list.files(localTreesDirectory1)
										for (m in 1:length(files2))
											{
												tab1 = read.csv(paste0(localTreesDirectory1,"/",files2[m]), head=T)
												nodeClades = matrix(nrow=dim(tab1), ncol=2); colnames(nodeClades) = c("startClade","endClade")
												nodeClades[which(!tab1[,"node2"]%in%tab1[,"node1"]),"endClade"] = "main"
												nodeClades[which(gsub("SC1_","",gsub("SC2_","",tab1[,"tipLabel"]))%in%humanPangolinSeqs),"endClade"] = "humanPangolin"
												tab1 = cbind(tab1, nodeClades); noRemainingCladeToIdentify = FALSE
												while (noRemainingCladeToIdentify == FALSE)
													{
														for (n in 1:dim(tab1)[1])
															{
																if (is.na(tab1[n,"endClade"]))
																	{
																		indices = which(tab1[,"node1"]==tab1[n,"node2"])
																		if (sum(!is.na(tab1[indices,"endClade"])) == 2)
																			{
																				if (tab1[indices[1],"endClade"] == tab1[indices[2],"endClade"])
																					{
																						tab1[indices[1],"startClade"] = tab1[indices[1],"endClade"]
																						tab1[indices[2],"startClade"] = tab1[indices[2],"endClade"]
																						tab1[n,"endClade"] = tab1[indices[1],"endClade"]
																					}	else	{
																						tab1[indices[1],"startClade"] = "main"
																						tab1[indices[2],"startClade"] = "main"
																						tab1[n,"endClade"] = "main"
																					}
																			}
																	}
															}
														if (sum(is.na(tab1[,"startClade"])) == 2)
															{
																tab1[which(!tab1[,"node1"]%in%tab1[,"node2"]),"startClade"] = "main"
															}				
														if (sum(is.na(tab1[,"startClade"])) == 0) noRemainingCladeToIdentify = TRUE
													}
												tab2 = tab1[-which(tab1[,"endClade"]=="humanPangolin"),]
												write.csv(tab2, paste0(localTreesDirectory2,"/",files2[m]), row.names=F, quote=F)
											}
									}
							}
					}
			}
	}

# 3. Investigating the profile of long-distance dispersal events

directories1 = list.files(analysis); lines = c()
for (i in 1:length(directories1))
	{
		directories2 = list.files(paste0(analysis,"/",directories1[i]))
		for (j in 1:length(directories2))
			{
				directories3 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j]))
				directories3 = directories3[which(grepl("transPoW",directories3))]
				for (k in 1:length(directories3))
					{
						files = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k]))
						files = files[which(grepl(".trees",files))]
						for (l in 1:length(files))
							{
								file = paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k],"/",files[l])
								trees = read.nexus(file)
								for (m in 1:length(trees))
									{
										n = sum(trees[[m]]$edge.length<0)
										if (n > 0)
											{
												lines = c(lines, paste(directories1[i],directories2[j],directories3[k],files[l],names(trees)[m],paste0(n," branch length(s) < 0"),sep=" - "))
											}
									}
							}
					}
			}
	}
if (writingFiles)
	{
		write(lines, "Branch_lengths_<_0.txt")
	}
blue1 = "#217B99"; blue2 = "#217B9950"; orange1 = "#FAA521"; orange2 = "#FAA52150"
directories = c("/SC1/noHumPanLoc/transPoW/","/SC2/noHumPanLoc/transPoW/")
for (h in 1:length(directories))
	{
		files = list.files(paste0(analysis,directories[h]))
		files = files[which(grepl("_3.csv",files))]
		for (i in 1:length(files))
			{
				mcc = read.csv(paste0(analysis,directories[h],files[i]), head=T)
				dists = rep(NA, dim(mcc)[1])
				for (j in 1:dim(mcc)[1])
					{
						x1 = cbind(mcc[j,"startLon"], mcc[j,"startLat"])
						x2 = cbind(mcc[j,"endLon"], mcc[j,"endLat"])
						dists[j] = rdist.earth(x1, x2, miles=F, R=NULL)
					}
				rS = cor(log(dists),log(mcc[,"length"]), method="spearman")		
				pdf(paste0(analysis,directories[h],gsub(".csv",".pdf",files[i])), width=5, height=4.4); par(mgp=c(0,0,0), oma=c(0,0,0,0), mar=c(3.3,3.5,2,2), lwd=0.2, col="gray30")
				# pdf(paste0(gsub(".csv",".pdf",files[i])), width=5, height=4.4); par(mgp=c(0,0,0), oma=c(0,0,0,0), mar=c(3.3,3.5,2,2), lwd=0.2, col="gray30")
					# plot(log(dists), log(mcc[,"length"]), col=blue2, pch=16, cex=0.8, axes=F, ann=F, frame=T); points(log(dists), log(mcc[,"length"]), col=blue1, pch=1, cex=0.8)
				# axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
				# axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.30,0), lwd=0.0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
				# title(ylab="phylogenetic branch durations (years, log-transformed)", cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
				# title(xlab="geographic distance (km, log-transformed)", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
				# title(main=paste0("r (Spearman) = ",round(rS,2)), cex.main=0.7, col.main="gray30", line=-1.5)			
					# plot(dists, dists/mcc[,"length"], col=blue2, pch=16, cex=0.8, axes=F, ann=F, frame=T); points(dists, dists/mcc[,"length"], col=blue1, pch=1, cex=0.8)
				# axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
				# axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.30,0), lwd=0.0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
				# title(ylab="lineage dispersal velocity (km/years)", cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
				# title(xlab="geographic distance (km)", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
					plot(log(dists), log(dists/mcc[,"length"]), col=blue2, pch=16, cex=0.8, axes=F, ann=F, frame=T)
				points(log(dists), log(dists/mcc[,"length"]), col=blue1, pch=1, cex=0.8)
				axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
				axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.30,0), lwd=0.0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
				title(ylab="lineage dispersal velocity (km/years, log-transformed)", cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
				title(xlab="geographic distance (km, log-transformed)", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
				dev.off()
			}
	}

# 4. Investigating the patterns of isolation-by-distance

directories = c("/SC1/noHumPanLoc/transPoW/","/SC2/noHumPanLoc/transPoW/"); rS_list_1 = list()
for (h in 1:length(directories)) # to compute the correlation between the dispersal duration and geographic distance associated with each branch (not used anymore, and was based on PoW-transformed trees)
	{
		folders = list.files(paste0(analysis,directories[h]))
		folders = folders[which(grepl("_ext3",folders))]
		rS_list_2 = list(); all_rSs = c()
		for (i in 1:length(folders))
			{
				files = list.files(paste0(analysis,directories[h],folders[i]))
				rSs = rep(NA, length(files))
				for (j in 1:length(files))
					{
						tab = read.csv(paste0(analysis,directories[h],folders[i],"/",files[j]), head=T)
						dists = rep(NA, dim(tab)[1])
						for (k in 1:dim(tab)[1])
							{
								x1 = cbind(tab[k,"startLon"], tab[k,"startLat"])
								x2 = cbind(tab[k,"endLon"], tab[k,"endLat"])
								dists[k] = rdist.earth(x1, x2, miles=F, R=NULL)
							}
						rSs[j] = cor(log(dists),log(tab[,"length"]), method="spearman")
					}
				rS_list_2[[i]] = rSs
			}
		for (i in 1:length(rS_list_2)) all_rSs = c(all_rSs, rS_list_2[[i]])
		HPD = round(HDInterval::hdi(all_rSs)[1:2],2)
		cat("\t",clades[h],": median rS = ",round(median(all_rSs),2),", 95% HPD = [",HPD[1],"-",HPD[2],"]",sep="")
		rS_list_1[[h]] = rS_list_2
	}

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
				tempT = read.nexus(paste0(analysis,gsub("timeTrees","transPoW",directories[h]),files[i]))
				tips = trees[[1]]$tip.label; rSs = rep(NA, length(files))
				for (j in 1:length(tempT))
					{
						vSs = list(); index = which(names(trees)==names(tempT)[j]); tree = trees[[index]]
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
		# SC1, cut-off = 50 years: 	median rS = 0.53, 95% HPD = [0.27-0.76]
		# SC2, cut-off = 50 years: 	median rS = 0.53, 95% HPD = [0.02-0.82]
		# SC1, cut-off = 75 years: 	median rS = 0.56, 95% HPD = [0.39-0.76]
		# SC2, cut-off = 75 years: 	median rS = 0.62, 95% HPD = [0.16-0.88]
		# SC1, no cut-off: 			median rS = 0.45, 95% HPD = [0.36-0.59]
		# SC2, no cut-off: 			median rS = 0.66, 95% HPD = [0.12-0.83]

		# SC1, cut-off = 50 years: 	median vS = 11.5, 95% HPD = [0-37.6]
		# SC2, cut-off = 50 years: 	median vS = 9.4, 95% HPD = [0-28.8]
		# SC1, cut-off = 75 years: 	median vS = 11.0, 95% HPD = [0-31.3]
		# SC2, cut-off = 75 years: 	median vS = 9.8, 95% HPD = [0-26.3]
		# SC1, no cut-off: 			median vS = 8.2, 95% HPD = [0-26.4]
		# SC2, no cut-off: 			median vS = 7.7, 95% HPD = [0-21.9]

# 5. Estimating some dispersal statistics for each reconstruction

directories1 = list.files(analysis); cutoffs = c(75,50)
for (i in 1:length(directories1))
	{
		directories2 = list.files(paste0(analysis,"/",directories1[i]))
		for (j in 1:length(directories2))
			{
				directories3 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j]))
				directories3 = directories3[which(grepl("timeTrees",directories3))]
				for (k in 1:length(directories3))
					{
						files = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k]))
						if (directories2[j] == "humPangLocSI") files = files[which(grepl("_ext1",files))]
						if (directories2[j] == "noHumPanLoc") files = files[which(grepl("_ext2",files))]
						for (l in 1:length(files))
							{
								file = paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k],"/",files[l])
								localTreesDirectory1 = file; nberOfExtractionFiles = length(list.files(localTreesDirectory1))
								for (m in 1:length(cutoffs))
									{
										if (directories2[j] == "humPangLocSI") localTreesDirectory2 = gsub("_ext1",paste0("_c",cutoffs[m]),file)
										if (directories2[j] == "noHumPanLoc") localTreesDirectory2 = gsub("_ext2",paste0("_c",cutoffs[m]),file)
										dir.create(file.path(localTreesDirectory2), showWarnings=F)
										for (n in 1:nberOfExtractionFiles)
											{
												tab1 = read.csv(paste0(localTreesDirectory1,"/TreeExtractions_",n,".csv"), head=T)
												tab2 = tab1[which(tab1[,"startYear"]>(max(tab1[,"endYear"])-cutoffs[m])),]
												write.csv(tab2, paste0(localTreesDirectory2,"/TreeExtractions_",n,".csv"), row.names=F, quote=F)
											}
									}
							}
					}
			}
	}
source("spreadStatistics_mod.r"); cutoffs = c(75,50)
for (i in 1:length(directories1))
	{
		directories2 = list.files(paste0(analysis,"/",directories1[i]))
		for (j in 1:length(directories2))
			{
				directories3 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j]))
				directories3 = directories3[which(grepl("timeTrees",directories3))]
				for (k in 1:length(directories3))
					{
						files = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k]))
						if (directories2[j] == "humPangLocSI") files = files[which(grepl("_ext1",files))]
						if (directories2[j] == "noHumPanLoc") files = files[which(grepl("_ext2",files))]
						for (l in 1:length(files))
							{
								file = paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k],"/",files[l])
								localTreesDirectory1 = file; nberOfExtractionFiles = length(list.files(localTreesDirectory1))
								for (m in 1:length(cutoffs))
									{
										if (directories2[j] == "humPangLocSI") localTreesDirectory2 = gsub("_ext1",paste0("_c",cutoffs[m]),file)
										if (directories2[j] == "humPangLocSI") statisticsDirectory = gsub("_ext1",paste0("_s",cutoffs[m]),file)
										if (directories2[j] == "noHumPanLoc") localTreesDirectory2 = gsub("_ext2",paste0("_c",cutoffs[m]),file)
										if (directories2[j] == "noHumPanLoc") statisticsDirectory = gsub("_ext2",paste0("_s",cutoffs[m]),file)
										dir.create(file.path(statisticsDirectory), showWarnings=F)
										timeSlices = 100; onlyTipBranches = FALSE; showingPlots = FALSE; nberOfCores = 10; slidingWindow = 1
										if (directories2[j] == "humPangLocSI") outputName = gsub("_ext1","",unlist(strsplit(file,"\\/"))[length(unlist(strsplit(file,"\\/")))])
										if (directories2[j] == "noHumPanLoc") outputName = gsub("_ext2","",unlist(strsplit(file,"\\/"))[length(unlist(strsplit(file,"\\/")))])
										outputName = paste0(statisticsDirectory,"/",outputName)
										spreadStatistics_mod(localTreesDirectory2, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
									}
							}
					}
			}
	}
cutoffs = c(75,50); minX = 9999; maxX = -9999; minY = 9999; maxY = -9999
for (i in 1:length(directories1)) # to estimate the maximum, mean, and median branch velocities (not considering human/pangolin branches/clades)
	{
		for (m in 1:length(cutoffs))
			{
				allVs = c(); allBranchVelocities = list(); maxBranchVelocities = c()
				meanBranchVelocities = c(); medianBranchVelocities = c()
				directories2 = list.files(paste0(analysis,"/",directories1[i]))
				directories2 = directories2[which(grepl("noHumPanLoc",directories2))]
				for (j in 1:length(directories2))
					{
						directories3 = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j]))
						directories3 = directories3[which(grepl("timeTrees",directories3))]
						for (k in 1:length(directories3))
							{
								files = list.files(paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k]))
								files = files[which(grepl("_ext2",files))]
								for (l in 1:length(files))
									{
										file = paste0(analysis,"/",directories1[i],"/",directories2[j],"/",directories3[k],"/",files[l])
										localTreesDirectory = gsub("_ext2",paste0("_c",cutoffs[m]),file)
										nberOfExtractionFiles = length(list.files(localTreesDirectory))
										line_allVs = c()
										line_maxVs = matrix(nrow=1, ncol=nberOfExtractionFiles)
										line_meanVs = matrix(nrow=1, ncol=nberOfExtractionFiles)
										line_medianVs = matrix(nrow=1, ncol=nberOfExtractionFiles)
										for (n in 1:nberOfExtractionFiles)
											{
												tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",n,".csv"), head=T)
												x1 = cbind(tab[,"startLon"], tab[,"startLat"])
												x2 = cbind(tab[,"endLon"], tab[,"endLat"])
												if (minX > x1[1,1]) minX = x1[1,1]
												if (maxX < x1[1,1]) maxX = x1[1,1]
												if (minY > x1[1,2]) minY = x1[1,2]
												if (maxY < x1[1,2]) maxY = x1[1,2]
												if (minX > x2[1,1]) minX = x2[1,1]
												if (maxX < x2[1,1]) maxX = x2[1,1]
												if (minY > x2[1,2]) minY = x2[1,2]
												if (maxY < x2[1,2]) maxY = x2[1,2]
												geoDists = rdist.earth(x1, x2, miles=F, R=NULL)
												if (length(diag(geoDists)) != length(tab[,"length"])) print(c(i,m,j,k,l))
												branchVelocities = diag(geoDists)/tab[,"length"]
												allVs = c(allVs, branchVelocities)
												line_allVs = c(line_allVs, branchVelocities)
												line_maxVs[1,n] = max(branchVelocities)
												line_meanVs[1,n] = mean(branchVelocities)
												line_medianVs[1,n] = median(branchVelocities)
											}
										line_allVs = as.matrix(t(line_allVs))
										row.names(line_allVs) = unlist(strsplit(file,"\\/"))[length(unlist(strsplit(file,"\\/")))]
										row.names(line_maxVs) = unlist(strsplit(file,"\\/"))[length(unlist(strsplit(file,"\\/")))]
										row.names(line_meanVs) = unlist(strsplit(file,"\\/"))[length(unlist(strsplit(file,"\\/")))]
										row.names(line_medianVs) = unlist(strsplit(file,"\\/"))[length(unlist(strsplit(file,"\\/")))]
										allBranchVelocities[[l]] = line_allVs
										maxBranchVelocities = rbind(maxBranchVelocities, line_maxVs)
										meanBranchVelocities = rbind(meanBranchVelocities, line_meanVs)
										medianBranchVelocities = rbind(medianBranchVelocities, line_medianVs)
									}
							}
					}
				HPD = round(HDInterval::hdi(allVs)[1:2],1)
				cat("\tMean branch velocities (",clades[i],", <=",cutoffs[m]," years) = ",round(mean(allVs),1),", 95% HPD = [",HPD[1],"-",HPD[2],"]\n",sep="")
				if (writingFiles)
					{
						sink(file=paste0(directories1[i],"_c",cutoffs[m],"_allBranchVelocities.txt"))
						for (l in 1:length(files))
							{
								cat(row.names(allBranchVelocities[[l]])," ",sep="")
								cat(allBranchVelocities[[l]])
								cat("\n")
							}
						sink(NULL)						
						write.table(maxBranchVelocities, paste0(directories1[i],"_c",cutoffs[m],"_maxBranchVelocities.csv"), sep=";", col.names=F)
						write.table(meanBranchVelocities, paste0(directories1[i],"_c",cutoffs[m],"_meanBranchVelocities.csv"), sep=";", col.names=F)
						write.table(medianBranchVelocities, paste0(directories1[i],"_c",cutoffs[m],"_medianBranchVelocities.csv"), sep=";", col.names=F)
					}
			}
	}

# 6. Visualising selected phylogenetic trees for SC1 and SC2

	# Preliminary step: re-organising the visualisation of both selected MCC trees ("increasing node order" option in FigTree)

logTransformation1 = FALSE; logTransformation2 = TRUE
directory2 = "humPangLocSI"; directory2 = "noHumPanLoc"
for (i in 1:length(clades))
	{
		if (i == 1) { hosts = c(); mccs = c() }
		mcc = read.csv(paste0(analysis,"/",clades[i],"/humPangLocSI/transPoW/",clades[i],"_",selected_NRRs[i],"_1.csv"), head=T)
		hosts = c(hosts, unique(mcc[,"host"])); mccs = rbind(mccs, mcc)
	}
hosts = unique(hosts); hosts = hosts[order(hosts)]; hosts = hosts[!is.na(hosts)]; counts = rep(NA, length(hosts))
for (i in 1:length(hosts)) counts[i] = sum(mccs[,"host"]==hosts[i], na.rm=T)
hosts = hosts[order(counts, hosts, decreasing=T)]
for (i in 1:length(clades))
	{
		tree = readAnnotatedNexus(paste0(analysis,"/",clades[i],"/",directory2,"/transPoW/",clades[i],"_",selected_NRRs[i],".tree"))
		mcc = read.csv(paste0(analysis,"/",clades[i],"/",directory2,"/transPoW/",clades[i],"_",selected_NRRs[i],"_1.csv"), head=T)
		mostRecentSamplingDatum = mostRecentSamplingDates[i]; correctedBranches = c()
		for (j in 1:length(tree$edge.length))
			{
				if (tree$edge.length[j] < 0)
					{
						buffer = tree$edge.length[j]; tree$edge.length[j] = 0
						branchToCorrect = which(tree$edge[,2]==tree$edge[j,1])
						if (!branchToCorrect%in%correctedBranches)
							{
								tree$edge.length[branchToCorrect] = tree$edge.length[branchToCorrect]+buffer
								correctedBranches = c(correctedBranches, branchToCorrect)
							}
					}
			}
		buffer = tree
		if (logTransformation1 == TRUE)
			{
				if (tree$root.annotation$`height_95%_HPD`[[1]] > 100)
					{
						buffer$root.annotation$`height_95%_HPD`[[1]] = 100+log(tree$root.annotation$`height_95%_HPD`[[1]]-100)
					}
				if (tree$root.annotation$`height_95%_HPD`[[2]] > 100)
					{
						buffer$root.annotation$`height_95%_HPD`[[2]] = 100+log(tree$root.annotation$`height_95%_HPD`[[2]]-100)
					}
			}
		if (logTransformation2 == TRUE)
			{
				buffer$root.annotation$`height_95%_HPD`[[1]] = log(tree$root.annotation$`height_95%_HPD`[[1]]+1)
				buffer$root.annotation$`height_95%_HPD`[[2]] = log(tree$root.annotation$`height_95%_HPD`[[2]]+1)
			}
		for (j in 1:dim(tree$edge)[1])
			{
				nodeAge1 = max(nodeHeights(tree))-nodeHeights(tree)[j,1]
				nodeAge2 = max(nodeHeights(tree))-nodeHeights(tree)[j,2]
				if (logTransformation2 == TRUE)
					{
						nodeAge1 = max(nodeHeights(tree))+1-nodeHeights(tree)[j,1]
						nodeAge2 = max(nodeHeights(tree))+1-nodeHeights(tree)[j,2]
					}
				if (logTransformation1 == TRUE)
					{
						if (nodeAge1 > 100) nodeAge1 = 100+log(nodeAge1-100)
						if (nodeAge2 > 100) nodeAge2 = 100+log(nodeAge2-100)
					}
				if (logTransformation2 == TRUE)
					{
						nodeAge1 = log(nodeAge1); nodeAge2 = log(nodeAge2)
					}
				buffer$edge.length[j] = nodeAge1-nodeAge2
				if (length(tree$annotations[[j]]$`height_95%_HPD`) > 1)
					{
						if (logTransformation1 == TRUE)
							{
								if (tree$annotations[[j]]$`height_95%_HPD`[[1]] > 100)
									{
										buffer$annotations[[j]]$`height_95%_HPD`[[1]] = 100+log(tree$annotations[[j]]$`height_95%_HPD`[[1]]-100)
									}
								if (tree$annotations[[j]]$`height_95%_HPD`[[2]] > 100)
									{
										buffer$annotations[[j]]$`height_95%_HPD`[[2]] = 100+log(tree$annotations[[j]]$`height_95%_HPD`[[2]]-100)
									}
							}
						if (logTransformation2 == TRUE)
							{
								buffer$annotations[[j]]$`height_95%_HPD`[[1]] = log(tree$annotations[[j]]$`height_95%_HPD`[[1]]+1)
								buffer$annotations[[j]]$`height_95%_HPD`[[2]] = log(tree$annotations[[j]]$`height_95%_HPD`[[2]]+1)
							}
					}	else	{
						# print(j)
					}
			}
		tree = buffer; rootHeight = max(nodeHeights(tree))
		root_time = mostRecentSamplingDatum-rootHeight; tree$tip.label = gsub("'","",tree$tip.label)
		minYear = mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]]; maxYear = mostRecentSamplingDatum
		if (clades[i] == "SC1") pdf(paste0(clades[i],"_",selected_NRRs[i],"_panel1.pdf"), width=4.0, height=4.0)
		if (clades[i] == "SC2") pdf(paste0(clades[i],"_",selected_NRRs[i],"_panel1.pdf"), width=4.0, height=2.0)
		par(mar=c(0.7,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.1,0), lwd=0.2, bty="o", col="gray30"); plottingRootNode = FALSE
		plot(tree, show.tip.label=F, show.node.label=F, edge.width=0.5, cex=0.6, align.tip.label=3, 
			 x.lim=c(minYear-(maxYear-max(nodeHeights(tree))), max(nodeHeights(tree))), col="gray30", edge.color="gray30")
		if (logTransformation1 == TRUE)
			{
				abline(v=(maxYear-75-root_time), lwd=0.3, lty=2)
				abline(v=(maxYear-50-root_time), lwd=0.3, lty=2)
			}
		if (logTransformation2 == TRUE)
			{
				abline(v=(maxYear-log(1+75)-root_time), lwd=0.3, lty=2)
				abline(v=(maxYear-log(1+50)-root_time), lwd=0.3, lty=2)	
				# abline(v=(maxYear-log(1+21.534)-root_time), lwd=0.3, lty=2) # 2000 (test)		
			}
		tree_obj = get("last_plot.phylo", envir=.PlotPhyloEnv); plottingRootBar = FALSE
		for (j in 1:dim(tree$edge)[1])
			{
				endYear = root_time+nodeHeights(tree)[j,2]
				endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
				endYear_colour = cols2[endYear_index]
				if ((tree$edge[j,2]%in%tree$edge[,1])&(length(tree$annotations[[j]]$`height_95%_HPD`) > 1))
					{
						x1 = (mostRecentSamplingDatum-tree$annotations[[j]]$`height_95%_HPD`[[2]])-root_time
						x2 = (mostRecentSamplingDatum-tree$annotations[[j]]$`height_95%_HPD`[[1]])-root_time
						lines(x=c(x1,x2), y=rep(tree_obj$yy[tree$edge[j,2]],2), lwd=2, lend=0, col=paste0(endYear_colour,"40"))
					}
				if ((plottingRootBar == FALSE)&&(!tree$edge[j,1]%in%tree$edge[,2]))
					{
						endYear = root_time+nodeHeights(tree)[j,1]
						endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
						endYear_colour = cols2[endYear_index]
						x1 = (mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]])-root_time
						x2 = (mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[1]])-root_time
						lines(x=c(x1,x2), y=rep(tree_obj$yy[tree$edge[j,1]],2), lwd=2, lend=0, col=paste0(endYear_colour,"40"))
						plottingRootBar = TRUE
					}				
			}
		for (j in 1:dim(tree$edge)[1])
			{
				endYear = root_time+nodeHeights(tree)[j,2]
				endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
				endYear_colour = cols2[endYear_index]
				if ((tree$edge[j,2]%in%tree$edge[,1])&&(tree$annotations[[j]]$posterior >= 0.95))
					{
						nodelabels(node=tree$edge[j,2], pch=16, cex=0.40, col=endYear_colour)
						nodelabels(node=tree$edge[j,2], pch=1, cex=0.40, col="gray30", lwd=0.2)
					}
				if (!tree$edge[j,2]%in%tree$edge[,1])
					{
						nodelabels(node=tree$edge[j,2], pch=15, cex=0.30, col=endYear_colour)
						nodelabels(node=tree$edge[j,2], pch=0, cex=0.30, col="gray30", lwd=0.2)
					}
				if ((plottingRootNode == FALSE)&&(!tree$edge[j,1]%in%tree$edge[,2]))
					{
						endYear = root_time+nodeHeights(tree)[j,1]
						endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
						endYear_colour = cols2[endYear_index]; plottingRootNode = TRUE		
						nodelabels(node=tree$edge[j,1], pch=16, cex=0.40, col=endYear_colour)
						nodelabels(node=tree$edge[j,1], pch=1, cex=0.40, col="gray30", lwd=0.2)
					}
			}
		for (j in 1:dim(tree$edge)[1])
			{
				if (!tree$edge[j,2]%in%tree$edge[,1])
					{
						index1 = which(mcc[,"tipLabel"]==tree$tip.label[tree$edge[j,2]])
						index2 = as.character(which(hosts==mcc[index1,"host"]))
						if (nchar(index2) == 1) index2 = paste0("     ",index2)
						if (nchar(index2) == 2) index2 = paste0("       ",index2)
						nodelabels(index2, node=tree$edge[j,2], cex=0.30, col="gray30", frame="none")
					}
			}
		selectedDates = c(-50000, -10000, -3000, 0, 1000, 1500, 1900, 2000, 2020); selectedLabels = selectedDates
		for (j in 1:length(selectedDates))
			{
				if (logTransformation1 == TRUE)
					{
						if ((mostRecentSamplingDatum-selectedDates[j]) > 100)
							{
								selectedDates[j] = mostRecentSamplingDatum-100-log((mostRecentSamplingDatum-100)-selectedDates[j])
							}
					}
				if (logTransformation2 == TRUE)
					{
						selectedDates[j] = mostRecentSamplingDatum-log(mostRecentSamplingDatum+1-selectedDates[j])
					}
			}
		selectedDates = c(minYear, selectedDates, maxYear); selectedLabels = c("", selectedLabels, "")
		axis(lwd=0.3, at=selectedDates-root_time, labels=selectedLabels, cex.axis=0.50, mgp=c(0,-0.20,-0.3), lwd.tick=0.3, 
			 col.lab="gray30", col="gray30", tck=-0.010, side=1)
		dev.off()
	}
cutoffs = c(75,50)
for (i in 1:length(clades)) # to plot the WLDV estimates based on the untransformed time-scaled trees
	{
		wldv_list1 = list()
		minX = 9999; maxX = -9999; minY = 9999; maxY = -9999
		for (j in 1:length(cutoffs))
			{
				wldv_list2 = list(); wldv_values = c()
				directories = list.files(paste0(analysis,"/",clades[i],"/",directory2,"/timeTrees/"))
				directories = directories[which((grepl(clades[i],directories))&(grepl(paste0("_s",cutoffs[j]),directories)))]
				for (k in 1:length(directories))
					{
						files = list.files(paste0(analysis,"/",clades[i],"/",directory2,"/timeTrees/",directories[k],"/"))
						file = files[which(grepl("dispersal_statistics",files))]
						tab = read.table(paste0(analysis,"/",clades[i],"/",directory2,"/timeTrees/",directories[k],"/",file), head=T)
						wldv_list2[[k]] = tab[,"weighted_branch_dispersal_velocity"]
						wldv_density = density(wldv_list2[[k]])
						if (minX > min(wldv_list2[[k]])) minX = min(wldv_list2[[k]])
						if (maxX < max(wldv_list2[[k]])) maxX = max(wldv_list2[[k]])
						if (minY > min(wldv_density$y)) minY = min(wldv_density$y)
						if (maxY < max(wldv_density$y)) maxY = max(wldv_density$y)
						wldv_values = c(wldv_values, wldv_list2[[k]])
					}
				wldv_list1[[j]] = wldv_list2; HPD = round(HDInterval::hdi(wldv_values)[1:2],1)
				cat("\tWLDV (",clades[i],", <=",cutoffs[j]," years) = ",round(median(wldv_values),1),", 95% HPD = [",HPD[1],"-",HPD[2],"]\n",sep="")
			}
		pdf(paste0(clades[i],"_",selected_NRRs[i],"_panel2_a.pdf"), width=3.5, height=2.2)
		par(mgp=c(0,0,0), oma=c(0,0,0,0), mar=c(2.1,2.3,0.5,0.5), lwd=0.3, col="gray30")
		for (j in 1:length(wldv_list1))
			{
				if (cutoffs[j] == 75) LTY = 1
				if (cutoffs[j] == 50) LTY = 2
				for (k in 1:length(wldv_list1[[j]]))
					{
						if ((j == 1)&(k == 1))
							{
								plot(density(wldv_list1[[j]][[k]]), xlim=c(0,50), ylim=c(minY,maxY), col=NA, axes=F, ann=F, frame=F)
							}
						if (cutoffs[j] == 75) polygon(density(wldv_list1[[j]][[k]]), lwd=0.01, border=NA, col=rgb(77,77,77,30,maxColorValue=255))
						if (cutoffs[j] == 50) lines(density(wldv_list1[[j]][[k]]), lwd=0.20, lty=LTY, col="gray30")
					}
			}
		axis(side=1, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,-0.14,0), lwd=0.3, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=2, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,0.18,0), lwd=0.3, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		title(xlab="weighted lineage dispersal velocity (WLDV)", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
		title(ylab="density", cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
		dev.off()
	}
cutoffs = c(75,50)
for (i in 1:length(clades)) # to report the WDC estimates based on the untransformed time-scaled trees
	{
		wdcs_list1 = list()
		minX = 9999; maxX = -9999; minY = 9999; maxY = -9999
		for (j in 1:length(cutoffs))
			{
				wdcs_list2 = list(); wdcs_values = c()
				directories = list.files(paste0(analysis,"/",clades[i],"/",directory2,"/timeTrees/"))
				directories = directories[which((grepl(clades[i],directories))&(grepl(paste0("_s",cutoffs[j]),directories)))]
				for (k in 1:length(directories))
					{
						files = list.files(paste0(analysis,"/",clades[i],"/",directory2,"/timeTrees/",directories[k],"/"))
						file = files[which(grepl("dispersal_statistics",files))]
						tab = read.table(paste0(analysis,"/",clades[i],"/",directory2,"/timeTrees/",directories[k],"/",file), head=T)
						wdcs_list2[[k]] = tab[,"weighted_diffusion_coefficient"]
						wdcs_density = density(wdcs_list2[[k]])
						if (minX > min(wdcs_list2[[k]])) minX = min(wdcs_list2[[k]])
						if (maxX < max(wdcs_list2[[k]])) maxX = max(wdcs_list2[[k]])
						if (minY > min(wdcs_density$y)) minY = min(wdcs_density$y)
						if (maxY < max(wdcs_density$y)) maxY = max(wdcs_density$y)
						wdcs_values = c(wdcs_values, wdcs_list2[[k]])
					}
				wdcs_list1[[j]] = wdcs_list2; HPD = round(HDInterval::hdi(wdcs_values)[1:2],1)
				cat("\tWDC (",clades[i],", >",cutoffs[j]," years) = ",round(median(wdcs_values),1),", 95% HPD = [",HPD[1],"-",HPD[2],"]\n",sep="")
			}
	}
cutoffs = c(75,50)
for (i in 1:length(clades)) # to plot the WDC estimates based on the untransformed time-scaled trees
	{
		wdcs_list1 = list(); all_values = list()
		minX = 9999; maxX = -9999; minY = 9999; maxY = -9999
		for (j in 1:length(cutoffs))
			{
				wdcs_list2 = list(); wdcs_values = c()
				directories = list.files(paste0(analysis,"/",clades[i],"/",directory2,"/timeTrees/"))
				directories = directories[which((grepl(clades[i],directories))&(grepl(paste0("_s",cutoffs[j]),directories)))]
				for (k in 1:length(directories))
					{
						files = list.files(paste0(analysis,"/",clades[i],"/",directory2,"/timeTrees/",directories[k],"/"))
						file = files[which(grepl("dispersal_statistics",files))]
						tab = read.table(paste0(analysis,"/",clades[i],"/",directory2,"/timeTrees/",directories[k],"/",file), head=T)
						wdcs_list2[[k]] = tab[,"weighted_diffusion_coefficient"]
						wdcs_density = density(wdcs_list2[[k]])
						if (minX > min(wdcs_list2[[k]])) minX = min(wdcs_list2[[k]])
						if (maxX < max(wdcs_list2[[k]])) maxX = max(wdcs_list2[[k]])
						if (minY > min(wdcs_density$y)) minY = min(wdcs_density$y)
						if (maxY < max(wdcs_density$y)) maxY = max(wdcs_density$y)
						wdcs_values = c(wdcs_values, wdcs_list2[[k]])
					}
				wdcs_list1[[j]] = wdcs_list2; all_values[[j]] = wdcs_values; HPD = round(HDInterval::hdi(wdcs_values)[1:2],1)
				cat("\tWDC (",clades[i],", <=",cutoffs[j]," years) = ",round(median(wdcs_values),1),", 95% HPD = [",HPD[1],"-",HPD[2],"]\n",sep="")
			}
		pdf(paste0(clades[i],"_",selected_NRRs[i],"_panel2_b1.pdf"), width=3.5, height=2.2)
		par(mgp=c(0,0,0), oma=c(0,0,0,0), mar=c(2.1,2.3,0.5,0.5), lwd=0.3, col="gray30")
		for (j in 1:length(wdcs_list1))
			{
				if (cutoffs[j] == 75) LTY = 1
				if (cutoffs[j] == 50) LTY = 2
				for (k in 1:length(wdcs_list1[[j]]))
					{
						if ((j == 1)&(k == 1))
							{
								plot(density(wdcs_list1[[j]][[k]]), xlim=c(0,13000), ylim=c(minY,0.0022), col=NA, axes=F, ann=F, frame=F)
							}
						if (cutoffs[j] == 75) polygon(density(wdcs_list1[[j]][[k]]), lwd=0.01, border=NA, col=rgb(77,77,77,30,maxColorValue=255))
						if (cutoffs[j] == 50) lines(density(wdcs_list1[[j]][[k]]), lwd=0.20, lty=LTY, col="gray30")
					}
			}
		axis(side=1, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,-0.14,0), lwd=0.3, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=2, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,0.18,0), lwd=0.3, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		title(xlab="weighted diffusion coefficient", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
		title(ylab="density", cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
		dev.off()
		pdf(paste0(clades[i],"_",selected_NRRs[i],"_panel2_b2.pdf"), width=3.5, height=2.2)
		par(mgp=c(0,0,0), oma=c(0,0,0,0), mar=c(2.1,2.3,0.5,0.5), lwd=0.3, col="gray30")
		for (j in 1:length(wdcs_list1))
			{
				if (cutoffs[j] == 75) LTY = 1
				if (cutoffs[j] == 50) LTY = 2
				if (j == 1)
						{
							plot(density(all_values[[j]]), xlim=c(0,10000), ylim=c(minY,0.0006), col=NA, axes=F, ann=F, frame=F)
						}
				if (cutoffs[j] == 75) polygon(density(all_values[[j]]), lwd=0.01, border=NA, col=rgb(77,77,77,30,maxColorValue=255))
				if (cutoffs[j] == 50) lines(density(all_values[[j]]), lwd=0.20, lty=LTY, col="gray30")
			}
		axis(side=1, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,-0.14,0), lwd=0.3, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=2, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,0.18,0), lwd=0.3, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		title(xlab="weighted diffusion coefficient", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
		title(ylab="density", cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
		dev.off()
	}
cutoffs = c(75,50)
for (i in 1:length(clades)) # to plot the WLDV estimates provided by TimeSlicer (not used anymore)
	{
		wldv_list1 = list()
		minX = 9999; maxX = -9999; minY = 9999; maxY = -9999
		for (j in 1:length(cutoffs))
			{
				wldv_list2 = list(); wldv_values = c()
				files = list.files(paste0(analysis,"/",clades[i],"/",directory2,"/statistics/"))
				files = files[which((grepl("WLDV",files))&(grepl(paste0("_",cutoffs[j]),files)))]
				for (k in 1:length(files))
					{
						tab = read.table(paste0(analysis,"/",clades[i],"/",directory2,"/statistics/",files[k]), head=F)
						wldv_list2[[k]] = tab[,1]; wldv_density = density(wldv_list2[[k]])
						if (minX > min(wldv_list2[[k]])) minX = min(wldv_list2[[k]])
						if (maxX < max(wldv_list2[[k]])) maxX = max(wldv_list2[[k]])
						if (minY > min(wldv_density$y)) minY = min(wldv_density$y)
						if (maxY < max(wldv_density$y)) maxY = max(wldv_density$y)
						wldv_values = c(wldv_values, wldv_list2[[k]])
					}
				wldv_list1[[j]] = wldv_list2; HPD = round(HDInterval::hdi(wldv_values)[1:2],1)
				cat("\tWLDV (",clades[i],", <=",cutoffs[j]," years) = ",round(median(wldv_values),1),", 95% HPD = [",HPD[1],"-",HPD[2],"]\n",sep="")
			}
		pdf(paste0(clades[i],"_",selected_NRRs[i],"_panel2c.pdf"), width=3.5, height=2.2)
		par(mgp=c(0,0,0), oma=c(0,0,0,0), mar=c(2.1,2.3,0.5,0.5), lwd=0.3, col="gray30")
		for (j in 1:length(wldv_list1))
			{
				if (cutoffs[j] == 75) LTY = 1
				if (cutoffs[j] == 50) LTY = 2
				for (k in 1:length(wldv_list1[[j]]))
					{
						if ((j == 1)&(k == 1))
							{
								plot(density(wldv_list1[[j]][[k]]), xlim=c(0,40), ylim=c(minY,maxY), col=NA, axes=F, ann=F, frame=F)
							}
						if (cutoffs[j] == 75) polygon(density(wldv_list1[[j]][[k]]), lwd=0.01, border=NA, col=rgb(77,77,77,30,maxColorValue=255))
						if (cutoffs[j] == 50) lines(density(wldv_list1[[j]][[k]]), lwd=0.20, lty=LTY, col="gray30")
					}
			}
		axis(side=1, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,-0.14,0), lwd=0.3, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=2, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,0.18,0), lwd=0.3, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		title(xlab="weighted lineage dispersal velocity (WLDV)", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
		title(ylab="density", cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
		dev.off()
	}
for (i in 1:length(clades)) # to plot the WDC estimates provided by TimeSlicer (not used anymore)
	{
		wdcs_list1 = list()
		minX = 9999; maxX = -9999; minY = 9999; maxY = -9999
		for (j in 1:length(cutoffs))
			{
				wdcs_list2 = list(); wdcs_values = c()
				files = list.files(paste0(analysis,"/",clades[i],"/",directory2,"/statistics/"))
				files = files[which((grepl("WDC",files))&(grepl(paste0("_",cutoffs[j]),files)))]
				for (k in 1:length(files))
					{
						tab = read.table(paste0(analysis,"/",clades[i],"/",directory2,"/statistics/",files[k]), head=F)
						wdcs_list2[[k]] = tab[,1]; wdcs_density = density(wdcs_list2[[k]])
						if (minX > min(wdcs_list2[[k]])) minX = min(wdcs_list2[[k]])
						if (maxX < max(wdcs_list2[[k]])) maxX = max(wdcs_list2[[k]])
						if (minY > min(wdcs_density$y)) minY = min(wdcs_density$y)
						if (maxY < max(wdcs_density$y)) maxY = max(wdcs_density$y)
						wdcs_values = c(wdcs_values, wdcs_list2[[k]])
					}
				wdcs_list1[[j]] = wdcs_list2; HPD = round(HDInterval::hdi(wdcs_values)[1:2],1)
				cat("\tWDC (",clades[i],", <=",cutoffs[j]," years) = ",round(median(wdcs_values),1),", 95% HPD = [",HPD[1],"-",HPD[2],"]\n",sep="")
			}
		pdf(paste0(clades[i],"_",selected_NRRs[i],"_panel2d.pdf"), width=3.5, height=2.2)
		par(mgp=c(0,0,0), oma=c(0,0,0,0), mar=c(2.1,2.3,0.5,0.5), lwd=0.3, col="gray30")
		for (j in 1:length(wdcs_list1))
			{
				if (cutoffs[j] == 75) LTY = 1
				if (cutoffs[j] == 50) LTY = 2
				for (k in 1:length(wdcs_list1[[j]]))
					{
						if ((j == 1)&(k == 1))
							{
								plot(density(wdcs_list1[[j]][[k]]), xlim=c(0,13000), ylim=c(minY,0.0015), col=NA, axes=F, ann=F, frame=F)
							}
						if (cutoffs[j] == 75) polygon(density(wdcs_list1[[j]][[k]]), lwd=0.01, border=NA, col=rgb(77,77,77,30,maxColorValue=255))
						if (cutoffs[j] == 50) lines(density(wdcs_list1[[j]][[k]]), lwd=0.20, lty=LTY, col="gray30")
					}
			}
		axis(side=1, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,-0.14,0), lwd=0.3, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=2, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,0.18,0), lwd=0.3, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		title(xlab="weighted diffusion coefficient (km2/year)", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
		title(ylab="density", cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
		dev.off()
	}

# 7. Visualising selected continuous phylogeographic reconstructions

logTransformation1 = FALSE; logTransformation2 = TRUE
directory2 = "humPangLocSI"; directory2 = "noHumPanLoc"
checkingSamplingCoordinates = FALSE
if (checkingSamplingCoordinates)
	{
		for (i in 1:length(clades))
			{
				files = list.files(paste0(analysis,"/",clades[i],"/",directory2,"/strictClock/"))
				files = files[which(grepl("_ext1",files))]
				NRRs = gsub("SC1_NRR","",gsub("SC2_NRR","",gsub("_b_ext1","",files)))
				for (j in 1:length(NRRs))
					{
						pdf(paste0("TEMP_",clades[i],"_NRR",NRRs[j],".pdf"))
						plot(background, col=cols1, axes=F, ann=F, box=F, legend=F)
						for (k in 1:100)
							{
								tab = read.csv(paste0(analysis,"/",clades[i],"/",directory2,"/strictClock/",clades[i],"_NRR",NRRs[j],"_b_ext1/TreeExtractions_",k,".csv"), head=T)
								points(tab[which(!tab[,"node2"]%in%tab[,"node1"]),c("endLon","endLat")], cex=0.2, col="gray30")
							}
						dev.off()
					}
			}
	}
for (i in 1:length(clades))
	{
		if (i == 1) { hosts = c(); mccs = c() }
		mcc = read.csv(paste0(analysis,"/",clades[i],"/humPangLocSI/transPoW/",clades[i],"_",selected_NRRs[i],"_1.csv"), head=T)
		hosts = c(hosts, unique(mcc[,"host"])); mccs = rbind(mccs, mcc)
	}
hosts = unique(hosts); hosts = hosts[order(hosts)]; hosts = hosts[!is.na(hosts)]; counts = rep(NA, length(hosts))
for (i in 1:length(hosts)) counts[i] = sum(mccs[,"host"]==hosts[i], na.rm=T)
hosts = hosts[order(counts, hosts, decreasing=T)]
croppingPolygons = FALSE; cutoffs = c(99999); polygons = list()
for (i in 1:length(clades))
	{
		tree = readAnnotatedNexus(paste0(analysis,"/",clades[i],"/noHumPanLoc/transPoW/",clades[i],"_",selected_NRRs[i],".tree")) # to get a fixed time-scale
		mostRecentSamplingDatum = mostRecentSamplingDates[i]; correctedBranches = c()
		for (j in 1:length(tree$edge.length))
			{
				if (tree$edge.length[j] < 0)
					{
						buffer = tree$edge.length[j]; tree$edge.length[j] = 0
						branchToCorrect = which(tree$edge[,2]==tree$edge[j,1])
						if (!branchToCorrect%in%correctedBranches)
							{
								tree$edge.length[branchToCorrect] = tree$edge.length[branchToCorrect]+buffer
								correctedBranches = c(correctedBranches, branchToCorrect)
							}
					}
			}
		buffer = tree
		if (logTransformation1 == TRUE)
			{
				if (tree$root.annotation$`height_95%_HPD`[[1]] > 100)
					{
						buffer$root.annotation$`height_95%_HPD`[[1]] = 100+log(tree$root.annotation$`height_95%_HPD`[[1]]-100)
					}
				if (tree$root.annotation$`height_95%_HPD`[[2]] > 100)
					{
						buffer$root.annotation$`height_95%_HPD`[[2]] = 100+log(tree$root.annotation$`height_95%_HPD`[[2]]-100)
					}
			}
		if (logTransformation2 == TRUE)
			{
				buffer$root.annotation$`height_95%_HPD`[[1]] = log(tree$root.annotation$`height_95%_HPD`[[1]]+1)
				buffer$root.annotation$`height_95%_HPD`[[2]] = log(tree$root.annotation$`height_95%_HPD`[[2]]+1)
			}
		for (j in 1:dim(tree$edge)[1])
			{
				nodeAge1 = max(nodeHeights(tree))-nodeHeights(tree)[j,1]
				nodeAge2 = max(nodeHeights(tree))-nodeHeights(tree)[j,2]
				if (logTransformation2 == TRUE)
					{
						nodeAge1 = max(nodeHeights(tree))+1-nodeHeights(tree)[j,1]
						nodeAge2 = max(nodeHeights(tree))+1-nodeHeights(tree)[j,2]
					}
				if (logTransformation1 == TRUE)
					{
						if (nodeAge1 > 100) nodeAge1 = 100+log(nodeAge1-100)
						if (nodeAge2 > 100) nodeAge2 = 100+log(nodeAge2-100)
					}
				if (logTransformation2 == TRUE)
					{
						nodeAge1 = log(nodeAge1); nodeAge2 = log(nodeAge2)
					}
				buffer$edge.length[j] = nodeAge1-nodeAge2
			}
		tree = buffer; rootHeight = max(nodeHeights(tree))
		root_time = mostRecentSamplingDatum-rootHeight; tree$tip.label = gsub("'","",tree$tip.label)
		minYear = mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]]; maxYear = mostRecentSamplingDatum
		mcc = read.csv(paste0(analysis,"/",clades[i],"/",directory2,"/transPoW/",clades[i],"_",selected_NRRs[i],"_2.csv"), head=T)		
		endYears_colours = rep(NA, dim(mcc)[1])
		for (j in 1:length(endYears_colours))
			{
				if (logTransformation1 == TRUE) nodeAge = maxYear-mcc[j,"endYear"]
				if (logTransformation2 == TRUE) nodeAge = maxYear+1-mcc[j,"endYear"]
				if ((logTransformation1 == TRUE)&(nodeAge > 100)) nodeAge = 100+log(nodeAge-100)
				if (logTransformation2 == TRUE) nodeAge = log(nodeAge)
				endYearMod = maxYear-nodeAge
				endYear_index = (((endYearMod-minYear)/(maxYear-minYear))*100)+1
				endYears_colours[j] = cols2[endYear_index]
			}
		localTreesDirectory = paste0(analysis,"/",clades[i],"/",directory2,"/transPoW/",clades[i],"_",selected_NRRs[i],"_ext2")
		nberOfExtractionFiles = length(list.files(localTreesDirectory)); prob = 0.80; precision = 50; startDatum = maxYear-1000
		polygons[[i]] = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
		polygons_colours = rep(NA, length(polygons[[i]]))
		for (j in 1:length(polygons[[i]]))
			{
				date = as.numeric(names(polygons[[i]][[j]])); polyAge = NULL
				if (logTransformation1 == TRUE) polyAge = maxYear-date
				if (logTransformation2 == TRUE) polyAge = maxYear+1-date
				if ((logTransformation1 == TRUE)&(nodeAge > 100)) polyAge = 100+log(polyAge-100)
				if (logTransformation2 == TRUE) polyAge = log(polyAge)
				endYearMod = maxYear-polyAge
				endYear_index = (((endYearMod-minYear)/(maxYear-minYear))*100)+1
				polygons_colours[j] = paste0(cols2[endYear_index],"20")
			}
		j = 1
		for (j in 1:length(cutoffs))
			{
				pdf(paste0(clades[i],"_",selected_NRRs[i],"_panel3.pdf"), width=7, height=4.8)
				par(oma=c(0,0,0,0), mar=c(0.0,2.5,0.0,0.0), lwd=0.2, col="gray30")
				if (plottingElevation == FALSE)
					{
						plot(background, col=cols1, axes=F, ann=F, box=F, legend=F)
					}	else	{
						plot(elevation_mod, col=cols3, axes=F, ann=F, box=F, legend=F)
					}
				lines(borders, col="white", lwd=0.3)
				if (plottingElevation == FALSE)
					{
						for (k in 1:length(polygons[[i]]))
							{
								if (as.numeric(names(polygons[[i]][[k]])) > (mostRecentSamplingDatum-cutoffs[j]))
									{
										for (l in 1:length(polygons[[i]][[k]]@polygons))
											{
												polygons[[i]][[k]]@polygons[[l]] = checkPolygonsHoles(polygons[[i]][[k]]@polygons[[l]])
											}
										pol = polygons[[i]][[k]]; crs(pol) = crs(background)
										if (croppingPolygons == TRUE) pol = crop(pol, provinces)
										plot(pol, axes=F, col=polygons_colours[k], add=T, border=NA)
									}
							}
					}
				for (k in 1:dim(mcc)[1])
					{
						if (mcc[k,"startYear"] > (mostRecentSamplingDatum-cutoffs[j]))
							{
								curvedarrow(cbind(mcc[k,"startLon"],mcc[k,"startLat"]), cbind(mcc[k,"endLon"],mcc[k,"endLat"]), arr.length=0,
						  		  		    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
						  	}
					}
				for (k in dim(mcc)[1]:1)
					{
						if ((mcc[k,"endYear"] > (mostRecentSamplingDatum-cutoffs[j]))&(!mcc[k,"node1"]%in%mcc[,"node2"]))
							{
								startYears_index = (((mcc[k,"startYear"]-minYear)/(maxYear-minYear))*100)+1
								points(mcc[k,"startLon"], mcc[k,"startLat"], pch=16, col=cols2[startYears_index], cex=0.4)
								points(mcc[k,"startLon"], mcc[k,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
							}
						if (mcc[k,"endYear"] > (mostRecentSamplingDatum-cutoffs[j]))
							{
								if (mcc[k,"node2"]%in%mcc[,"node1"])
									{	
										points(mcc[k,"endLon"], mcc[k,"endLat"], pch=16, col=endYears_colours[k], cex=0.4)
										points(mcc[k,"endLon"], mcc[k,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
									}	else		{
										points(mcc[k,"endLon"], mcc[k,"endLat"], pch=23, bg=endYears_colours[k], cex=0.35, col=NA)
										points(mcc[k,"endLon"], mcc[k,"endLat"], pch=5, col="gray30", lwd=0.2, cex=0.35)
									}
							}
					}
				positions1 = rbind(mcc[,c("endLon","endLat")]); colnames(positions1) = c("x2","y2"); positions2 = c()
				for (k in dim(mcc)[1]:1)
					{
						if (mcc[k,"endYear"] > (mostRecentSamplingDatum-cutoffs[j]))
							{
								if (!mcc[k,"node2"]%in%mcc[,"node1"])
									{
										if (!is.na(mcc[k,"host"]))
											{
												index = which(hosts==mcc[k,"host"]); x1 = mcc[k,"endLon"]; y1 = mcc[k,"endLat"]
												if (nchar(as.character(index)) == 1)
													{
														x2s = c(x1+0.4, x1+0.4, x1-0.4, x1-0.4); y2s = c(y1+0.4, y1-0.4, y1+0.4, y1-0.4)
													}												
												if (nchar(as.character(index)) == 2)
													{
														x2s = c(x1+0.5, x1+0.5, x1-0.5, x1-0.5); y2s = c(y1+0.5, y1-0.5, y1+0.5, y1-0.5)
													}
												dS = rep(NA, length(x2s))
												for (l in 1:length(x2s))
													{
														dS[l] = min(sqrt(((x2s[l]-positions1[,1])^2)+((y2s[l]-positions1[,2])^2)))
													}
												x2 = x2s[which(dS==max(dS))[1]]; y2 = y2s[which(dS==max(dS))[1]]
												if (is.null(positions2))
													{
														text(x2, y2, as.character(index), cex=0.4, col="gray30")
														positions1 = rbind(positions1, cbind(x2, y2))
														positions2 = rbind(positions2, cbind(x1, y1, index))														
													}	else	{
														if (!index%in%positions2[,3])
															{
																text(x2, y2, as.character(index), cex=0.4, col="gray30")
																positions1 = rbind(positions1, cbind(x2, y2))
																positions2 = rbind(positions2, cbind(x1, y1, index))														
															}	else	{
																indices = which(positions2[,3]==index)
																d = min(sqrt(((x1-positions2[indices,1])^2)+((y1-positions2[indices,2])^2)))
																if (d > 0.5)
																	{
																		text(x2, y2, as.character(index), cex=0.4, col="gray30")
																		positions1 = rbind(positions1, cbind(x2, y2))
																		positions2 = rbind(positions2, cbind(x1, y1, index))
																	}
															}
													}
											}
									}
							}
					}
				highlightingZhoushanIsland = FALSE
				if (highlightingZhoushanIsland)
					{
						points(cbind(122.1417,30.1441), pch=16, col="red")
					}
				selectedDates = c(-10000, -3000, 0, 1000, 1500, 1900, 2000, 2020); selectedLabels = selectedDates
				for (j in 1:length(selectedDates))
					{
						if (logTransformation1 == TRUE)
							{
								if ((mostRecentSamplingDatum-selectedDates[j]) > 100)
									{
										selectedDates[j] = maxYear-100-log((mostRecentSamplingDatum-100)-selectedDates[j])
									}
							}
						if (logTransformation2 == TRUE)
							{
								selectedDates[j] = maxYear-log(maxYear+1-selectedDates[j])
							}
					}
				if (plottingElevation == TRUE)
					{
						plot(elevation_mod, legend.only=T, add=T, col=cols3, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.853,0.859,0.090,0.905),
	 						 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=F,
							 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, tck=-1.2, col.tick="gray30", col="gray30", col.axis="gray30", line=0, mgp=c(0,0.45,0)))						
					}	else	{
						rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minYear; rast[2] = maxYear
						plot(rast, legend.only=T, add=T, col=cols2, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.853,0.859,0.090,0.905),
			 				 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=F,
					 		 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, tck=-1.2, col.tick="gray30", col="gray30", col.axis="gray30", line=0, mgp=c(0,0.45,0), at=selectedDates, labels=selectedLabels))
					}
				dev.off()
			}
	}

# 8. Displaying and analysing the estimated position of the human ancestor

	# 8.1. Figures based on estimated HPD polygons (not used anymore)

usePreviousScript = FALSE
if (usePreviousScript)
	{
		logTransformation1 = FALSE; logTransformation2 = TRUE
		directory2 = "humPangLocSI"; directory2 = "noHumPanLoc"
		directory3 = "strictClock"; directory3 = "transPoW"
		human_samples = c("SC1_AY394995.1","SC2_MN908947.3")
		croppingPolygons = FALSE; cutoffs = c(99999)
		ancestor_pols_list = list(); mcc_times_list = list()
		for (i in 1:length(clades))
			{
				ancestor_pols = list(); mcc_times = list()
				human_sample = human_samples[which(clades==clades[i])]
				files = list.files(paste0(analysis,"/",clades[i],"/",directory2,"/",directory3))
				files = gsub(".csv","",files[which(grepl(".csv",files))])
				files = files[order(as.numeric(gsub("SC1_NRR","",gsub("SC2_NRR","",files))))]
				for (j in 1:length(files))
					{
						localTreesDirectory = paste0(analysis,"/",clades[i],"/",directory2,"/",directory3,"/",files[j],"_ext1")
						mcc = read.csv(paste0(analysis,"/",clades[i],"/",directory2,"/",directory3,"/",files[j],"_1.csv"), head=T)
						nberOfExtractionFiles = length(list.files(localTreesDirectory))
						positions = matrix(nrow=nberOfExtractionFiles, ncol=2)
						colnames(positions) = c("longitude","latitude")
						for (k in 1:nberOfExtractionFiles)
							{
								tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",k,".csv"), head=T)
								positions[k,"longitude"] = tab[which(tab[,"tipLabel"]==human_sample),"startLon"]
								positions[k,"latitude"] = tab[which(tab[,"tipLabel"]==human_sample),"startLat"]
							}
						mcc_times[[j]] = mcc[which(mcc[,"tipLabel"]==human_sample),"startYear"]
						H = Hpi(cbind(positions[,"longitude"],positions[,"latitude"]))
						kde = kde(cbind(positions[,"longitude"],positions[,"latitude"]), H=H, compute.cont=T, gridsize=c(1000,1000))
						prob = 0.80; contourLevel = contourLevels(kde, prob=(1-prob)); pols = list()
						contourLines = contourLines(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate, level=contourLevel)
						for (k in 1:length(contourLines)) pols[[k]] = Polygon(cbind(contourLines[[k]]$x,contourLines[[k]]$y))
						ps = Polygons(pols,1); contourPolygons = SpatialPolygons(list(ps)); # plot(contourPolygons, add=T)
						spdf = SpatialPolygonsDataFrame(contourPolygons, data.frame(ID=1:length(contourPolygons)))
						names(spdf) = mcc_times[[j]]; ancestor_pols[[j]] = spdf
					}
				ancestor_pols_list[[i]] = ancestor_pols; mcc_times_list[[i]] = mcc_times
			}
		for (i in 1:length(clades))
			{
				tree = readAnnotatedNexus(paste0(analysis,"/",clades[i],"/",directory2,"/",directory3,"/",clades[i],"_",selected_NRRs[i],".tree"))
				mostRecentSamplingDatum = mostRecentSamplingDates[i]; correctedBranches = c()
				for (j in 1:length(tree$edge.length))
					{
						if (tree$edge.length[j] < 0)
							{
								buffer = tree$edge.length[j]; tree$edge.length[j] = 0
								branchToCorrect = which(tree$edge[,2]==tree$edge[j,1])
								if (!branchToCorrect%in%correctedBranches)
									{
										tree$edge.length[branchToCorrect] = tree$edge.length[branchToCorrect]+buffer
										correctedBranches = c(correctedBranches, branchToCorrect)
									}
							}
					}
				buffer = tree
				if (logTransformation1 == TRUE)
					{
						if (tree$root.annotation$`height_95%_HPD`[[1]] > 100)
							{
								buffer$root.annotation$`height_95%_HPD`[[1]] = 100+log(tree$root.annotation$`height_95%_HPD`[[1]]-100)
							}
						if (tree$root.annotation$`height_95%_HPD`[[2]] > 100)
							{
								buffer$root.annotation$`height_95%_HPD`[[2]] = 100+log(tree$root.annotation$`height_95%_HPD`[[2]]-100)
							}
					}
				if (logTransformation2 == TRUE)
					{
						buffer$root.annotation$`height_95%_HPD`[[1]] = log(tree$root.annotation$`height_95%_HPD`[[1]]+1)
						buffer$root.annotation$`height_95%_HPD`[[2]] = log(tree$root.annotation$`height_95%_HPD`[[2]]+1)
					}
				tree = buffer; rootHeight = max(nodeHeights(tree))
				root_time = mostRecentSamplingDatum-rootHeight; tree$tip.label = gsub("'","",tree$tip.label)
				minYear = mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]]; maxYear = mostRecentSamplingDatum
				polygons = list(); c = 0
				for (j in 1:length(ancestor_pols_list[[i]]))
					{
						c = c+1; polygons[[c]] = ancestor_pols_list[[i]][[j]]
					}
				polygons_colours = rep(NA, length(polygons))
				for (j in 1:length(polygons))
					{
						date = as.numeric(names(polygons[[j]])); polyAge = maxYear-date
						if (logTransformation1 == TRUE) polyAge = maxYear-date
						if (logTransformation2 == TRUE) polyAge = maxYear+1-date
						if ((logTransformation1 == TRUE)&(polyAge > 100)) polyAge = 100+log(polyAge-100)
						if (logTransformation2 == TRUE) polyAge = log(polyAge)
						endYearMod = maxYear-polyAge
						endYear_index = (((endYearMod-minYear)/(maxYear-minYear))*100)+1
						polygons_colours[j] = paste0(cols2[endYear_index],"20")
					}
				for (j in 1:length(cutoffs))
					{
						pdf(paste0(clades[i],"_panel6.pdf"), width=7, height=4.8)
						par(oma=c(0,0,0,0), mar=c(0.0,2.5,0.0,0.0), lwd=0.2, col="gray30")
						plot(background, col=cols1, axes=F, ann=F, box=F, legend=F)
						# lines(borders, col="white", lwd=0.3)
						for (k in 1:length(polygons))
							{
								if (as.numeric(names(polygons[[k]])) > (mostRecentSamplingDatum-cutoffs[j]))
									{
										for (l in 1:length(polygons[[k]]@polygons))
											{
												polygons[[k]]@polygons[[l]] = checkPolygonsHoles(polygons[[k]]@polygons[[l]])
											}
										pol = polygons[[k]]; crs(pol) = crs(background)
										if (croppingPolygons == TRUE) pol = crop(pol, provinces)
										plot(pol, axes=F, col=polygons_colours[k], add=T, border=NA)
									}
							}
						rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minYear; rast[2] = mostRecentSamplingDatum
						plot(rast, legend.only=T, add=T, col=cols2, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.853,0.859,0.090,0.905),
			 				 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=F,
							 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, tck=-1.2, col.tick="gray30", col="gray30", col.axis="gray30", line=0, mgp=c(0,0.45,0)))
						points(cbind(113.2587,23.1414), pch=16, cex=0.6, col="gray30") # Guangzhou
						text(cbind(113.2587,24.1414), "G.", cex=0.4, col="gray30")
						points(cbind(114.3090,30.5980), pch=16, cex=0.6, col="gray30") # Wuhan
						text(cbind(114.3090,31.5980), "W.", cex=0.4, col="gray30")
						dev.off()
					}
			}
	}

	# 8.2. Figures based on recCA positions shared by Philippe

directory2 = "humPangLocSI"; directory2 = "noHumPanLoc"
probabilities = c(0.95,0.75,0.50); red_cols = list()
red_cols[[1]] = rgb(204,0,0,60,maxColorValue=255)
red_cols[[2]] = rgb(204,0,0,100,maxColorValue=255)
red_cols[[3]] = rgb(204,0,0,140,maxColorValue=255)
croppingPolygons = FALSE
for (i in 1:length(clades))
	{
		positions = read.table(paste0(analysis,"/",clades[i],"/",directory2,"/parentL.txt"), head=T)[,1:2]
		H = Hpi(cbind(positions[,"longitude"],positions[,"latitude"]), pilot="samse"); polygons = list()
		kde = kde(cbind(positions[,"longitude"],positions[,"latitude"]), H=H, compute.cont=T, gridsize=c(1000,1000))
		for (j in 1:length(probabilities))
			{
				prob = probabilities[j]; contourLevel = contourLevels(kde, prob=(1-prob)); pols = list()
				contourLines = contourLines(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate, level=contourLevel)
				for (k in 1:length(contourLines)) pols[[k]] = Polygon(cbind(contourLines[[k]]$x,contourLines[[k]]$y))
				ps = Polygons(pols,1); contourPolygons = SpatialPolygons(list(ps)); # plot(contourPolygons, add=T)
				polygons[[j]] = SpatialPolygonsDataFrame(contourPolygons, data.frame(ID=1:length(contourPolygons)))
			}
		pdf(paste0(clades[i],"_panel5.pdf"), width=7, height=4.8)
		par(oma=c(0,0,0,0), mar=c(0.0,2.5,0.0,0.0), lwd=0.2, col="gray30")
		plot(background, col=cols1, axes=F, ann=F, box=F, legend=F)
		lines(borders, col="white", lwd=0.3)
		for (j in 1:length(polygons))
			{
				for (k in 1:length(polygons[[j]]@polygons))
					{
						polygons[[j]]@polygons[[k]] = checkPolygonsHoles(polygons[[j]]@polygons[[k]])
					}
				pol = polygons[[j]]; crs(pol) = crs(background)
				if (croppingPolygons == TRUE) pol = crop(pol, provinces)
				plot(pol, axes=F, col=red_cols[[j]], add=T, border=NA)
			}
			
		points(cbind(113.2587,23.1414), pch=16, cex=0.6, col="gray30") # Guangzhou
		text(cbind(113.2587,24.1414), "G.", cex=0.4, col="gray30")
		points(cbind(114.3090,30.5980), pch=16, cex=0.6, col="gray30") # Wuhan
		text(cbind(114.3090,31.5980), "W.", cex=0.4, col="gray30")
		dev.off()
	}

	# 8.3. Analysing the distances between recCA and some provinces

hubei = subset(provinces, NAME=="Hubei")@polygons[[1]]@Polygons[[1]]@coords
guangdong = subset(provinces, NAME=="Guangdong")@polygons[[1]]@Polygons[[1]]@coords
for (i in 1:length(clades))
	{
		if (clades[i] == "SC1") { province = guangdong; provinceName = "Guangdong" }
		if (clades[i] == "SC2") { province = hubei; provinceName = "Hubei" }
		positions = read.table(paste0(analysis,"/",clades[i],"/noHumPanLoc/parentL.txt"), head=T)[,1:2]
		minDistances = matrix(nrow=dim(positions)[1], ncol=1)
		for (j in 1:dim(positions)[1])
			{
				x1 = cbind(positions[j,"longitude"], positions[j,"latitude"])
				x2 = cbind(province[,1], province[,2])
				minDistances[j,1] = min(rdist.earth(x1, x2, miles=F, R=NULL))
			}
		if (writingFiles)
			{
				write.table(minDistances, paste0(clades[i],"_minDistance_recCA_",provinceName), row.names=F, col.names=F, quote=F)
			}
		medianV = round(median(minDistances),2); HPD = round(HDInterval::hdi(minDistances)[1:2],2)
		cat("\t",clades[i],": median of the min. distances between recCA positions and the ",provinceName," = ",medianV," km, 95% HPD = [",HPD[1],"-",HPD[2],"]\n",sep="")
	}

# 9. Visualising all continuous phylogeographic reconstructions

logTransformation1 = FALSE; logTransformation2 = TRUE
directory2 = "humPangLocSI"; directory2 = "noHumPanLoc"
for (i in 1:length(clades))
	{
		if (i == 1) { hosts = c(); mccs = c() }
		mcc = read.csv(paste0(analysis,"/",clades[i],"/humPangLocSI/transPoW/",clades[i],"_",selected_NRRs[i],"_1.csv"), head=T)
		hosts = c(hosts, unique(mcc[,"host"])); mccs = rbind(mccs, mcc)
	}
hosts = unique(hosts); hosts = hosts[order(hosts)]; hosts = hosts[!is.na(hosts)]; counts = rep(NA, length(hosts))
for (i in 1:length(hosts)) counts[i] = sum(mccs[,"host"]==hosts[i], na.rm=T)
hosts = hosts[order(counts, hosts, decreasing=T)]
croppingPolygons = FALSE; polygons_list2 = list()
for (i in 1:length(clades))
	{
		polygons_list1 = list(); mostRecentSamplingDatum = mostRecentSamplingDates[i]
		files = list.files(paste0(analysis,"/",clades[i],"/",directory2,"/transPoW"))
		files = gsub("_2.csv","",files[which(grepl("_2.csv",files))])
		files = files[order(as.numeric(gsub("SC1_NRR","",gsub("SC2_NRR","",files))))]
		for (j in 1:length(files))
			{
				localTreesDirectory = paste0(analysis,"/",clades[i],"/",directory2,"/transPoW/",files[j],"_ext2")
				nberOfExtractionFiles = length(list.files(localTreesDirectory)); prob = 0.80; precision = 50;
				# tree = readAnnotatedNexus(paste0(analysis,"/",clades[i],"/",directory2,"/transPoW/",files[j],".tree"))
				startDatum = mostRecentSamplingDatum-1000 # mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]]
				polygons_list1[[j]] = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
			}
		polygons_list2[[i]] = polygons_list1
	}
saveRDS(polygons_list2, "Analyses_11012023.rds"); polygons_list = polygons_list2
polygons_list = readRDS("Analyses_11012023.rds"); croppingPolygons = FALSE; reportingHosts = FALSE
for (g in 1:2) { for (h in 1:length(clades)) {		
		pdf(paste0(clades[h],"_all_NRRs_fig",g,".pdf"), width=8.0, height=9.0) # dev.new(width=8.0, height=9.0)
		par(mfrow=c(5,3), oma=c(0.0,0.7,0.0,0.0), mar=c(0.0,0.0,0.0,0.0), lwd=0.2, col="gray30")
		if ((h==1)&(g==1)) { iS = c(1:13,15:16) }
		if ((h==1)&(g==2)) { iS = 17:31 }
		if ((h==2)&(g==1)) { iS = c(1:2,4:16) }
		if ((h==2)&(g==2)) { iS = 17:27 }
		for (i in iS) # for SC2 (part 2)
			{
				tree = readAnnotatedNexus(paste0(analysis,"/",clades[h],"/",directory2,"/transPoW/",clades[h],"_NRR",i,".tree"))
				mostRecentSamplingDatum = mostRecentSamplingDates[h]; correctedBranches = c()
				startDatum = mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]]
				for (j in 1:length(tree$edge.length))
					{
						if (tree$edge.length[j] < 0)
							{
								buffer = tree$edge.length[j]; tree$edge.length[j] = 0
								branchToCorrect = which(tree$edge[,2]==tree$edge[j,1])
								if (!branchToCorrect%in%correctedBranches)
									{
										tree$edge.length[branchToCorrect] = tree$edge.length[branchToCorrect]+buffer
										correctedBranches = c(correctedBranches, branchToCorrect)
									}
							}
					}
				buffer = tree
				if (logTransformation1 == TRUE)
					{
						if (tree$root.annotation$`height_95%_HPD`[[1]] > 100)
							{
								buffer$root.annotation$`height_95%_HPD`[[1]] = 100+log(tree$root.annotation$`height_95%_HPD`[[1]]-100)
							}
						if (tree$root.annotation$`height_95%_HPD`[[2]] > 100)
							{
								buffer$root.annotation$`height_95%_HPD`[[2]] = 100+log(tree$root.annotation$`height_95%_HPD`[[2]]-100)
							}
					}
				if (logTransformation2 == TRUE)
					{
						buffer$root.annotation$`height_95%_HPD`[[1]] = log(tree$root.annotation$`height_95%_HPD`[[1]]+1)
						buffer$root.annotation$`height_95%_HPD`[[2]] = log(tree$root.annotation$`height_95%_HPD`[[2]]+1)
					}
				for (j in 1:dim(tree$edge)[1])
					{
						nodeAge1 = max(nodeHeights(tree))-nodeHeights(tree)[j,1]
						nodeAge2 = max(nodeHeights(tree))-nodeHeights(tree)[j,2]
						if (logTransformation2 == TRUE)
							{
								nodeAge1 = max(nodeHeights(tree))+1-nodeHeights(tree)[j,1]
								nodeAge2 = max(nodeHeights(tree))+1-nodeHeights(tree)[j,2]
							}
						if (logTransformation1 == TRUE)
							{
								if (nodeAge1 > 100) nodeAge1 = 100+log(nodeAge1-100)
								if (nodeAge2 > 100) nodeAge2 = 100+log(nodeAge2-100)
							}
						if (logTransformation2 == TRUE)
							{
								nodeAge1 = log(nodeAge1); nodeAge2 = log(nodeAge2)
							}
						buffer$edge.length[j] = nodeAge1-nodeAge2
					}
				tree = buffer; rootHeight = max(nodeHeights(tree))
				root_time = mostRecentSamplingDatum-rootHeight; tree$tip.label = gsub("'","",tree$tip.label)
				minYear = mostRecentSamplingDatum-tree$root.annotation$`height_95%_HPD`[[2]]; maxYear = mostRecentSamplingDatum
				mcc = read.csv(paste0(analysis,"/",clades[h],"/",directory2,"/transPoW/",clades[h],"_NRR",i,"_2.csv"), head=T)
				endYears_colours = rep(NA, dim(mcc)[1])
				for (j in 1:length(endYears_colours))
					{
						if (logTransformation1 == TRUE) nodeAge = maxYear-mcc[j,"endYear"]
						if (logTransformation2 == TRUE) nodeAge = maxYear+1-mcc[j,"endYear"]
						if ((logTransformation1 == TRUE)&(nodeAge > 100)) nodeAge = 100+log(nodeAge-100)
						if (logTransformation2 == TRUE) nodeAge = log(nodeAge)
						endYearMod = maxYear-nodeAge
						endYear_index = (((endYearMod-minYear)/(maxYear-minYear))*100)+1
						endYears_colours[j] = cols2[endYear_index]
					}
				polygons = polygons_list[[h]][[i]]; polygons_colours = rep(NA, length(polygons))
				for (j in 1:length(polygons))
					{
						date = as.numeric(names(polygons[[j]])); polyAge = NULL
						if (logTransformation1 == TRUE) polyAge = maxYear-date
						if (logTransformation2 == TRUE) polyAge = maxYear+1-date
						if ((logTransformation1 == TRUE)&(nodeAge > 100)) polyAge = 100+log(polyAge-100)
						if (logTransformation2 == TRUE) polyAge = log(polyAge)
						endYearMod = maxYear-polyAge
						endYear_index = (((endYearMod-minYear)/(maxYear-minYear))*100)+1
						polygons_colours[j] = paste0(cols2[endYear_index],"20")
					}
				plot(background, col=cols1, axes=F, ann=F, box=F, legend=F)
				lines(borders, col="white", lwd=0.3)
				mtext(paste0(clades[h]," - NRR",i), side=3, line=-2.0, at=103, col="gray30", cex=0.5)
				for (j in 1:length(polygons))
					{
						for (k in 1:length(polygons[[j]]@polygons))
							{
								polygons[[j]]@polygons[[k]] = checkPolygonsHoles(polygons[[j]]@polygons[[k]])
							}
						pol = polygons[[j]]; crs(pol) = crs(background)
						if (croppingPolygons == TRUE) pol = crop(pol, provinces)
						plot(pol, axes=F, col=polygons_colours[j], add=T, border=NA)
					}
				for (j in 1:dim(mcc)[1])
					{
						curvedarrow(cbind(mcc[j,"startLon"],mcc[j,"startLat"]), cbind(mcc[j,"endLon"],mcc[j,"endLat"]), arr.length=0,
						  		  	arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
					}
				for (j in dim(mcc)[1]:1)
					{
						if (!mcc[j,"node1"]%in%mcc[,"node2"])
							{
								startYears_index = (((mcc[j,"startYear"]-minYear)/(maxYear-minYear))*100)+1
								points(mcc[j,"startLon"], mcc[j,"startLat"], pch=16, col=cols2[startYears_index], cex=0.4)
								points(mcc[j,"startLon"], mcc[j,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
							}
						if (mcc[j,"node2"]%in%mcc[,"node1"])
							{	
								points(mcc[j,"endLon"], mcc[j,"endLat"], pch=16, col=endYears_colours[j], cex=0.4)
								points(mcc[j,"endLon"], mcc[j,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
							}	else		{
								points(mcc[j,"endLon"], mcc[j,"endLat"], pch=23, bg=endYears_colours[j], cex=0.35, col=NA)
								points(mcc[j,"endLon"], mcc[j,"endLat"], pch=5, col="gray30", lwd=0.2, cex=0.35)
							}
					}
				selectedDates = c(-10000, -3000, 0, 1000, 1500, 1900, 2000, 2020); selectedLabels = selectedDates
				for (j in 1:length(selectedDates))
					{
						if (logTransformation1 == TRUE)
							{
								if ((mostRecentSamplingDatum-selectedDates[j]) > 100)
									{
										selectedDates[j] = mostRecentSamplingDatum-100-log((mostRecentSamplingDatum-100)-selectedDates[j])
									}
							}
						if (logTransformation2 == TRUE)
							{
								selectedDates[j] = mostRecentSamplingDatum-log(mostRecentSamplingDatum+1-selectedDates[j])
							}
					}
				rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minYear; rast[2] = maxYear
				plot(rast, legend.only=T, add=T, col=cols2, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.840,0.850,0.090,0.905), # smallplot=c(0.853,0.859,0.090,0.905)
			 		 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=F,
					 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, tck=-1.2, col.tick="gray30", col="gray30", col.axis="gray30",
					 				line=0, mgp=c(0,0.45,0), at=selectedDates, labels=selectedLabels))
			}
		dev.off()
	}}

# 10. Mapping the bat species richness map on the topographic background

background_world = crop(raster("All_Natural_Earth_files/World_background.tif"), extent(-170,180,-63,90))
background_world[background_world[]==106] = NA; r = background_world
cols4 = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]
pdf(paste0("SC1-2_NRR_panel8.pdf"), width=7, height=3.5)
par(oma=c(0,0,0,0), mar=c(3,3,3,0))
plot(background_world, col=cols4, box=F, axes=F, legend=F)
plot(study_area, col="gray60", lwd=1.0, border=NA, add=T)
dev.off()
pdf(paste0("SC1-2_NRR_panel9.pdf"), width=7, height=4.8)
par(oma=c(0,0,0,0), mar=c(3,3,3,0), lwd=0.2, col="gray30")
srm = raster("Bat_species_richness/Bat_species_richness.tif")
srm[srm[]==0] = NA; cols5 = colorRampPalette(brewer.pal(9,"YlGn"))(121)[20:121]
plot(background, col=cols1, axes=F, ann=F, box=F, legend=F)
plot(srm, col=cols5, add=T, axes=F, ann=F, box=F, legend=F)
plot(srm, legend.only=T, add=T, col=cols5, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.700,0.708,0.15,0.45),
	 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=F,
	 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, tck=-1.2, col.tick="gray30", col="gray30", col.axis="gray30", line=0, mgp=c(0,0.45,0)))						
dev.off()

