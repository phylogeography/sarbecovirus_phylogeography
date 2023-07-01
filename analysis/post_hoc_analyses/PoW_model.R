###Libraries required for the analyses
library("ape")
library("stats")
library("nleqslv")
library("scales")
library("treeio")
library("ggtree")
require(tidytree)
##########################################

#######
#path to log file for the inferred short-term substitution rate (*.log output file from BEAST v1.10.4)
#this should be the inferred rate for either the SC1 and SC2 datasets
filepath_rates <- "/Path/to/directory/SC1&2_POW/SC2/rate/SC2.log"

#goes from 1 to 31 for SC1 and 27 for SC2, corresponding to the number of non-recombinant segments for each dataset
for (file_num in 1:27){
  print(paste("file_number:",file_num))
  #path to ultrametric distance trees (*.trees output file from BEAST v1.10.4)
  filepath_distances <- paste("/Path/to/directory/SC1&2_POW/SC2/nohumPan/SC2_NRR", as.character(file_num), ".trees", sep="")

  #######
  
  #######
  #con reads each row in the log file, removing the first 10% of rows as burn-in
  con = readLines(paste(filepath_rates, sep = ""))
  con <- con[c(5:length(con))]
  con <- con[c((round(length(con)/10)+1):length(con))]
  #######
  
  #sample_size sets the number of samples taken at random from the post-burn-in posterior rate distribution and distance trees.
  sample_size=100
  
  #######
  #selecting states from the log file at random
  sampled_states <- sample(con,sample_size)
  #######
  
  #######
  #parameters of the substitution model
  
  #uncomment below lines for SC1
  #pA=0.25
  #pC=0.226
  #pG=0.169
  #pT=0.355
  #kap=5.89
  
  #uncomment below for SC2
  pA=0.28
  pC=0.208
  pG=0.178
  pT=0.335
  kap=4.254
  
  
  #######
  #reading the mean.Rate column for each segment 
  #uncomment below for SC1
  #meanrate_column=154+file_num
  #uncomment below for SC2
  meanrate_column=142+file_num
  
  rates <- c()
  for(i in 1:length(sampled_states)){
    rates <- c(rates,as.numeric(strsplit(sampled_states[[i]],"\t")[[1]][[meanrate_column]]))
  }
  #######
  
  #use the best fit value of muMax in RNA viruses for PoW transformation of every tree
  muMax = 3.65*10**(-2)
  muMaxs=rep(muMax,sample_size)
  
  #uncomment below if you like to add variation around the value of muMax (not recommended without prior info)
  #muMax = 3.65*10**(-2)
  #muMax_stdev=10.0*10**-3
  #muMaxs=rnorm(sample_size,muMax,muMax_stdev)
  
  #con2 reads the distance tree file
  con2 = readLines(paste(filepath_distances, sep = ""))
  
  for(i in 1:length(con2)){
    if(paste(strsplit(con2[[i]],'')[[1]][c(1:4)],collapse = '')=='tree'){
      break
    }
  }
  
  #tree_list includes a list of produced trees from BEAST, removing the first 10% as burn-in
  tree_list <- con2[c(i:(length(con2)-1))] 
  tree_list <- tree_list[c((round(length(tree_list)/10)+1):length(tree_list))]
  sampled_trees <- sample(tree_list,sample_size)
  
  
  #create a new tree file based on the subsampling in the previous step 
  recreated_treefile <- c(con2[c(1:(i-1))],sampled_trees,'End;')
  write.table(recreated_treefile, '/Path/to/directory/sampled_trees.trees', sep = "\n", quote = FALSE, row.names = FALSE)
  
  #import the tree file created in the previous step
  filename <- "/Path/to/directory/sampled_trees.trees"
  importedtree <- read.nexus(filename, tree.names = NULL, force.multi = FALSE)
  
  #transform each of the sampled distance trees using the PoW method 
  for(el in 1:sample_size){
    print(el)
    muMax = muMaxs[[el]]
    meanRate = rates[[el]]
    
    pY = pT + pC
    pR = pA + pG
    delta = 0.1
    betaMax = muMax/(2*((pT*pC + pA*pG)*kap + pY*pR))
    steps = seq(0,log10(betaMax)+9,delta)
    steps2 = seq(1,length(steps),1)
    lambda <- function(l){
      final <- sum(exp(l*(1/delta*steps+1))*(2*((pT*pC + pA*pG)*kap + pY*pR))*10**(-9+steps))/sum(exp(l*steps2)) - meanRate
      return (final)
    }
    
    LAMBDA = nleqslv(0,lambda)$x
    proportions = exp(LAMBDA*(1/delta*steps+1))/sum(exp(LAMBDA*steps2))
    ratespread = 10**(-9+steps)
    
    e2 <- function(t) exp(-ratespread*t)
    e3 <- function(t) exp(-(pR*kap+pY)*ratespread*t)
    e4 <- function(t) exp(-(pY*kap+pR)*ratespread*t)
    pTC <- function(t) pC+(pC*pR)/pY*e2(t)-pC/pY*e4(t)
    pAG <- function(t) pG+(pG*pY)/pR*e2(t)-pG/pR*e3(t)
    pTA <- function(t) pA*(1-e2(t))
    pTG <- function(t) pG*(1-e2(t))
    pCA <- function(t) pA*(1-e2(t))
    pCG <- function(t) pG*(1-e2(t))
    s1 <- function(t) sum(proportions*2*pT*pTC(t))
    s2 <- function(t) sum(proportions*2*pA*pAG(t))
    v <- function(t) sum(proportions*(2*pT*pTA(t)+2*pT*pTG(t)+2*pC*pCA(t)+2*pC*pCG(t)))
    a1 <- function(t) -log(1-(pY*s1(t))/(2*pT*pC)-v(t)/(2*pY))
    a2 <- function(t) -log(1-(pR*s2(t))/(2*pA*pG)-v(t)/(2*pR))
    b <- function(t) -log(1-v(t)/(2*pY*pR))
    k1 <- function(t) (a1(t)-pR*b(t))/(pY*b(t))
    k2 <- function(t) (a2(t)-pY*b(t))/(pR*b(t))
    func <- function(t){
      final <- 2*pT*pC/pY*(a1(t)-pR*b(t))+(2*pA*pG/pR)*(a2(t)-pY*b(t))+2*pY*pR*b(t) - totalBranchLength
      return(final)
    }
    
    
    allBranches=mapply(c,importedtree[[el]]$edge[,1],importedtree[[el]]$edge[,2],(importedtree[[el]]$edge.length),SIMPLIFY = FALSE)
    
    Allpaths <- list(list())
    for(j in 0:importedtree[[el]]$Nnode+1){
      branches <- list()
      for(i in 1:(length(nodepath(importedtree[[el]],j,importedtree[[el]]$Nnode+2))-1)){
        branches[[i]] <- c(nodepath(importedtree[[el]],j,importedtree[[el]]$Nnode+2)[i+1],nodepath(importedtree[[el]],j,importedtree[[el]]$Nnode+2)[i])
      }
      Allpaths[[j]] <- c(branches)
    }
    
    trackBranches <- list()
    globcount=0
    for(k in 1:length(Allpaths)){
      totalBranchLength=0
      count=0
      convertedHeight<-list()
      for(i in 1:length(Allpaths[[k]])){
        for(j in 1:length(allBranches)){
          if(sum(allBranches[[j]] %in% Allpaths[[k]][[i]],na.rm = TRUE)==2){
            matches=0
            count = count + 1
            globcount = globcount + 1
            trackBranches[[globcount]] <- c(allBranches[[j]][c(1:2)])
            for(m in 1:length(trackBranches)){
              if(sum(trackBranches[[m]] %in% Allpaths[[k]][[i]],na.rm = TRUE)==2){
                matches = matches + 1
              }
            }
            if(count==1 && matches==1){
              totalBranchLength = totalBranchLength + allBranches[[j]][[3]]
              convertedHeight[[count]] <- c(nleqslv(0,func)$x)
              importedtree[[el]]$edge.length[[j]]=convertedHeight[[count]]
            }
            if(count>1 && matches==1){
              totalBranchLength = totalBranchLength + allBranches[[j]][[3]]
              convertedHeight[[count]] <- c(nleqslv(0,func)$x)
              importedtree[[el]]$edge.length[[j]]=convertedHeight[[count]]-convertedHeight[[count-1]]
            }
            if(count>1 && matches>1){
              TRUE
            }
          }
        }
      }
    }
  }
  #save the PoW-transformed trees for each non-recombinant segment in a local directory
  #this file can then be directly imported to TreeAnnotator to construct a PoW maximum clade credibility tree
  ape::write.nexus(importedtree, file=paste("/Path/to/directory/SC1&2_POW/SC2/nohumPan/PoWTransformed_SC2_NRR", as.character(file_num), ".trees", sep=""))
}

