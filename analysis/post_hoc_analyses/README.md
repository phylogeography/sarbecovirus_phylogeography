Code to generate Figures 1 and 2 is available in the jupyter notebooks. 

The rest of the post-hoc files are for the phylogeography. `PoW_model.R` subsamples trees outputted by the XML files in the `PoW-transformed_phylogeography` subdirectories and then transforms them using the PoW model from Ghafari et al. (Current Biology, 2021). 

The main R script for the subsequent phylogeographic analysis is in the file `R_post_hoc_analyses.r`. This script is divided into nine blocks of code and performs the following steps:

1. Preparing the GIS files and the different colour scales
2. Extracting the spatio-temporal information embedded in annotated trees
3. Investigating the profile of long-distance dispersal events
4. Investigating the patterns of isolation-by-distance
5. Estimating some dispersal statistics for each reconstruction (not used)
6. Visualising selected phylogenetic trees for SC1 and SC2
7. Visualising selected continuous phylogeographic reconstructions
8. Displaying the estimated position of the human ancestor
9. Visualising all continuous phylogeographic reconstructions

The R script requires the preliminary installation of the following R packages: "adephylo", "diagram", "HDInterval", "lubridate", "maptools", "rgdal", and "seraphim".

The R scripts "Tree_data_extraction1.r" and "Tree_data_extraction2.r" respectively contain a R function to extract spatio-temporal information embedded in annotated maximum clade credibility (MCC) tree and posterior trees inferred by a continuous phylogeographic analysis.
