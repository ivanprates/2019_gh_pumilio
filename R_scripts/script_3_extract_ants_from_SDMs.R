############################################
## R script by Ivan Prates (ivanprates.org), Smithsonian National Museum of Natural History, Washington, DC, USA, June 2019.
## This script does the following:
## Extracts a presence/absence matrix for ant species in each site sampled for Oophaga pumilio alkaloids.
## Estimates a a Sorensen dissimilarity matrix for the ant data.
## Ant presence is based on species distribution models for each ant species, developed in Maxent.

## Part 1: Getting ready:

## Installing packages:
#install.packages("fossil")
#install.packages("raster")

## Loading packages:
library(fossil)
library(raster)

## Selecting path depending on which computer I'm using:
path = "~/Dropbox (Smithsonian)/2018_pumilio/2018-07/"

## Making folder to save the files we'll generate:
dir.create(paste0(path, "DATA_dissimilarity_matrices"))

## PART 2. Extracting a presence/absence matrix for ant species in each site sampled for Oophaga pumilio alkaloids:

## Importing Oophaga pumilio site information:
sites_extract = read.csv(path, "DATA_alkaloids_saporito2007/2018-07_sites_46pixels_1km-pixels_corrected_for_land.csv", header = TRUE, row.names = 1)
coord_extract = sites_extract[, 2:3]
coord_extract

## Defining threshold for conversion of continuous species distribution model into binary surfaces:
threshold = "T10" # 10th percentile presence threshold.
#threshold = "MPT" # Minimum training presence (not used).

## Listing ant species distributions models:
ant_models = list.files(path = paste0("~/Dropbox (Smithsonian)/Models/models_", threshold, "_correct"), pattern = '*.asc', full.names = TRUE)

## Stacking layers from extraction:
ant_layers = raster::stack(ant_models) # Careful with other files in the same folder with "asc.aux.xml" format!

## Extracting raster to points to get ant presence/absence data:
ant_presabs = raster::extract(ant_layers, coord_extract) # Important to specify package as there are many functions 'extract' out there!
ant_presabs[is.na(ant_presabs)] = 0 # Replacing NA for 0 (absences).
ant_presabs

## Combining ant data to site data:
ant_data = cbind(sites_extract, ant_presabs)

## Saving the resulting presence/absence ant matrix:
write.csv(ant_data, file = paste0(path, "DATA_alkaloid_bearing_ants/2018-07_presabs_ants_", threshold, ".csv"))

## PART 3. Estimating a a Sorensen dissimilarity matrix for the ant data:

## Importing ant matrix again (if necessary):
ant_import = read.csv(file = paste0(path, "DATA_alkaloid_bearing_ants/2018-07_presabs_ants_", threshold, ".csv"), header = TRUE)
ant_presabs = ant_import[, 7:length(ant_import)]

## Transposing matrix as required by the fossil package for estimation of a distance matrix:
ant_presabs = t(ant_presabs)
#colnames(ant_presabs) = sites_extract$alk_siteID # Renaming columns as alk_siteID
colnames(ant_presabs) = ant_import$alk_siteID # Renaming columns as alk_siteID
ant_presabs

## Estimating a Sorensen distance matrix:
pair_Sorensen = ecol.dist(ant_presabs, method = sorenson, type = "dis")
pair_Sorensen = as.matrix(pair_Sorensen)
pair_Sorensen
write.csv(pair_Sorensen, file = paste0(path, "DATA_dissimilarity_matrices/2018-07_Sorensen_ants_", threshold, ".csv"))

## Done!