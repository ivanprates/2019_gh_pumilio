###################################
## Script by Ivan Prates, July 2018
## ivanprates.org

## This script does two things:
## PART 1. Extracts a presence/absence matrix for ant species in each site sampled for Oophaga pumilio alkaloids;
## PART 2. Estimates a a Sorensen dissimilarity matrix for the ant data.
## Ant presence is based on species distribution models for each ant species, developed in Maxent.

# Getting ready
#install.packages("fossil")
#install.packages("raster")
library(fossil)
library(raster)

# Selecting path depending on which computer I'm using
#path = "~/Dropbox (Smithsonian)/2018_pumilio/2018-07/"
path = "C:/Users/RhinellaX/Dropbox/Science/MYPAPERS_ongoing/2018_pumilio/2018-07/"
path = "~/Dropbox/science/MYPAPERS_ongoing/2018_pumilio/2018-07/"


# Making folder to save the files we'll generate
dir.create(paste0(path, "DATA_dissimilarity_matrices"))

## PART 1. Extracting a presence/absence matrix for ant species in each site sampled for Oophaga pumilio alkaloids

# Importing Oophaga pumilio site information
sites_extract = read.csv(path, "DATA_alkaloids_saporito2007/2018-07_sites_46pixels_1km-pixels_corrected_for_land.csv", header = TRUE, row.names = 1)
coord_extract = sites_extract[, 2:3]
coord_extract

# Defining thershold for conversion of continuous species distribution model to binary surfaces
threshold = "T10" # 10th percentile presence threshold
#threshold = "MPT" # Minimum training presence

# List ant species distributions models
ant_models = list.files(path = paste0("~/Dropbox (Smithsonian)/Models/models_", threshold, "_correct"), pattern = '*.asc', full.names = TRUE)

# Stack layers from extraction
ant_layers = raster::stack(ant_models) # Careful with other files in the same folder with "asc.aux.xml" format!

# Extract raster to points to get ant presence/absence data
ant_presabs = raster::extract(ant_layers, coord_extract) # important to specify package as there are many functions 'extract' out there!
ant_presabs[is.na(ant_presabs)] = 0 # replacing NA for 0 (abscences)
ant_presabs

# Combine ant data to site data
ant_data = cbind(sites_extract, ant_presabs)

# Save the presence/absence ant matrix
write.csv(ant_data, file = paste0(path, "DATA_alkaloid_bearing_ants/2018-07_presabs_ants_", threshold, ".csv"))

## PART 2. Estimating a a Sorensen dissimilarity matrix for the ant data

# To import ant matrix again:
ant_import = read.csv(file = paste0(path, "DATA_alkaloid_bearing_ants/2018-07_presabs_ants_", threshold, "_with_tetramorium.csv"), header = TRUE)
ant_presabs = ant_import[, 7:length(ant_import)]

# Transposing as required by the fossil package for estimation of a distance matrix
ant_presabs = t(ant_presabs)
#colnames(ant_presabs) = sites_extract$alk_siteID # Renaming columns as alk_siteID
colnames(ant_presabs) = ant_import$alk_siteID # Renaming columns as alk_siteID
ant_presabs

# Estimating a sorensen distance matrix
pair_Sorensen = ecol.dist(ant_presabs, method = sorenson, type = "dis")
pair_Sorensen = as.matrix(pair_Sorensen)
pair_Sorensen
#write.csv(pair_Sorensen, file = paste0(path, "DATA_dissimilarity_matrices/2018-07_Sorensen_ants_", threshold, ".csv"))
write.csv(pair_Sorensen, file = paste0(path, "DATA_dissimilarity_matrices/2018-07_Sorensen_ants_", threshold, "_with_tetramorium.csv"))
