############################################
## R script by Ivan Prates (ivanprates.org), Smithsonian National Museum of Natural History, Washington, DC, USA, June 2019.
## This script does the following:
## Estimates a Sorensen dissimilarity matrix for the the complete Oophaga pumilio alkaloid dataset (both for individual alkaloids and alkaloid classes).
## Selects only those alkaloids and alkaloid classes known to occur in ants, and then estimates a Sorensen dissimilarity matrix.
## Estimates a matrix of geographic distances between sites sampled for Oophaga pumilio alkaloids.
## Assigns Oophaga pumilio cytochrome B sequences to sites sampled for alkaloids (as informed by Voronoi polygons estimated in ArcGIS).

## PART 1: Getting ready:

## Installing packages:
#install.packages("fossil")
#install.packages("magrittr")
#install.packages("plyr")
#install.packages("reshape2")

## Loading packages:
library(fossil)
library(magrittr)
library(plyr)
library(reshape2)

## Selecting path depending on which computer I'm using:
path = "~/Dropbox/2018_pumilio/2018-07/"

## Creating folders to save the files we'll generate:
dir.create(paste0(path, "DATA_alkaloids_in_ants"))
dir.create(paste0(path, "DATA_dissimilarity_matrices"))

## Selecting alkaloid data type (i.e., individual alkaloids versus alkaloid classes):
data_type = "alkaloid"
#data_type = "alkaloid_class"

## PART 2: Estimating a Sorensen dissimilarity matrix for the alkaloid data:

## Importing presence/absence matrix of alkaloid data:
input_data = read.csv(file = paste0(path, "DATA_alkaloids_saporito2007/2018-07_presabs_", data_type, ".csv"), header = TRUE, row.names = 1) # read presence/abscence matrix
alk_data = input_data[, 6:dim(input_data)[2]] # Selecting columns corresponding to alkaloid data
rownames(alk_data) = input_data$alk_siteID

## Transposing thius matrix as required by the fossil package for estimation of a distance (Sorensen) matrix:
alk_data = t(alk_data)
alk_data

## Estimating a Sorensen distance matrix for the alkaloid data:
pair_sorensen = ecol.dist(alk_data, method = sorenson, type = "dis")
pair_sorensen = as.matrix(pair_sorensen)
pair_sorensen
write.csv(pair_sorensen, file = paste0(path, "DATA_dissimilarity_matrices/2018-07_Sorensen_", data_type, ".csv")) # Writing this matrix as a csv file.

## PART 3: Selecting only those alkaloids from alkaloid classes known to occur in ants:

## Converting a presence/absence matrix to a list: 
input_data = read.csv(file = paste0(path, "DATA_alkaloids_saporito2007/2018-07_presabs_", data_type, ".csv"), header = TRUE, row.names = 1) # reads presence/abscence matrix
molten  = melt(input_data, id.vars = c("saporito_siteID", "latitude",	"longitude", "alk_siteID", "Nskins"))
molten$value
colnames(molten) # lists column names
only_presence = molten[molten$value == 1,]
only_presence = only_presence[, -7]
only_presence$variable %<>% gsub(pattern = "X", replacement = "") # Eliminating the "X" that R adds to headers when they're numbers.
write.csv(only_presence, file = paste0(path, "DATA_alkaloids_saporito2007/2018-07_list_", data_type, ".csv"), row.names = FALSE) # Saving.

## Importing list of alkaloids known to occur in ants (from literature review) to act as criteria for alkaloid selection:
criteria_list = read.csv(file = paste0(path, "DATA_alkaloid_bearing_ants/2018-07_list_only_alkaloids_found_in_ants_including_tricyclics.csv"), header = TRUE)

## Creating a list that will serve as a 'mask' to select only those alkaloids under criteria:
criteria = criteria_list[, colnames(criteria_list) == data_type]
length(unique(criteria))

## Using this mask, extract only those alkaloids (or alkaloid classes) under criteria:
selected = only_presence[only_presence$variable %in% criteria, ]
length(selected$variable)

## Saving final alkaloid list:
write.csv(selected, file = paste0(path, "DATA_alkaloids_in_ants/2018-07_list_", data_type, "_in_ants_including_tricyclics.csv"))

## Converting from list to presence/absence matrix:
matrix = melt(selected, id.var = "variable", measure.var = "alk_siteID")
presabs = with(matrix, table(variable, value))
presabs = t(presabs)
write.csv(presabs, file = paste0(path, "DATA_alkaloids_in_ants/2018-07_presabs_", data_type, "_in_ants_including_tricyclics.csv"))

## Transposing matrix as required by the fossil package for estimation of a distance matrix:
presabs = t(presabs)
presabs

## Estimating a Sorensen distance matrix:
pair_Sorensen = ecol.dist(presabs, method = sorenson, type = "dis")
pair_Sorensen = as.matrix(pair_Sorensen)
pair_Sorensen
write.csv(pair_Sorensen, file = paste0(path, "DATA_dissimilarity_matrices/2018-07_Sorensen_", data_type, "_in_ants_including_tricyclics.csv")) # Saving.

## PART 4: Estimating a matrix of geographic distances between sites sampled for Oophaga pumilio alkaloids:

## Importing site data:
sites = read.csv(file = paste0(path, "DATA_alkaloids_saporito2007/2018-07_sites_46pixels_1km-pixels_corrected_for_land.csv"), header = TRUE,  row.names = 1)
lat_long = sites[, 2:3] # Column order matters: Has to be longitude, latitude
lat_long

## Estimating geo dist matrix:
geo = earth.dist(lat_long, dist = FALSE) # geo dist matrix
colnames(geo) = input_data$alk_siteID
rownames(geo) = input_data$alk_siteID
geo
write.csv(geo, file = paste0(path, "DATA_dissimilarity_matrices/2018-07_geographic_distance_matrix.csv")) # Saving.

## PART 5. Assigning cytochrome B sequences to sites sampled for alkaloids:

## Importing genetic data and locality data for matching:
gen_data = read.csv(file = paste0(path, "DATA_genetics/Hauswaldt_et_al_2011_data/pumilio_cytB_Hauswaldt_et_al_2011.csv"), header = TRUE) # Reads genetic data.
sites = read.csv(file = paste0(path, "DATA_alkaloids_saporito2007/2018-07_sites_46pixels_1km-pixels_corrected_for_land.csv"), header = TRUE,  row.names = 1)

## Creating an empty list to populate with sequence data:
seq_list = vector("list", length(sites$gen_siteID)) # It will have the same length as the number of alkaloid sites (n = 46).
names(seq_list) = sites$alk_siteID # Changing names.
names(seq_list)

## Loop: Matching DNA sequences with sites sampled for alkaloids (this was informed by Voronoi polygons in ArcGIS):
for(i in 1:nrow(sites)){
  gen_siteID = sites$gen_siteID[i] # Picking the name of the site sampled for genetics.
  seq_list[[i]] = gen_data$cytb_seq[gen_data$gen_siteID == gen_siteID] # Selecting sequences that match sites sampled for alkaloids.
  } # Done
pre_alignment = ldply(seq_list, data.frame) # Converting list into a data frame.
colnames(pre_alignment) = c("alk_siteID", "cytb_seq") # Changing column names.

## Now we'll name individual sequences matching the corresponding sites sampled for alkaloids.
add_ids = vector("list", length(sites$gen_siteID)) # creating empty vector to receive individual number.

## Now preparing alignments:

## Option 1: Fasta format:

## Loop: Changing names of sequences to add individual number and ">" in preparation for fasta file:
for(i in 1:nrow(sites)){
  temp_slot = paste0(">", unique(pre_alignment$alk_siteID)[i], "_", 1:table(pre_alignment$alk_siteID)[i]) # Creating numbers for each individual sampled for DNA.
  add_ids[[i]] = temp_slot # We're using a temporary object to store the names.
  rm(temp_slot) # Removing temporary object.
  } # Done.

add_ids = ldply(add_ids, data.frame) # Converting list into a data frame.
pre_alignment$alk_siteID = add_ids[, 1] # Replacing alkaloid site ID to include also individual sequence names.
colnames(pre_alignment) = c("sample", "cytb_seq") # Changing column names.
pre_alignment = pre_alignment[, 1:2] # Using first (individual name) and second (DNA sequence) columns to make alignment.

# Saving as a final alignment:
write.table(pre_alignment, file = paste0(path, "DATA_genetics/pumilio_cytb_matching_alkaloid_sites.fasta"), row.names = FALSE, quote = FALSE, sep = "\n", eol = " \n", col.names = FALSE)

## Option 2: Mega format:

## Loop: Changing names of sequences to add individual number and ">" in preparation for fasta file:
for(i in 1:nrow(sites)){
  temp_slot = paste0("#", unique(pre_alignment$alk_siteID)[i], "_", 1:table(pre_alignment$alk_siteID)[i], "{", unique(pre_alignment$alk_siteID)[i], "}") # Creating numbers for each individual sampled for DNA.
  add_ids[[i]] = temp_slot # We're using a temporary object to store the names.
  rm(temp_slot) # Removing temporary object.
} # Done.

add_ids = ldply(add_ids, data.frame) # Converting list into a data frame.
pre_alignment$alk_siteID = add_ids[, 1] # Replacing alkaloid site ID to include also individual sequence names.
colnames(pre_alignment) = c("sample", "cytb_seq") # Changing column names.
pre_alignment = pre_alignment[, 1:2] # Using first (individual name) and second (DNA sequence) columns to make alignment.

## Saving as a final alignment:
write.table(pre_alignment, file = paste0(path, "DATA_genetics/pumilio_cytb_matching_alkaloid_sites.meg"), row.names = FALSE, quote = FALSE, sep = "\n", eol = " \n", col.names = FALSE)

## Add Mega header manually to the output file, as follows:
## "#Mega"
## "!Title pumilio_cytb_matching_alkaloid_sites;"

## Making the genetic distance matrix generated in MEGA symetrical for use in MMRR:

## Reading half-matrix output by MEGA:
half_matrix = read.csv(file = "C:/Users/RhinellaX/Dropbox/Science/MYPAPERS_ongoing/2018_pumilio/2018-07/DATA_genetics/pumilio_cytb_dist_p-distances_half-matrix.csv", header = TRUE, row.names = 1)
dim(half_matrix)

## Creating a square matrix with the dimensions of the final (symetrical) matrix:
full_matrix = matrix(data = 0, nrow = nrow(half_matrix), ncol = ncol(half_matrix))
dim(full_matrix)

## Making new matrix extract upper half from half matrix:
full_matrix[lower.tri(full_matrix)] = half_matrix[lower.tri(half_matrix)]

full_matrix == half_matrix # Checking whether diagonals match.
full_matrix[upper.tri(half_matrix)] = t(half_matrix)[upper.tri(t(half_matrix))]
isSymmetric(full_matrix) # Check: Did it work (i.e., is matrix now symetrical)?

## Adding column and row names:
row.names(full_matrix) = row.names(half_matrix)
colnames(full_matrix) = row.names(half_matrix)

# Saving final matrix:
write.csv(full_matrix, file = "C:/Users/RhinellaX/Dropbox/Science/MYPAPERS_ongoing/2018_pumilio/2018-07/DATA_genetics/2018-07_pumilio_cytb_dist_p-distances.csv", row.names = TRUE)

# Done!
