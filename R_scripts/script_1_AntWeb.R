############################################
## R script by Ivan Prates (ivanprates.org), Smithsonian National Museum of Natural History, Washington, DC, USA, June 2019.
## This script does the following:
## Retrieves from the AntWeb database which ant genera occur sympatrically with Ooophaga pumilio.
## Dowloads records for each ant species in these target ant genera.
## This script was used to obtain data from the AntWeb database in June 2017. Functionality not tested ever since.

## Part 1: Getting ready:

## Installing packages:
#install.packages("AntWeb", dependencies = TRUE)
#install.packages("devtools")
#install.packages("reshape")
#install.packages("data.table")

## Loading packages:
library(AntWeb)
library(devtools)
library(reshape)
library(data.table)

## Setting working directory:
setwd("~/Dropbox/ivan_lab/2017_pumilio/DATA_alkaloid-bearing_ants")

## PART 2: Retrieve from AntWeb which ant genera occur sympatrically with Ooophaga pumilio:

## Reading list of alkaloid-bearing ant genera:
ant_list = read.table("ant_genera.csv", sep = ",", head = TRUE)
ant_list$taxon = as.character(ant_list$taxon)
str(ant_list)

## Creating an empty list to populate with AntWeb data:
aw_list = vector("list", dim(ant_list)[1])
names(aw_list) = ant_list[,1]
names(aw_list) 

## Getting data from AntWeb as a loop:
## For Costa Rica:
for(i in 1: dim(ant_list)[1]){
  temp = as.data.frame(aw_data(country = "Costa Rica", limit = 20000, genus = ant_list$taxon[i])$data)
  aw_list[[i]] = temp
  rm(temp)
  rm(i)
}
results_CR = rbindlist(l = aw_list, fill = T)
results_CR

## For Nicaragua:
for(i in 1: dim(ant_list)[1]){
  temp = as.data.frame(aw_data(country = "Nicaragua", limit = 20000, genus = ant_list$taxon[i])$data)
  aw_list[[i]] = temp
  rm(temp)
  rm(i)
}
results_NI = rbindlist(l = aw_list, fill = T)
results_NI

## For Panama:
for(i in 1: dim(ant_list)[1]){
  temp = as.data.frame(aw_data(country = "Panama", limit = 20000, genus = ant_list$taxon[i])$data)
  aw_list[[i]] = temp
  rm(temp)
  rm(i)
}
results_PA = rbindlist(l = aw_list, fill = T)
results_PA

## Combining the data:
joined_list = list(results_CR, results_NI, results_PA)
ant_records = rbindlist(l = joined_list, fill = T)
write.csv(ant_records, file = "ant_species.csv")

## PART 3: Getting records for each ant species in target ant genera:

#ant_records = read.table("ant_species.csv", sep = ",", head = TRUE)

## Changing your list from factor (default) to character (required by AntWeb::aw_data): 
ant_records$scientific_name = as.character(ant_records$scientific_name)

## Checking if all looks good:
ant_species = unique(ant_records$scientific_name)

## Creating a list (object) that will receive the downloaded data:
aw_list = vector("list", length(ant_species)[1]) # dim = dimension of the object. The "[1]" part corresponds the number of rows, that is, taxa.

## Naming the elements in this list - each corresponds to one target taxon:
names(aw_list) = ant_species

## Now getting data from AntWeb as a loop:
for(i in 1: length(ant_species)[1]){ # Defines the sequence of target taxa: from the first (1) to the last, tells R to do what follows the {
  ## We'll create a temporary object to store the data coming in every loop. For each species/genera, it will save $data, corresponding to the data table provided by AntWeb.
  temp = as.data.frame(aw_data(georeferenced = TRUE, limit = 20000, bbox = '40,-140,-40,-30', scientific_name = ant_species[i])$data) # Here, change settings as needed (see function aw_data documentation).
  aw_list[[i]] = temp # Place temporary data for each taxon in it's corresponding slot in the list we made before (aw_list).
  rm(temp) # Removing the temporary object for the next round (avoid cluttering of the environment space).
  rm(i) # Removing the order index.
} # Done!

## Now merging data from all taxa into a single table, using rbind:
results = rbindlist(l = aw_list, fill = TRUE)
colnames(results)

## Saving results as csv:
write.csv(results, file = "ant_records.csv")
save.image(file = "AntWeb.RData")

## Done!
