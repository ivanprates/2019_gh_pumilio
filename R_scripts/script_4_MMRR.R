############################################
## R script by Ivan Prates (ivanprates.org), Smithsonian National Museum of Natural History, Washington, DC, USA, June 2019.
## This script does the following:
## Performs Multiple Matrix Regression with Randomization (MMRR) analyses based on ant, alkaloid, genetic, and geographic distance data.

## Part 1: Getting ready:

## Installing packages:
#install.packages("car", dependencies = TRUE)
#install.packages("fossil")
#install.packages("tseries") # The read.matrix function requires {tseries} package to be installed and loaded.

## Loading packages:
library(car)
library(fossil)
library(tseries)

## Selecting path depending on which computer I'm using:
path = "~/Dropbox (Smithsonian)/2018_pumilio/2018-07/"

## Making folder to save the files we'll generate:
dir.create(paste0(path, "OUTPUTS_MMRR/"))

## Part 2: Importing and setting up distance matrixes to be used in MMRR analyses:

## Defining threshold for conversion of continuous species distribution model into binary surfaces:
threshold = "T10" # 10th percentile presence threshold.
#threshold = "MPT" # Minimum training presence (not used).

## Importing matrices for ant, alkaloid, alkaloid class, and geographic distances:
gen = as.matrix(read.csv(file = paste0(path, "DATA_genetics/2018-07_pumilio_cytb_dist_p-distances.csv"), header = TRUE,  row.names = 1))
geo = as.matrix(read.csv(file = paste0(path, "DATA_dissimilarity_matrices/2018-07_geographic_distance_matrix.csv"), header = TRUE,  row.names = 1))
ants = as.matrix(read.csv(file = paste0(path, "DATA_dissimilarity_matrices/2018-07_Sorensen_ants_", threshold, ".csv"), header = TRUE,  row.names = 1))
alkaloid = as.matrix(read.csv(file = paste0(path, "DATA_dissimilarity_matrices/2018-07_Sorensen_alkaloid.csv"), header = TRUE, row.names = 1))
alkaloid_in_ants = as.matrix(read.csv(file = paste0(path, "DATA_dissimilarity_matrices/2018-07_Sorensen_alkaloid_in_ants_with_tricyclics.csv"), header = TRUE, row.names = 1))
alkaloid_class = as.matrix(read.csv(file = paste0(path, "DATA_dissimilarity_matrices/2018-07_Sorensen_alkaloid_class.csv"), header = TRUE,  row.names = 1))
alkaloid_class_in_ants = as.matrix(read.csv(file = paste0(path, "DATA_dissimilarity_matrices/2018-07_Sorensen_alkaloid_class_in_ants_with_tricyclics.csv"), header = TRUE,  row.names = 1))

## Setting matrix diagonals as NA:
diag(gen) = NA
diag(geo) = NA
diag(ants) = NA
diag(alkaloid) = NA
diag(alkaloid_in_ants) = NA
diag(alkaloid_class) = NA
diag(alkaloid_class_in_ants) = NA

## Extracting minimum values:
min(ants, na.rm = TRUE)
min(alkaloid, na.rm = TRUE)
min(alkaloid_in_ants, na.rm = TRUE)
min(alkaloid_class, na.rm = TRUE)
min(alkaloid_class_in_ants, na.rm = TRUE)

## Extracting maximum values:
max(ants, na.rm = TRUE)
max(alkaloid, na.rm = TRUE)
max(alkaloid_in_ants, na.rm = TRUE)
max(alkaloid_class, na.rm = TRUE)
max(alkaloid_class_in_ants, na.rm = TRUE)

## Part 3: Selecting predictor variables for MMRR analyses:

## Make a list of the explanatory (X) matrices. Order doesn't matter. Can include more than two matrices, if desired.

## Using ants (Sorensen) only:
predictor = ants
predictor_list = list(ants = ants)
predictor_name = "ants"

## Using geographic distances only:
predictor = geo
predictor_list = list(geo = geo)
predictor_name = "geo"

## Using genetic distances only:
predictor = gen
predictor_list = list(gen = gen)
predictor_name = "gen"

## Using ants and genetic distances jointly:
predictor_list = list(ants = ants, gen = gen)
predictor_name = "ant_gen"

## Part 4: Selecting dependent variables for MMRR analyses:

## Setting individual alkaloids as the dependent variable:
dependent_variable = alkaloid
dependent_variable_name = "alkaloid"

## Setting individual alkaloids as the dependent variable:
dependent_variable = alkaloid_in_ants
dependent_variable_name = "alkaloid_in_ants"

## Setting alkaloid classes as the dependent variable:
dependent_variable = alkaloid_class
dependent_variable_name = "alkaloid_class"

## Setting alkaloid classes in ants, with tricyclics as the dependent variable:
dependent_variable = alkaloid_class_in_ants
dependent_variable_name = "alkaloid_class_in_ants"

## Setting ant composition as the dependent variable:
dependent_variable = ants
dependent_variable_name = "ants"

## Setting genetic distance as the dependent variable:
dependent_variable = gen
dependent_variable_name = "gen"

## Part 5: Running MMRR analyses:

## nperm does not need to be specified, default is nperm = 999):
MMRR_output = MMRR(dependent_variable, predictor_list, nperm = 10000)

## Saving MMRR outputs as a text file:
capture.output(MMRR_output, file = paste0(path, "OUTPUTS_MMRR/MMRR_output_", dependent_variable_name, "_by_", predictor_name, ".txt"))

## Part 6: Plotting pairwise distances and regression line:
## The line is for visualization purposes only; significance was estimated based on MMRR, above.

## Creating .png file to save plot into:
png(filename = paste0(path, "OUTPUTS_MMRR/MMRR_output_", dependent_variable_name, "_by_", predictor_name, ".png"))

## Plotting:
plot = matplot(predictor, dependent_variable, type = "p", lty = 1, lwd = 1, lend = par("lend"), pch = 20, col = 1:1, cex = NULL, bg = NA, xlab = predictor_name, ylab = dependent_variable_name)
dependent_variable_dist = as.dist(dependent_variable) # setting dependent variable as dist to plot line
predictor_dist = as.dist(predictor) # setting predictor variable as dist to plot line
regLine(lm(dependent_variable_dist ~ predictor_dist, data = plot), lwd = 3, col = palette()[1]) # plotting line

## Finish plotting and save into .pdf file:
dev.off()

## Done!
