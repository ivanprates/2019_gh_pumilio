############################################
## R script by Ivan Prates (ivanprates.org), Smithsonian National Museum of Natural History, Washington, DC, USA, June 2019.
## This script does the following:
## Performs a linear regression on alkaloid (or alkaloid class) richness versus ant richness.

## PART 1: Getting ready:

## Loading packages:
library(car)

## Selecting path depending on which computer I'm working from:
path = "~/Dropbox (Smithsonian)/2018_pumilio/2018-07/"

## Making folder to save the files we'll generate:
dir.create(paste0(path, "OUTPUTS_linear_regression/"))

## Importing counts for alkaloid data:
alkaloid_counts = read.csv(file = paste0(path, "DATA_alkaloids_saporito2007/2018-07_alkaloid_counts.csv"), header = TRUE, row.names = "alk_siteID")
alkaloid = alkaloid_counts[1:46, 1]
alkaloid_class = alkaloid_counts[1:46, 2]

## Importing counts for ant data:
ant_counts = read.csv(file = paste0(path, "DATA_alkaloid_bearing_ants/2018-07_presabs_ants_T10_counts.csv"), header = TRUE, row.names = "alk_siteID")
ants = ant_counts[1:46, 1]

## Part 2: Selecting variables:

## Using ants as an independent variable:
predictor = ants
predictor_name = "ants"

## Using alkaloid classes as an independent variable:
#predictor = alkaloid_class
#predictor_name = "alkaloid_class"

## Setting individual alkaloids as the dependent variable:
dependent_variable = alkaloid
dependent_variable_name = "alkaloid"

## Setting alkaloid classes as the dependent variable:
dependent_variable = alkaloid_class
dependent_variable_name = "alkaloid_class"

## Part 3: Running linear regressions:

## Run the linear regression:
linear_regression = lm(dependent_variable ~ predictor)
results = summary(linear_regression)

## Saving outputs as a text file:
capture.output(results, file = paste0(path, "OUTPUTS_linear_regression/linear_regression_", dependent_variable_name, "_by_", predictor_name, ".txt"))

## Creating .png file to save plot into:
png(filename = paste0(path, "OUTPUTS_linear_regression/linear_regression_", dependent_variable_name, "_by_", predictor_name, ".png"))

## Plotting:
plot = matplot(predictor, dependent_variable, type = "p", lty = 1, lwd = 1, lend = par("lend"), pch = 20, col = 1:1, cex = NULL, bg = NA, xlab = predictor_name, ylab = dependent_variable_name)
regLine(lm(dependent_variable ~ predictor, data = plot), lwd = 3, col = palette()[1]) # plotting line

## Finish plotting and save into .pdf file:
dev.off()

## Done!
