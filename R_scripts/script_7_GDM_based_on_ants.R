### Script to generate General Dissimilarity Model (GDM) based on ant species composition for Prates at al. 2019 Ecology & Evolution
### By Andrea Paz, 2019 #####


library(gdm)
library(rgdal)
setwd("Pumilio_ants/")
###load alkaloid data with geographic information (pixel coordinates)
alks<-read.csv("4_matrix_46pixels_1km-pixels_EXTRACT_ANTS_FOR_THESE.csv")
sppTab <- alks[,c(1:3,6:length(alks))]

##environmental data from rasters (ant models)
# load cont models that are pre arranged to pumilio distribution size
setwd("Models/")
model_list<-list.files(path="./",pattern='correct_')
envRast <- stack(model_list)
###load oophaga pumilio distribution and confirm extents are correct and NA values are in correct format
setwd("pumilio_distribution/")
pumilio_iucn<-readOGR(("."),"oophaga_pumilio")
envRast<-crop(envRast,pumilio_iucn)
envRast<-mask(envRast,pumilio_iucn)
for (i in 1:nlayers(envRast)){
listOfDataFrames<-extract(envRast[[i]], pumilio_iucn, cellnumbers = TRUE)
na.cells <- as.data.frame(do.call("rbind", listOfDataFrames))
envRast[[i]][na.cells[is.na(na.cells$value),]$cell] <- 0
}
##site pair formatting for GDM package
gdmTab <- formatsitepair(sppTab, bioFormat=1, XColumn="longitude", YColumn="latitude", siteCol="saporito_siteID", predData=envRast)
gdmTab[1:3,]
# make sure there are no NA values
gdmTab <- na.omit(gdmTab)
## Create a GDM for alkaloids based on ant distributions
gdm.1 <- gdm(gdmTab, geo=T)
gdm.1.pred <- predict(gdm.1, gdmTab)
plot(gdmTab$distance, gdm.1.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))



###WTo plot the GDM on a map
tabTrans <- gdm.transform(gdm.1, envRast)
rastDat <- na.omit(getValues(tabTrans))
pcaSamp <- prcomp(rastDat)
pcaRast <- predict(tabTrans, pcaSamp, index=1:3)
pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255
pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255
pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255
plotRGB(pcaRast, r=1, g=3, b=2)

###show model predictive power
str(gdm.1)

###individual variable contribution
varImp<-gdm.varImp(gdmTab,geo=F,fullModelOnly=F)

varImp[[2]]
varImp[[1]]

###Here based on the results run a final time with just the selected variables modifying the working directory with the models to:
###setwd("Models_contrib/")
