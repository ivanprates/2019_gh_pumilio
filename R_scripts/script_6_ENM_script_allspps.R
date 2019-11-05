
### Script to make SDM's of ants for Prates et al 2019 Ecology and Evolution 
### By Andrea Paz, 2019 ####
                      
options(java.parameters = "-Xmx14000m") ##Check this in every computer and set accoridng to available RAM
###Load required packages
library(ENMeval)
library(spThin)
library(rgeos)
library(sp)
library(rgdal)
library(raster)
#Setwd to the project folder, this must be adjusted in different computers
setwd("~/Andrea_Paz_data/Pumilio_models/")

#Load and stack predictor variables , in this case 19 bioclim variables from the WorldClim database
predictors <- stack(list.files(path="./Bioclim/",pattern='asc', full.names=TRUE))

###Load and prepare occurence data (long,lat for maxent) , 
    ##Must have three columns named Species, Latitude and Longitude in that order

species_localities<-read.csv("ants.csv",h=T)
species_list<-unique(species_localities$Species)

###loop through the species list creating a model for each one if possible n>5 else write the name of species to file
for (i in species_list){
  ##Create a shapefile for each species with the occurrence points and write to a file
Points<-subset(species_localities,Species==i)
ID<-c(1:length(Points$Latitude))
Points_occ<-data.frame(ID,Points[,c(3,2)])
dir<-getwd()
coordinates(Points_occ)<-c('Longitude','Latitude')
writeOGR(Points_occ, dir,layer=paste("Points",i,sep="_"), driver="ESRI Shapefile")
occ<-readOGR(dir, paste("Points",i,sep="_"))
occ<-as.data.frame(occ)
Species<-rep(i,length(occ[,1]))
occ<-data.frame(Species,occ[,c(2,3)])

#Locality data thinning using spThin package thin.par is distance in km
presencia_thinned<-thin(occ, lat.col = "coords.x2", long.col = "coords.x1", spec.col = "Species",
                       thin.par=5, reps=5, locs.thinned.list.return = TRUE, write.files = TRUE,
                      max.files = 5, out.dir=paste(i,"thinned_data_test",sep="_"), out.base = "thinned_data",
                       write.log.file = TRUE, log.file = "spatial_thin_log.txt",
                      verbose = TRUE)
occ1<-presencia_thinned[[1]]
if (length(occ1$Longitude)>=5){ ##for les than 5 points no model is created just a file with the name of the species
  tryCatch({
  #define background as MCP 
##MCP only  
backgd<-convHull(occ1)@polygons
#with buffer first project
		##Define projection and project
crs(backgd)<- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
new_coord<-"+proj=eqdc +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs "
backgd1<-spTransform(backgd, CRS(new_coord))
		##then buffer at 100km
backgd1 <- gBuffer(backgd1, width=100000)
		##then reproject
backgd<-spTransform(backgd1, crs(backgd))
bkg <- rasterize(backgd, predictors, field=1)
bkg1<-mask(bkg,predictors$bio_1)
bkg2<-crop(bkg1,backgd)
#Create background points
backg <- randomPoints(bkg2, n=10000) 
####Crop predictors to species specific background area
predictors1<-crop(predictors,backgd)
predictors1<-mask(predictors1,backgd)

###Generate models using ENMeval with different parameter combinations using randomkfold for n>20 else jackknife
if (length(occ1$Longitude)>=20){ 
  k<-5
modelo<-ENMevaluate(occ1,predictors1 , bg.coords = backg, occ.grp = NULL,
            bg.grp = NULL, RMvalues = seq(0.5, 4, 0.5),
            fc = c("L", "LQ", "H", "LQH", "LQHP"),
            categoricals = NULL, n.bg = NULL, method = "randomkfold",
            overlap = FALSE, aggregation.factor = c(2, 2),
            kfolds = k, bin.output = FALSE, clamp = TRUE,
            rasterPreds = T, parallel = T, numCores =3,algorithm="maxent.jar")
save(modelo,file=paste('ENM_eval_blocks_R',i,sep="_"))
}
else {
  modelo<-ENMevaluate(occ1,predictors1 , bg.coords = backg, occ.grp = NULL,
                      bg.grp = NULL, RMvalues = seq(0.5, 4, 0.5),
                      fc = c("L", "LQ", "H", "LQH", "LQHP"),
                      categoricals = NULL, n.bg = NULL, method = "jackknife",
                      overlap = FALSE, aggregation.factor = c(2, 2),
                      kfolds = k, bin.output = FALSE, clamp = TRUE,
                      rasterPreds = T, parallel = T, numCores = 3,algorithm="maxent.jar")
  save(modelo,file=paste('ENM_eval_jakk_R',i,sep="_"))
}


###pick final set of parameters producing model with lowest AIC 
best_model<-which (modelo@results$delta.AICc == 0)
best_model<-best_model[1]
model_pred<-predict(modelo@models[[best_model]],predictors1,args=c("outputformat=logistic"))
plot(model_pred)
points(occ1,pch=16,cex=0.6)
    ###Write model to file
writeRaster(model_pred,paste(i,"AIC.asc",sep="_"),format="ascii")


####Thresholding with T10
thresholds<-modelo@models[[best_model]]@results
t10<-thresholds['X10.percentile.training.presence.logistic.threshold',]
reclass_mat_t10<-matrix(c(0,t10,0,t10,1,1),ncol=3,byrow=T)
pred_thresholded_t10<-reclassify(model_pred,reclass_mat_t10)
    ###Write thresholded model to file 
writeRaster(pred_thresholded_t10,paste(i,"AIC_T10.asc",sep="_"),format="ascii")

####Thresholding with MPT
thresholds<-modelo@models[[best_model]]@results
mpt<-thresholds['Minimum.training.presence.logistic.threshold',]
reclass_mat_mpt<-matrix(c(0,mpt,0,mpt,1,1),ncol=3,byrow=T)
pred_thresholded_mpt<-reclassify(model_pred,reclass_mat_mpt)
###Write thresholded model to file 
writeRaster(pred_thresholded_mpt,paste(i,"AIC_MPT.asc",sep="_"),format="ascii")

  },error=function(e){cat("ERROR:",i,conditionMessage(e),"\n")})
}
else{
  print (paste(i,"cannot be modeled less than 5 points"))
  bla<-paste(i,"cannot be modeled less than 5 points")
  write.table(bla,"errors.txt",append=T)
  }
}

####After modelling all maps must be set to the same extent using the IUCN distribution map for Oophaga pumilio
 ###Set working directory to folder with shapefile of interest
setwd("pumilio_distribution/")

pumilio<-readOGR(dsn=".",layer="oophaga_pumilio")
##set working directory to folder with continiuous models, 
      
setwd("Models/")
model_list<-list.files(path="./",pattern='AIC.asc')
  ##extend or crop model to match the desired extent
for(i in model_list){
  model<-raster(i)
  model = crop(model, pumilio)
  model = extend(model, pumilio)
  ##Write model to ASCII file
  writeRaster(model,paste("correct",i,sep="_"),format="ascii")
}

      ### Compute richness , map and write raster
###tapply the same procedure (before) to thresholded ones
####Load and stack thresholded models of correct extent
setwd("Models_thresholded/")
model_list<-list.files(path="./",pattern='correct_')
models <- stack(model_list)
richness<-sum(models,na.rm=T)
writeRaster(richness,"ant_richness_crop_T10.asc",format="ascii")
plot(richness)

####  To create a matrix with alkaloids and ants for same pixels.



