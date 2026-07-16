


files <- list.files(
  path = "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/4_Pred_ROMs/",
  pattern = "xml$",
  recursive = TRUE,
  full.names = TRUE
)

unlink(files)





# load fuctions
rm(list=ls())
library(rlist)
library(gbm)
library(dismo)
library(rgdal)
library(maptools)
library(maps)
library(mapdata)
library(raster) 
library(ncdf4)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)

range01 <- function(r){
  r.min = cellStats(r, "min")
  r.max = cellStats(r, "max")
  r.scale <- ((r - r.min) / (r.max - r.min))
  return(r.scale) #(r-rmin)/(rmax-rmin)
}
setwd("/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/4_Pred_ROMs/")   
template=raster("/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/4_Pred_ROMs/2002-01-01/tos.grd")
plot(template)

outDir="/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/5_Predictions_ROMs/"




#read model output
Species="SWOR"
modrep=readRDS("/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/3_Output/SWOR.res1.tc3.lr03.single_gridded.rds") # 10 unique models
spname="swor"


#predict function
predCIs_ROMS<-function(get_date,spname,modrep,stack,template,outDir,studyarea,droppath){
  
  spDir=paste0(outDir,spname,"/")
  if(!file.exists(spDir)){dir.create(spDir)}
  
  
  # stack=create_ROMS_daily_stack(get_date = get_date, template=template,staticDir=staticDir, predDir = predDir)
  #stack=as.data.frame(stackfile[[1]],stringsAsFactors=F) ## predicts on data.frame, not on a stack
  #head(stack)
  #colnames(stack) <- c("chl", "deptho", "o2", "so_cglo", "thetao_cglo", "zos_cglo")
 
  #mod_pred10 <- predict.gbm(modrep,newdata=stack,n.trees=1000,type='response')
  #mod_pred10s <- do.call(cbind,lapply(mod_pred10,data.frame,stringsAsFactors=FALSE))
  #colnames(mod_pred10s) <- as.character(seq(1,1,by=1))
  ## mean prediction over 10 runs
  #meanPred <- mod_pred10
  
  #stack$pred=meanPred
  #head(stack)
  
  
  
  names(stackfile)[1] <- 'bbv'
  names(stackfile)[2] <- 'deptho_sd'
  names(stackfile)[3] <- 'deptho'
  names(stackfile)[4] <- 'eke'
  names(stackfile)[5] <- 'ild'
  names(stackfile)[6] <- 'moon_phase'
  names(stackfile)[7] <- 'sos'
  names(stackfile)[8] <- 'ssh_sd'
  names(stackfile)[9] <- 'ssh'
  names(stackfile)[10] <- 'ssu_rotate'
  names(stackfile)[11] <- 'ssv_rotate'
  names(stackfile)[12] <- 'tos_sd'
  names(stackfile)[13] <- 'tos'
  
  pred.raster <- predict(stackfile,  modrep , type = "response")
  plot(pred.raster)
  maps::map('world2',add=TRUE,col=grey(0.7),fill=T)
  
  # Load the shapefile as a SpatVector
  #shp_vect <- vect(yft)
  # Convert the SpatVector to Spatial object
  #shp_spatial <- as(shp_vect, "Spatial")
  #shp_spatial <- spTransform(shp_spatial, crs(pred.raster))

 
  # Mask the raster with the shapefile
  #pred.raster <- mask(pred.raster, shp_spatial)
  
 
 
  

  ## standard error
  #sdPred <- apply(mod_pred10s,1,sd)
  #sePred <- sdPred/sqrt(ncol(mod_pred10s))
  
  ## confidence intervals
  #lowPred <- meanPred - 1.96*sePred
  #highPred <- meanPred + 1.96*sePred
  
 
  ## make rasters 
 
  #meanPredR <- setValues(template,meanPred)
  #plot(meanPredR)
  #lowPredR <- setValues(template,lowPred)%>%mask(.,studyarea)
  #highPredR <-setValues(template,highPred)%>%mask(.,studyarea)
  #sePredR <- setValues(template,sePred)%>%mask(.,studyarea)

  
  ## write rasters 
  writeRaster(pred.raster,paste0(spDir,spname,"_",get_date,"_mean"),outDir=outDir,overwrite=T)
  #writeRaster(lowPredR,paste0(spDir,spname,"_",get_date,"_lowCI"),overwrite=T)
  #writeRaster(highPredR,paste0(spDir,spname,"_",get_date,"_highCI"),overwrite=T)
  #writeRaster(sePredR,paste0(spDir,spname,"_",get_date,"_se"),overwrite=T)
  
  ## write PNGS
  make_png_operationalization(r=pred.raster,spname=spname,get_date=get_date,outDir=outDir,type="mean")
  #make_png_operationalization(r=lowPredR,spname=spname,get_date=get_date,outDir=outDir,type="lowCI")
  #make_png_operationalization(r=highPredR,spname=spname,get_date=get_date,outDir=outDir,type="highCI")
  #make_png_operationalization_SE(r=sePredR,spname=spname,get_date=get_date,outDir=outDir,type="se")
  
} # -------------------> final ecoroms predicion script

#get dates and predict
get_date<-list.dirs(getwd(),recursive = F, full.names=F)
flist <- list()
stackfile <- list()
for (i in 1:length(get_date)){  
  flist[[i]] <-  list.files(get_date[i], recursive = TRUE, full.names = T,  pattern = "tos.grd|tos_sd.grd|ssh.grd|ssh_sd.grd|sos.grd|deptho.grd|deptho_sd.grd|ssu_rotate.grd|ssv_rotate.grd|ild.grd|bbv.grd|EKE.grd|lunar_illumination.grd$")  
  
  #stackfile[[i]] <- stack(flist[[i]])
  stackfile <- stack(flist[[i]])
 
  
  predCIs_ROMS(get_date = get_date[i],modrep=modrep,spname = spname,studyarea=studyarea,stack=stack,template=template,outDir = outDir,droppath=droppath)
  
  print(get_date[i])
}





