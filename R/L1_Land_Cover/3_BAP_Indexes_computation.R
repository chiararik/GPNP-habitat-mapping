# BAP images processing:
# 1. create spectral indexes


rm(list = ls())
# Packages ----
list.of.packages <- c(
  "devtools",
  "terra",
  "data.table",
  "ff",
  "lubridate",
  "dplyr",
  "ggplot2",
  "ggspatial",
  "gridExtra",
  "ggpubr",
  "tidyterra",
  "RStoolbox",
  "rasterdiv",
  "snow",
  'doParallel'
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(package.i,character.only = TRUE))
}


# Folders ----
sito <- "PNGP"
root <- file.path(paste0("D:/ABRESO/",sito))
bap.folder <- file.path(root,"Landsat/level2/BAP/SG")
l3.folder <- file.path(root,"Landsat/level3/BAP_GF")

# Extract indexes ----
setwd(bap.folder)
imgs.list <- list.files(pattern="*_GS_195028_02_T1_SR_BAP_SWIR2*")
f <- function(x){as.Date(substr(x,1,4),format="%Y")}
yy.list <- lapply(imgs.list,FUN = f)
days.list <- sort(as.Date(unlist(yy.list)))
first.date <- min(days.list)
last.date <- max(days.list)
yy <- (strftime(first.date, format = "%Y"):strftime(last.date, format = "%Y"))
#cond <- sapply(yy, function(x) { substr(tail(x,1),nchar(tail(x,1)),nchar(tail(x,1))) == 0 | 
    #substr(tail(x,1),nchar(tail(x,1)),nchar(tail(x,1))) == 5})
#yy <- yy[cond]

for (y in yy) {
  
  # Growing season ----
  setwd(bap.folder)
  blue <- rast(list.files(pattern=glob2rx(paste0(y,"*GS*blue*"))))/10000
  green <- rast(list.files(pattern=glob2rx(paste0(y,"*GS*green*"))))/10000
  red <- rast(list.files(pattern=glob2rx(paste0(y,"*GS*red*"))))/10000
  NIR <- rast(list.files(pattern=glob2rx(paste0(y,"*GS*NIR*"))))/10000
  SWIR1 <- rast(list.files(pattern=glob2rx(paste0(y,"*GS*SWIR1*"))))/10000
  SWIR2 <- rast(list.files(pattern=glob2rx(paste0(y,"*GS*SWIR2*"))))/10000
  
  setwd(l3.folder)
  
  CVI<-(NIR*(red/((green)^2)))*1000
  #CVI[CVI < -1000]=-1000
  #CVI[CVI > 1000]=1000
  writeRaster(CVI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_CVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  CI<-((red-blue)/red)*10000
  CI[CI < -10000]=-10000
  CI[CI > 10000]=10000
  writeRaster(CI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_CI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  EVI<-(2.5*((NIR-red)/(NIR+6*red-7.5*blue+1)))*10000
  EVI[EVI < -10000]=-10000
  EVI[EVI > 10000]=10000
  writeRaster(EVI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_EVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  GVMI<-(((NIR+0.1)-(SWIR2+0.02))/((NIR+0.1)+(SWIR2+0.02)))*10000
  GVMI[GVMI < -10000]=-10000
  GVMI[GVMI > 10000]=10000
  writeRaster(GVMI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_GVMI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  GLI<-((2*green-red-blue)/(2*green+red+blue))*10000
  GLI[GLI < -10000]=-10000
  GLI[GLI > 10000]=10000
  writeRaster(GLI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_GLI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  mARI<- ((green^(-1)-red^(-1))*NIR)*10000
  mARI[mARI < -10000]=-10000
  mARI[mARI > 10000]=10000 
  writeRaster(mARI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_mARI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  MCARI1<- (1.2*(2.5*(NIR-red)-1.3*(NIR-green)))*10000
  MCARI1[MCARI1 < -10000]=-10000
  MCARI1[MCARI1 > 10000]=10000
  writeRaster(MCARI1,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_MCARI1.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  MSAVI <- ((2*NIR+1-sqrt((2*NIR+1)^2-8*(NIR-red)))/2)*10000
  MSAVI[MSAVI< -10000]=-10000
  MSAVI[MSAVI> 10000]=10000
  writeRaster(MSAVI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_MSAVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NDVI<- ((NIR-red)/(NIR+red))*10000
  NDVI[NDVI < -10000]=-10000
  NDVI[NDVI > 10000]=10000
  writeRaster(NDVI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_NDVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  GNDVI<- ((NIR-green)/(NIR+green))*10000
  GNDVI[GNDVI< -10000]=-10000
  GNDVI[GNDVI> 10000]=10000
  writeRaster(GNDVI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_GNDVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NBR<- ((NIR-SWIR1)/(NIR+SWIR1))*10000
  NBR[NBR< -10000]=-10000
  NBR[NBR> 10000]=10000
  writeRaster(NBR,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_NBR.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  SAVI<- (((NIR-red)/(NIR+red+0.5))*(1+0.5))*10000
  SAVI[SAVI< -10000]=-10000
  SAVI[SAVI> 10000]=10000
  writeRaster(SAVI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_SAVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  TSAVI<-((0.7*(NIR-0.7*red-(-0.32)))/(-0.32*NIR+red-(-0.32)*0.7+1.5*(1+(0.7)^2)))*10000
  TSAVI[TSAVI< -10000]=-10000
  TSAVI[TSAVI> 10000]=10000
  writeRaster(TSAVI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_TSAVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  SPVI<-(0.4*(3.7*(NIR-red)-1.2*abs(green-red)))*10000
  writeRaster(SPVI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_SPVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  SBI<-(0.3037*blue+0.2793*green+0.4743*red+0.5585*NIR+0.5082*SWIR1+0.1863*SWIR2)*10000
  writeRaster(SBI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_SBI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  GVI<-(-0.2848*blue-0.2435*green-0.5436*red+0.7243*NIR+0.0840*SWIR1-0.1800*SWIR2)*10000
  writeRaster(GVI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_GVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  WET<-(0.1509*blue+0.1973*green+0.3279*red+0.3406*NIR-0.7112*SWIR1-0.4572*SWIR2)*10000
  writeRaster(WET,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_WET.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  WDRVI<- ((0.1*NIR-red)/(0.1*NIR+red))*10000
  WDRVI[WDRVI>10000]<- 10000
  WDRVI[WDRVI< -10000]<- -10000
  writeRaster(WDRVI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_WDRVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  RDI<- (SWIR2/NIR)*10000
  #RDI[RDI>10000]<- 10000
  RDI[RDI<0]<-0
  writeRaster(RDI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_RDI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NDBI<- ((SWIR1 - NIR)/(SWIR1 + NIR))*10000
  NDBI[NDBI>10000]<- 10000
  NDBI[NDBI< -10000]<- -10000
  writeRaster(NDBI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_NDBI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  MAVI<- ((NIR-red)/(NIR+red+SWIR1))*10000
  MAVI[MAVI>10000]<- 10000
  MAVI[MAVI< -10000]<- -10000
  writeRaster(MAVI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_MAVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  ARVI<- ((NIR-red-1*(red-blue))/(NIR+red-1*(red-blue)))*10000
  ARVI[ARVI>10000]<- 10000
  ARVI[ARVI< -10000]<- -10000
  writeRaster(ARVI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_ARVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NDWI<- ((green - NIR)/(green + NIR))*10000
  NDWI[NDWI>10000]<- 10000
  NDWI[NDWI< -10000]<- -10000
  writeRaster(NDWI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_NDWI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  BSI<- (((SWIR2+red)-(NIR+blue))/((SWIR2+red)+(NIR+blue)))*10000
  BSI[BSI>10000]<- 10000
  BSI[BSI< -10000]<- -10000
  writeRaster(BSI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_BSI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NDMI<- ((NIR - SWIR2)/(NIR + SWIR2))*10000
  NDMI[NDMI>10000]<- 10000
  NDMI[NDMI< -10000]<- -10000
  writeRaster(NDMI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_NDMI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NBR2<- ((SWIR1-SWIR2) / (SWIR1 + SWIR2))*10000
  NBR2[NBR2>10000]<- 10000
  NBR2[NBR2< -10000]<- -10000
  writeRaster(NBR2,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_NBR2.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NDSI<- ((green-SWIR1) / (green + SWIR1))*10000
  NDSI[NDSI>10000]<- 10000
  NDSI[NDSI< -10000]<- -10000
  writeRaster(NDMI,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_NDSI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  
  rm(blue,red,green,NIR,SWIR1,SWIR2,
     CVI,CI,EVI,GVMI,GLI,mARI,MCARI1,
     MSAVI,NDVI,GNDVI,NBR,SAVI,TSAVI,SPVI,SBI,
     GVI,WET,WDRVI,RDI,NDBI,MAVI,ARVI,NDWI,BSI,
     NDMI,NBR2,NDSI)
  
  # Senescence season ----
  setwd(bap.folder)
  blue <- rast(list.files(pattern=glob2rx(paste0(y,"*SS*blue*"))))/10000
  green <- rast(list.files(pattern=glob2rx(paste0(y,"*SS*green*"))))/10000
  red <- rast(list.files(pattern=glob2rx(paste0(y,"*SS*red*"))))/10000
  NIR <- rast(list.files(pattern=glob2rx(paste0(y,"*SS*NIR*"))))/10000
  SWIR1 <- rast(list.files(pattern=glob2rx(paste0(y,"*SS*SWIR1*"))))/10000
  SWIR2 <- rast(list.files(pattern=glob2rx(paste0(y,"*SS*SWIR2*"))))/10000
  
  setwd(l3.folder)
  
  CVI<-(NIR*(red/((green)^2)))*1000
  #CVI[CVI < -1000]=-1000
  #CVI[CVI > 1000]=1000
  writeRaster(CVI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_CVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  CI<-((red-blue)/red)*10000
  CI[CI < -10000]=-10000
  CI[CI > 10000]=10000
  writeRaster(CI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_CI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  EVI<-(2.5*((NIR-red)/(NIR+6*red-7.5*blue+1)))*10000
  EVI[EVI < -10000]=-10000
  EVI[EVI > 10000]=10000
  writeRaster(EVI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_EVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  GVMI<-(((NIR+0.1)-(SWIR2+0.02))/((NIR+0.1)+(SWIR2+0.02)))*10000
  GVMI[GVMI < -10000]=-10000
  GVMI[GVMI > 10000]=10000
  writeRaster(GVMI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_GVMI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  GLI<-((2*green-red-blue)/(2*green+red+blue))*10000
  GLI[GLI < -10000]=-10000
  GLI[GLI > 10000]=10000
  writeRaster(GLI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_GLI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  mARI<- ((green^(-1)-red^(-1))*NIR)*10000
  mARI[mARI < -10000]=-10000
  mARI[mARI > 10000]=10000 
  writeRaster(mARI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_mARI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  MCARI1<- (1.2*(2.5*(NIR-red)-1.3*(NIR-green)))*10000
  MCARI1[MCARI1 < -10000]=-10000
  MCARI1[MCARI1 > 10000]=10000
  writeRaster(MCARI1,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_MCARI1.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  MSAVI <- ((2*NIR+1-sqrt((2*NIR+1)^2-8*(NIR-red)))/2)*10000
  MSAVI[MSAVI< -10000]=-10000
  MSAVI[MSAVI> 10000]=10000
  writeRaster(MSAVI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_MSAVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NDVI<- ((NIR-red)/(NIR+red))*10000
  NDVI[NDVI < -10000]=-10000
  NDVI[NDVI > 10000]=10000
  writeRaster(NDVI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_NDVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  GNDVI<- ((NIR-green)/(NIR+green))*10000
  GNDVI[GNDVI< -10000]=-10000
  GNDVI[GNDVI> 10000]=10000
  writeRaster(GNDVI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_GNDVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NBR<- ((NIR-SWIR1)/(NIR+SWIR1))*10000
  NBR[NBR< -10000]=-10000
  NBR[NBR> 10000]=10000
  writeRaster(NBR,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_NBR.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  SAVI<- (((NIR-red)/(NIR+red+0.5))*(1+0.5))*10000
  SAVI[SAVI< -10000]=-10000
  SAVI[SAVI> 10000]=10000
  writeRaster(SAVI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_SAVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  TSAVI<-((0.7*(NIR-0.7*red-(-0.32)))/(-0.32*NIR+red-(-0.32)*0.7+1.5*(1+(0.7)^2)))*10000
  TSAVI[TSAVI< -10000]=-10000
  TSAVI[TSAVI> 10000]=10000
  writeRaster(TSAVI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_TSAVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  SPVI<-(0.4*(3.7*(NIR-red)-1.2*abs(green-red)))*10000
  writeRaster(SPVI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_SPVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  SBI<-(0.3037*blue+0.2793*green+0.4743*red+0.5585*NIR+0.5082*SWIR1+0.1863*SWIR2)*10000
  writeRaster(SBI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_SBI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  GVI<-(-0.2848*blue-0.2435*green-0.5436*red+0.7243*NIR+0.0840*SWIR1-0.1800*SWIR2)*10000
  writeRaster(GVI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_GVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  WET<-(0.1509*blue+0.1973*green+0.3279*red+0.3406*NIR-0.7112*SWIR1-0.4572*SWIR2)*10000
  writeRaster(WET,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_WET.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  WDRVI<- ((0.1*NIR-red)/(0.1*NIR+red))*10000
  WDRVI[WDRVI>10000]<- 10000
  WDRVI[WDRVI< -10000]<- -10000
  writeRaster(WDRVI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_WDRVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  RDI<- (SWIR2/NIR)*10000
  #RDI[RDI>10000]<- 10000
  RDI[RDI<0]<-0
  writeRaster(RDI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_RDI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NDBI<- ((SWIR1 - NIR)/(SWIR1 + NIR))*10000
  NDBI[NDBI>10000]<- 10000
  NDBI[NDBI< -10000]<- -10000
  writeRaster(NDBI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_NDBI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  MAVI<- ((NIR-red)/(NIR+red+SWIR1))*10000
  MAVI[MAVI>10000]<- 10000
  MAVI[MAVI< -10000]<- -10000
  writeRaster(MAVI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_MAVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  ARVI<- ((NIR-red-1*(red-blue))/(NIR+red-1*(red-blue)))*10000
  ARVI[ARVI>10000]<- 10000
  ARVI[ARVI< -10000]<- -10000
  writeRaster(ARVI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_ARVI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NDWI<- ((green - NIR)/(green + NIR))*10000
  NDWI[NDWI>10000]<- 10000
  NDWI[NDWI< -10000]<- -10000
  writeRaster(NDWI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_NDWI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  BSI<- (((SWIR2+red)-(NIR+blue))/((SWIR2+red)+(NIR+blue)))*10000
  BSI[BSI>10000]<- 10000
  BSI[BSI< -10000]<- -10000
  writeRaster(BSI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_BSI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NDMI<- ((NIR - SWIR2)/(NIR + SWIR2))*10000
  NDMI[NDMI>10000]<- 10000
  NDMI[NDMI< -10000]<- -10000
  writeRaster(NDMI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_NDMI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NBR2<- ((SWIR1-SWIR2) / (SWIR1 + SWIR2))*10000
  NBR2[NBR2>10000]<- 10000
  NBR2[NBR2< -10000]<- -10000
  writeRaster(NBR2,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_NBR2.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  NDSI<- ((green-SWIR1) / (green + SWIR1))*10000
  NDSI[NDSI>10000]<- 10000
  NDSI[NDSI< -10000]<- -10000
  writeRaster(NDMI,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_NDSI.tif"),overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
  
  
  rm(blue,red,green,NIR,SWIR1,SWIR2,
     CVI,CI,EVI,GVMI,GLI,mARI,MCARI1,
     MSAVI,NDVI,GNDVI,NBR,SAVI,TSAVI,SPVI,SBI,
     GVI,WET,WDRVI,RDI,NDBI,MAVI,ARVI,NDWI,BSI,
     NDMI,NBR2,NDSI)
  
}

