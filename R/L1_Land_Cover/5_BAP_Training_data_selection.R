# Training data building and predictors selection

# Cleaning of the LULC (smaller patches and classes)
# Selection of the purest pixels through z-statistics
# Tests for predictor selection based on the sample size selected as training
# All computed at Landsat resolution (30 m)

################################ Libraries ################################
library(randomForest)
library(psych)
library(sp)
#library(rgdal)
library(terra)
library(e1071)
library(raster)
library(lattice)
library(ggplot2)
library(caret)
library(Boruta)
library(corrplot)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(data.table)


# Folders ##################################################################
sito <- "PNGP"
root <- file.path(paste0("D:/ABRESO/",sito))
bap.folder <- file.path(root,"Landsat/level2/BAP/SG")
L2folder <- file.path(root,"Landsat/level2/AC")
l3.folder <- file.path(root,"Landsat/level3/BAP_GF")

#### Retrieve available images and scores  ----
imgs.list <- list.files(path = L2folder,pattern="T1")
f <- function(x){as.Date(substr(x,18,25),format="%Y%m%d")}
days.list <- lapply(imgs.list,FUN = f)
days.list <- sort(as.Date(unlist(days.list)))

first.date <- min(days.list)
last.date <- max(days.list)

yy <- (strftime(first.date, format = "%Y"):strftime(last.date, format = "%Y"))


sel.pred <- c("L3","GS_BSI","GS_CI","GS_GLI","GS_mARI","GS_NBR2","GS_NDBI","GS_NDWI",
              "GS_SBI","SS_ARVI","SS_BSI","SS_CI","SS_CVI","SS_EVI","SS_GLI",
              "SS_GVI","SS_GVMI","SS_mARI","SS_MCARI1","SS_NBR2","SS_NDBI",
              "SS_WET","GS_NIR","GS_SWIR1","SS_blue","SS_NIR","SS_SWIR1",
              "elev","roughness","slope","twi","Geologia_5m")

sel.pred.names <- c("LULC","GS_BSI","GS_CI","GS_GLI","GS_mARI","GS_NBR2","GS_NDBI","GS_NDWI",
                    "GS_SBI","SS_ARVI","SS_BSI","SS_CI","SS_CVI","SS_EVI","SS_GLI",
                    "SS_GVI","SS_GVMI","SS_mARI","SS_MCARI1","SS_NBR2","SS_NDBI",
                    "SS_WET","GS_NIR","GS_SWIR1","SS_blue","SS_NIR","SS_SWIR1",
                    "Elevation","Roughness","Slope","TWI","Geology")


aoi <- vect("D:/ABRESO/PNGP/Shapefile/PNGP.shp")

for (y in yy) {
  
  ################################ Import data ################################
  #### Import spectral predictors ####
  setwd(l3.folder)
  indlist <- list.files(getwd(), pattern=glob2rx(paste0(y,'*.tif$')), all.files=TRUE, full.names=TRUE)
  ind.names <- lapply(indlist,function(x) {paste0(substr(x,43,44),"_",substr(x,66,nchar(x)-4))})
  indices<- rast(indlist) 
  names(indices) <- unlist(ind.names)
  
  setwd(bap.folder)
  bandslist <- list.files(getwd(), pattern=glob2rx(paste0(y,'*.tif$')), all.files=TRUE, full.names=TRUE)
  #bandslist <- bandslist[-c(5,12)] # To remove scoring 
  bands.names <- lapply(bandslist,function(x) {paste0(substr(x,43,44),"_",substr(x,66,nchar(x)-11))})
  bands<- rast(bandslist)
  
  #### Import ground references ####
  lulc <- rast("D:/ABRESO/PNGP/LULC/Training/LULC_L3_5m.tif") %>% 
    mask(aoi) %>% resample(indices, method="near") 
  lulc[lulc==0] <- NA
  
  ### Import static predictors ####
  setwd("D:/ABRESO/PNGP/Ancillary/5m")
  predlist <- list.files(getwd(), pattern=glob2rx(paste0('*.tif$')), all.files=TRUE, full.names=TRUE)
  #predlist <- predlist[-6] # TO remove Geology
  pred.names <- lapply(predlist,function(x) {substr(x,29,nchar(x)-4)})
  predictors <- rast(predlist) %>% resample(indices, method="near") 
  
  all.names <- c("L3",ind.names,bands.names,pred.names)
  all.pred <- c(lulc,indices,bands,predictors)
  names(all.pred) <- unlist(all.names)
  
  all.pred <- subset(all.pred, sel.pred) %>% mask(aoi)
  
  ######################### Subset training #################################
  
  # Compute z-statistics ----
  # https://doi.org/10.1016/j.rse.2015.12.031
  
  df <- as.data.frame(all.pred)
  df <- na.omit(df)
  
  df$L3 <- as.factor(df$L3)
  gc()
  
  df <- df[,1:31]
  
  # for each class compute statistics 
  for (f in levels(df$L3)[1:length(levels(df$L3))]) {
    
    df.sub <- subset(df, df$L3==f)
    df.sub <- df.sub[,2:ncol(df)]
    
    data_long <- tidyr::gather(df.sub, factor_key=TRUE)
    values <- data_long %>% group_by(key)%>%
      summarise(mean= mean(value), sd= sd(value)) #, max = max(value),min = min(value)
    
    lulc.f <- all.pred[[1]]
    lulc.f[lulc.f !=as.numeric(f)] <- NA
    
    coverage <- nrow(df.sub)*(xres(lulc))^2
    print(paste("class",f,"=",coverage,"m2 -",ncell(lulc.f[lulc.f==as.numeric(f)]),"pixels"))
    tot.cell <- ncell(lulc.f[lulc.f==as.numeric(f)])
    
    n <- 2
    r <- 1
    for (l in 1:(nlyr(all.pred)-1)) {
      l <- all.pred[[n]]
      m <- as.numeric(values[r,2])
      ds <- as.numeric(values[r,3])
      z.stat.part <- ((l-m)/ds)^2
      z.stat.part <- mask(z.stat.part,lulc.f)
      writeRaster(z.stat.part,paste0(root,"/Landsat/z_statistics/",y,"_class_",f,"_variable_", names(l),"_zstat_30m_8cl.tif"),
                  overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
      n <- n+1
      r <- r+1
      rm(l,z.stat.part)
    }
    
    # compute final z-stat
    setwd(paste0(root,"/Landsat/z_statistics"))
    z.stats <- rast(list.files(pattern=glob2rx(paste0(y,"_class_",f,"_variable_*_zstat_30m_8cl.tif"))))
    f.z.stat <- function(x){return(sqrt(sum(x,na.rm=TRUE)))}
    z.stats <- app(z.stats, fun=f.z.stat, cores=8)
    z.stats <- mask(z.stats,lulc.f)
    z.stats[is.infinite(z.stats)] <- NA
    writeRaster(z.stats, filename=paste0(root,"/Landsat/z_statistics/",y,"_class_",f,"_TOT_zstat_30m_8cl.tif"), 
                overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
    
    # Analyse each class
    sample.size <- ncell(z.stats[!is.na(z.stats)])
    hist(z.stats,
         main = paste("Class",f),
         xlab = "z-value", ylab = "Frequency",
         col = "springgreen")
    #df.stat <- as.data.frame(z.stats) %>% 
    #na.omit()
    #define histogram
    #my_plot <- ggplot(df.stat, aes(x=lyr.1)) + geom_histogram()
    #define table
    #stat.sum <- tableGrob(terra::summary(z.stats, size=ncell(z.stats)))
    #create scatterplot and add table underneath it
    #info <- grid.arrange(my_plot, stat.sum, top=textGrob(paste("Class",f,"- sample size:",sample.size,"pixels")))
    #info <- grid.arrange(arrangeGrob(my_plot, stat.sum, top=(paste("Class",f,"- sample size:",sample.size,"pixels")), ncol=2))
    #ggsave(paste(y,"_class",f,"_30m.jpeg"), info)
    
    # Select pixels falling in μ ± σ of z-score
    
    if (sample.size > 50) {
      μ <- as.numeric(global(z.stats, fun="mean", na.rm=TRUE))
      σ <- as.numeric(global(z.stats, fun="sd", na.rm=TRUE))
      t1 <- μ - σ
      t2 <- μ + σ
      training <- app(z.stats, fun = function(x) {return(x > t1 & x < t2)})
      training <- mask(training, lulc.f)
      sel.cell <- ncell(training[training==1])
      sel.percentage <- round(sel.cell*100/tot.cell,0)
      print(paste("training class",f,"-",sel.cell,"pixels, manteined", sel.percentage,"%"))
      training[training==1] <- as.numeric(f)
      training[training==0] <- NA
      writeRaster(training, paste0(root,"/Landsat/z_statistics/Training/",y,"_class_",f,"_training_30m_8cl.tif"), datatype='INT2S', filetype = "GTiff", overwrite=TRUE)
      
    } else {
      print("Sample size too small")
      training <- app(z.stats, fun = function(x) {return(!is.na(x))})
      training[training==1] <- as.numeric(f)
      training[training==0] <- NA
      writeRaster(training, paste0(root,"/Landsat/z_statistics/Training/",y,"_class_",f,"_training_30m_8cl.tif"), datatype='INT2S', filetype = "GTiff", overwrite=TRUE)
    }
    rm(z.stats,df.sub,data_long,values,lulc.f,training,μ,σ,t1,t2)
    gc()
  }
  
}
