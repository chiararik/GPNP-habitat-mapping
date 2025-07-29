# Composing of the BAP 
# as per the following publications:
# http://dx.doi.org/10.1080/07038992.2014.945827
# https://doi.org/10.1016/j.rse.2018.10.031

rm(list = ls())
# Libraries ----
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
  "RStoolbox"
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


# Parameters ----

sito <- "PNGP"    #"BROCON" or "VALGRANDE" or "PNPG"

# Folders ----
root <- file.path(paste0("D:/ABRESO/",sito))
L1folder <- file.path(root,"Landsat/level1")
L2folder <- file.path(root,"Landsat/level2/AC")
L2TC.folder <- file.path(root,"Landsat/level2/ATC/improved_cosine")
L3folder <- file.path(root,"Landsat/level3")
bap.folder <- file.path(root,"Landsat/level2/BAP")

#### Retrieve available images and scores  ----
imgs.list <- list.files(path = L2TC.folder,pattern="T1$")
f <- function(x){as.Date(substr(x,18,25),format="%Y%m%d")}
days.list <- lapply(imgs.list,FUN = f)
days.list <- sort(as.Date(unlist(days.list)))

first.date <- min(days.list)
last.date <- max(days.list)

yy <- (strftime(first.date, format = "%Y"):strftime(last.date, format = "%Y"))
#cond <- sapply(yy, function(x) { substr(tail(x,1),nchar(tail(x,1)),nchar(tail(x,1))) == 0 | 
    #substr(tail(x,1),nchar(tail(x,1)),nchar(tail(x,1))) == 5})
#yy <- yy[cond]

aoi <- vect(paste0(root,"/TAGLIO_",sito,"_T32TLR_S2_32632.shp"))
s1.w <- 1     # Distance to clouds and snow
s2.w <- 0.5   # Sensor 
s3.w <- 0.25  # Distance to central date
s4.w <- 0.5   # Haze-Optimized Transformation
s5.w <- 0.75  # Coverage
s6.w <- 1     # Snow free area

for (y in yy){
  
  # Growing season ----
  season1 <- interval(ymd(paste0(y,'-06-15')), ymd(paste0(y,'-08-31')))
  season1.imgs.list <- list()
  
  for (i in imgs.list){
    data.img <- as.Date(substr(i,18,25),format="%Y%m%d")
    if (data.img %within% season1) { 
      
      # variables ----
      sensore <- substr(i,1,2)
      satellite <- substr(i,3,4)
      level <- substr(i,6,9)
      path <- substr(i,11,13)
      row <- substr(i,14,16)
      year <- substr(i,18,21)
      month <- substr(i,22,23)
      day <- substr(i,24,25)
      cc <- substr(i,36,37)
      acquisitionTime <- substr(i,18,25)
      processingTime <- substr(i,27,34)
      tx <- substr(i,39,40)
      id <- substr(i,1,40)
      endname <- substr(i,10,40)
      
      setwd(file.path(L2TC.folder,i))
      
      if (file.exists(paste0("Snow&Cloud_mask_30m_",i,".tif"))==FALSE) {next}
      
      
      # Visualize data ----
      setwd(file.path(L2folder,i))
      
      QA <- rast(list.files(pattern=glob2rx(paste0("*QA_PIXEL*"))))
      
      # Visualize data ----
      
      
      if (satellite == "08" | satellite == "09") {
        rgb <- rast(c(list.files(pattern=glob2rx("*_B4.*")), 
                      list.files(pattern=glob2rx("*_B3.*")), 
                      list.files(pattern=glob2rx("*_B2.*")))) %>% 
          crop(aoi)
      } else {
        rgb <- rast(c(list.files(pattern=glob2rx("*_B3.*")), 
                      list.files(pattern=glob2rx("*_B2.*")), 
                      list.files(pattern=glob2rx("*_B1.*")))) %>% 
          crop(aoi)
      }
      
      plotRGB(rgb,r=1,g=2,b=3,axes=TRUE,mar = 2,main=i, scale=65535) #
      #rm(rgb)
      
      ### s1: distance to clouds/shadows & snow score: ----
      setwd(file.path(L2TC.folder,i))
      clouds <- rast(list.files(pattern=glob2rx(paste0("Snow&Cloud_mask_30m_",i,"*")))) %>% 
        app(fun = function(x) {return(x == 205)}) # retrieve cloud mask
      snow <- rast(list.files(pattern=glob2rx(paste0("Snow&Cloud_mask_30m_",i,"*")))) %>% 
        app(fun = function(x) {return(x == 100)}) # retrieve snow mask
      
      # Add also QA clouds
      #QA <- crop(QA,clouds)
      
      #if (satellite == "08" | satellite == "09"){
        #cloudShadowvalues <-c(21826,21890,22080,22144,22280,23888,23952,24088,24216,24344,24472)
        #cloudShadows <- app(QA,fun = function(x) {return(x %in% cloudShadowvalues)})
        
        #cirrusvalues <-c(54596,54852,55052)
        #cirrus <- app(QA,fun = function(x) {return(x %in% cirrusvalues)}) %>% 
          #mask(aoi)
        
        #f2<-function(x,y){return(x == 1 | y == 1)} 
        #cloud_pass0 <- lapp(c(cloudShadows,cirrus),fun=f2) %>% 
          #mask(aoi)
      #}
      
      #if (satellite == "04" | satellite == "05" | satellite == "07") {
        #cloudShadowvalues <-c(5442,5506,56965760,5896,7440,7568,7696,7824,7960,8088)
        #cloud_pass0 <- app(QA,fun = function(x) {return(x %in% cloudShadowvalues)}) %>% 
          #mask(aoi)
      #} 
      
      clouds.snow <- clouds+snow #cloud_pass0+
      clouds.snow[clouds.snow==2]<-1
      
      pixel_size = xres(clouds.snow)
      cdst_req <- 50
      
      background <- clouds.snow
      background[background == 0] <- NA
      dist_min <- distance(background,unit="m")/pixel_size # compute distance to clouds and snow for each pixel
      
      s1 <- dist_min
      s1[s1 < cdst_req | s1 > 0] <- 1./(1.+exp(-0.01*(s1-cdst_req/2))) # distances lower than 10 pixels -> function
      s1[s1 >= cdst_req] <- 1   # distance greater than 10 pixels -> 1
      s1[s1 == 0] <- 0    # cloud and shadows -> 0
      
      rm(background, dist_min, clouds.snow)
      
      ### Create cloud & NA mask ----
      na.mask <- rgb[[3]]
      na.mask <- terra::resample(na.mask,clouds)
      na.mask[!is.na(na.mask)] <- 0
      rm(rgb)
      
      mask <- clouds
      mask[mask != 1] <- 0
      mask[mask == 1] <- NA
      
      mask <- mask+na.mask
      
      ### s2: sensor score ----
      if (year>2003 & sensore=="07"){ sv <- 0.5 } else { sv <- 1 }
      s2 <- rast(ncol=ncol(s1), nrow=nrow(s1), ext=ext(s1), crs = crs(s1)) #create an empty rast
      values(s2) <- sv    
      
      ### s3: doy score ----
      t.doy <- 210 # end of July    -- 196 # mid July
      sl <- 9.3
      doy <- as.numeric(strftime(data.img, format = "%j"))
      v <- (1/(sl*sqrt(2*3.14)))*exp(-0.5*((doy-t.doy)/sl)^2)
      v.max <- (1/(sl*sqrt(2*3.14)))*exp(-0.5*((t.doy-t.doy)/sl)^2)
      v.doy <- v/v.max
      
      s3 <- rast(ncol=ncol(s1), nrow=nrow(s1), ext=ext(s1), crs = crs(s1)) #create an empty rast
      values(s3) <- v.doy   
      
      ### s4: haze-Optimized Transformation (HOT) ----
      # Get bands 
      setwd(file.path(L2TC.folder,i))
      if (satellite == "04" | satellite == "05" | satellite == "07") {
        blue <- rast(list.files(pattern=glob2rx("*B1.*"))) 
        red <- rast(list.files(pattern=glob2rx("*B3.*")))
      } else if (satellite == "08" | satellite == "09") {
        blue <- rast(list.files(pattern=glob2rx('*B2.*')))
        red <- rast(list.files(pattern=glob2rx('*B4.*')))
      }
      
      HOT <- blue-0.5*red-0.08
      s4 <- 1/(1+exp((10/0.02)*(HOT+0.075)))
      
      ### s5: coverage ----
      cfc <- ncell(mask[mask==0])/ncell(mask)
      s5 <- rast(ncol=ncol(s1), nrow=nrow(s1), ext=ext(s1), crs = crs(s1)) #create an empty rast
      values(s5) <- cfc
      
      ### s6: snow area ----
      sc <- ncell(snow[snow!=1])/ncell(snow)
      s6 <- rast(ncol=ncol(s1), nrow=nrow(s1), ext=ext(s1), crs = crs(s1)) #create an empty rast
      values(s6) <- sc
      
      
      ### Scores sum ----
      score.sum <- ((s1*s1.w)+(s2*s2.w)+(s3*s3.w)+#(s4*s4.w)+
                      (s5*s5.w)+(s6*s6.w)/(s1.w+s2.w+s3.w+s5.w+s6.w)) %>% #+s4.w
        mask(mask) 
      #score.sum[is.na(score.sum)] <- 0
      
      rm(mask)
      # Save score img
      score.sum <- round(score.sum*1000,0)
      #plot(score.sum)
      setwd(L2TC.folder)
      writeRaster(score.sum,paste0(i,"_SCORE.tif"),datatype="INT2S",filetype="GTiff",overwrite=TRUE)
      
      # Add to a list 
      #season1.imgs.list <- append(season1.imgs.list,i)
      season1.imgs.list[[length(season1.imgs.list) + 1]] = i
      
      rm(score.sum,s1,s2,s3,s4,s5) 
    }
  }
  
  setwd(L2TC.folder)
  write.table(season1.imgs.list, file = paste0(y,"_season1.imgs.list.txt"), sep = ",")
  
  # Senescence season ----
  season2 <- interval(ymd(paste0(y,'-09-15')), ymd(paste0(y,'-11-30')))
  season2.imgs.list <- list()
  
  for (i in imgs.list){
    data.img <- as.Date(substr(i,18,25),format="%Y%m%d")
    if (data.img %within% season2) { 
      # variables ----
      sensore <- substr(i,1,2)
      satellite <- substr(i,3,4)
      level <- substr(i,6,9)
      path <- substr(i,11,13)
      row <- substr(i,14,16)
      year <- substr(i,18,21)
      month <- substr(i,22,23)
      day <- substr(i,24,25)
      cc <- substr(i,36,37)
      acquisitionTime <- substr(i,18,25)
      processingTime <- substr(i,27,34)
      tx <- substr(i,39,40)
      id <- substr(i,1,40)
      endname <- substr(i,10,40)
      
      setwd(file.path((L2TC.folder),i))
      if (file.exists(paste0("Snow&Cloud_mask_30m_",i,".tif"))==FALSE) {next}
      
      # Visualize data ----
      setwd(file.path(L2folder,i))
      
      QA <- rast(list.files(pattern=glob2rx(paste0("*QA_PIXEL*"))))
      
      # Visualize data ----
      
      
      if (satellite == "08" | satellite == "09") {
        rgb <- rast(c(list.files(pattern=glob2rx("*_B4.*")), 
                      list.files(pattern=glob2rx("*_B3.*")), 
                      list.files(pattern=glob2rx("*_B2.*")))) %>% 
          crop(aoi)
      } else {
        rgb <- rast(c(list.files(pattern=glob2rx("*_B3.*")), 
                      list.files(pattern=glob2rx("*_B2.*")), 
                      list.files(pattern=glob2rx("*_B1.*")))) %>% 
          crop(aoi)
      }
      
      plotRGB(rgb,r=1,g=2,b=3,axes=TRUE,mar = 2,main=i, scale=65535) #
      #rm(rgb)
      
      ### s1: distance to clouds/shadows & snow score: ----
      setwd(file.path((L2TC.folder),i))
      clouds <- rast(list.files(pattern=glob2rx(paste0("Snow&Cloud_mask_30m_",i,"*")))) %>% 
        app(fun = function(x) {return(x == 205)}) # retrieve cloud mask
      snow <- rast(list.files(pattern=glob2rx(paste0("Snow&Cloud_mask_30m_",i,"*")))) %>% 
        app(fun = function(x) {return(x == 100)}) # retrieve snow mask
      
      clouds.snow <- clouds+snow 
      clouds.snow[clouds.snow==2]<-1
      
      pixel_size = xres(clouds.snow)
      cdst_req <- 50
      
      background <- clouds.snow
      background[background == 0] <- NA
      dist_min <- distance(background,unit="m")/pixel_size # compute distance to clouds for each pixel
      
      s1 <- dist_min
      s1[s1 < cdst_req | s1 > 0] <- 1./(1.+exp(-0.01*(s1-cdst_req/2))) # distances lower than 10 pixels -> function
      s1[s1 >= cdst_req] <- 1   # distance greater than 10 pixels -> 1
      s1[s1 == 0] <- 0    # cloud and shadows -> 0
      
      rm(background, dist_min,clouds.snow)
      
      ### Create cloud & NA mask ----
      na.mask <- rgb[[3]]
      na.mask <- terra::resample(na.mask,clouds)
      na.mask[!is.na(na.mask)] <- 0
      rm(rgb)
      
      mask <- clouds
      mask[mask != 1] <- 0
      mask[mask == 1] <- NA
      
      mask <- mask+na.mask
      
      ### s2: sensor score ----
      if (year>2003 & sensore=="07"){ sv <- 0.5 } else { sv <- 1 }
      s2 <- rast(ncol=ncol(s1), nrow=nrow(s1), ext=ext(s1), crs = crs(s1)) #create an empty rast
      values(s2) <- sv    
      
      ### s3: doy score ----
      t.doy <- 285 # mid October
      sl <- 9.3
      doy <- as.numeric(strftime(data.img, format = "%j"))
      v <- (1/(sl*sqrt(2*3.14)))*exp(-0.5*((doy-t.doy)/sl)^2)
      v.max <- (1/(sl*sqrt(2*3.14)))*exp(-0.5*((t.doy-t.doy)/sl)^2)
      v.doy <- v/v.max
      
      s3 <- rast(ncol=ncol(s1), nrow=nrow(s1), ext=ext(s1), crs = crs(s1)) #create an empty rast
      values(s3) <- v.doy   
      
      ### s4: haze-Optimized Transformation (HOT) ----
      # Get bands 
      setwd(file.path(L2TC.folder,i))
      if (satellite == "04" | satellite == "05" | satellite == "07") {
        blue <- rast(list.files(pattern=glob2rx("*B1.*"))) 
        red <- rast(list.files(pattern=glob2rx("*B3.*")))
      } else if (satellite == "08" | satellite == "09") {
        blue <- rast(list.files(pattern=glob2rx('*B2.*')))
        red <- rast(list.files(pattern=glob2rx('*B4.*')))
      }
      
      HOT <- blue-0.5*red-0.08
      s4 <- 1/(1+exp((10/0.02)*(HOT+0.075)))
      
      ### s5: coverage ----
      cfc <- ncell(mask[mask==0])/ncell(mask)
      s5 <- rast(ncol=ncol(s1), nrow=nrow(s1), ext=ext(s1), crs = crs(s1)) #create an empty rast
      values(s5) <- cfc
      
      ### s6: snow free area ----
      sc <- ncell(snow[snow!=1])/ncell(snow)
      s6 <- rast(ncol=ncol(s1), nrow=nrow(s1), ext=ext(s1), crs = crs(s1)) #create an empty rast
      values(s6) <- sc
      
      ### Scores sum ----
      score.sum <- ((s1*s1.w)+(s2*s2.w)+(s3*s3.w)+#(s4*s4.w)+
                      (s5*s5.w)+(s6*s6.w)/(s1.w+s2.w+s3.w+s5.w+s6.w))  %>% #s4.w+
        mask(mask)
      #score.sum[is.na(score.sum)] <- 0
      
      rm(mask)
      # Save score img
      score.sum <- round(score.sum*1000,0)
      setwd(L2TC.folder)
      writeRaster(score.sum,paste0(i,"_SCORE.tif"),datatype="INT2S",filetype="GTiff",overwrite=TRUE)
      
      # Add to a list 
      #season2.imgs.list <- append(season2.imgs.list,i)
      season2.imgs.list[[length(season2.imgs.list) + 1]] = i
      
      rm(score.sum,s1,s2,s3,s4,s5)
    }
    
  }
  setwd(L2TC.folder)
  write.table(season2.imgs.list, file = paste0(y,"_season2.imgs.list.txt"), sep = ",")
  
  print(paste("year",y,"done"))
  gc()
}


#### Create BAP images ----

for (y in yy){
  
  # Growing season ----
  setwd(L2TC.folder)
  season1.imgs.list <- as.list(read.table(paste0(y,"_season1.imgs.list.txt"), sep = ","))
  season1.imgs <- c()
  blue.imgs <- c() 
  green.imgs <- c()
  red.imgs <- c()
  NIR.imgs <- c()
  SWIR1.imgs <- c()
  SWIR2.imgs <- c()
  
  for (i in 1:length(season1.imgs.list)){
    img <- season1.imgs.list[[i]]
    satellite <- substr(img,3,4)
    setwd(L2TC.folder)
    img.score <- rast(paste0(img,"_SCORE.tif"))
    season1.imgs <- c(season1.imgs, img.score)
    
    setwd(file.path(L2TC.folder,img))
    if (satellite == "04" | satellite == "05" | satellite == "07") {
      blue <- rast(list.files(pattern=glob2rx("*B1.*"))) %>% crop(aoi) %>% resample(img.score)
      green <- rast(list.files(pattern=glob2rx("*B2.*"))) %>% crop(aoi) %>% resample(img.score)
      red <- rast(list.files(pattern=glob2rx("*B3.*"))) %>% crop(aoi) %>% resample(img.score)
      NIR <- rast(list.files(pattern=glob2rx("*B4.*"))) %>% crop(aoi) %>% resample(img.score)
      SWIR1 <- rast(list.files(pattern=glob2rx("*B5.*"))) %>% crop(aoi) %>% resample(img.score)
      SWIR2 <- rast(list.files(pattern=glob2rx("*B7.*"))) %>% crop(aoi) %>% resample(img.score)
    } else if (satellite == "08" | satellite == "09") {
      blue <- rast(list.files(pattern=glob2rx('*B2.*'))) %>% crop(aoi) %>% resample(img.score)
      green <- rast(list.files(pattern=glob2rx('*B3.*'))) %>% crop(aoi) %>% resample(img.score)
      red <- rast(list.files(pattern=glob2rx('*B4.*'))) %>% crop(aoi) %>% resample(img.score)
      NIR <- rast(list.files(pattern=glob2rx('*B5.*'))) %>% crop(aoi) %>% resample(img.score)
      SWIR1 <- rast(list.files(pattern=glob2rx('*B6.*'))) %>% crop(aoi) %>% resample(img.score)
      SWIR2 <- rast(list.files(pattern=glob2rx('*B7.*'))) %>% crop(aoi) %>% resample(img.score)
    }
    
    blue.imgs <- c(blue.imgs,blue) 
    green.imgs <- c(green.imgs,green)
    red.imgs <- c(red.imgs,red)
    NIR.imgs <- c(NIR.imgs,NIR)
    SWIR1.imgs <- c(SWIR1.imgs,SWIR1)
    SWIR2.imgs <- c(SWIR2.imgs,SWIR2)
    
  }
  
  season1.imgs <- rast(season1.imgs)
  blue.imgs <- rast(blue.imgs) 
  green.imgs <- rast(green.imgs)
  red.imgs <- rast(red.imgs)
  NIR.imgs <- rast(NIR.imgs)
  SWIR1.imgs <- rast(SWIR1.imgs)
  SWIR2.imgs <- rast(SWIR2.imgs)
  
  score.matrix <- as.matrix(season1.imgs) %>% 
    as.data.frame() 
  columns <- unlist(season1.imgs.list)
  colnames(score.matrix) = columns
  score.matrix[is.na(score.matrix)] <- 0
  score.matrix <- score.matrix %>% 
    mutate_if(is.integer, as.numeric)
  score.matrix$BAP_Column <- apply(score.matrix,1,which.max)
  #score.matrix$BAP <- substr(as.character(score.matrix$BAP_Column),45,nchar(as.character(score.matrix$BAP_Column))-1)
  score.matrix$BAP <- as.numeric(score.matrix$BAP_Column)
  
  # Create score img
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- score.matrix$BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_SCORING.tif"),overwrite=TRUE,datatype="INT2S",filetype="GTiff")
  rm(img)
  
  
  # blue 
  blue.imgs.matrix <- as.matrix(blue.imgs) %>% 
    as.data.frame() %>% 
    cbind(score.matrix$BAP)
  blue.BAP <- cbind(c3 = blue.imgs.matrix[cbind(seq_along(blue.imgs.matrix$`score.matrix$BAP`), 
                                                blue.imgs.matrix$`score.matrix$BAP`)])
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- blue.BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_blue.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  # green
  green.imgs.matrix <- as.matrix(green.imgs) %>% 
    as.data.frame() %>% 
    cbind(score.matrix$BAP)
  green.BAP <- cbind(c3 = green.imgs.matrix[cbind(seq_along(green.imgs.matrix$`score.matrix$BAP`), green.imgs.matrix$`score.matrix$BAP`)])
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- green.BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_green.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  # red
  red.imgs.matrix <- as.matrix(red.imgs) %>% 
    as.data.frame() %>% 
    cbind(score.matrix$BAP)
  red.BAP <- cbind(c3 = red.imgs.matrix[cbind(seq_along(red.imgs.matrix$`score.matrix$BAP`), red.imgs.matrix$`score.matrix$BAP`)])
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- red.BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_red.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  # NIR
  NIR.imgs.matrix <- as.matrix(NIR.imgs) %>% 
    as.data.frame() %>% 
    cbind(score.matrix$BAP)
  NIR.BAP <- cbind(c3 = NIR.imgs.matrix[cbind(seq_along(NIR.imgs.matrix$`score.matrix$BAP`), NIR.imgs.matrix$`score.matrix$BAP`)])
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- NIR.BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_NIR.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  # SWIR1
  SWIR1.imgs.matrix <- as.matrix(SWIR1.imgs) %>% 
    as.data.frame() %>% 
    cbind(score.matrix$BAP)
  SWIR1.BAP <- cbind(c3 = SWIR1.imgs.matrix[cbind(seq_along(SWIR1.imgs.matrix$`score.matrix$BAP`), SWIR1.imgs.matrix$`score.matrix$BAP`)])
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- SWIR1.BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_SWIR1.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  # SWIR2
  SWIR2.imgs.matrix <- as.matrix(SWIR2.imgs) %>% 
    as.data.frame() %>% 
    cbind(score.matrix$BAP)
  SWIR2.BAP <- cbind(c3 = SWIR2.imgs.matrix[cbind(seq_along(SWIR2.imgs.matrix$`score.matrix$BAP`), SWIR2.imgs.matrix$`score.matrix$BAP`)])
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- SWIR2.BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_GS_195028_02_T1_SR_BAP_SWIR2.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  rm(season1.imgs,season1.imgs.list,blue.imgs,green.imgs,red.imgs,NIR.imgs,SWIR1.imgs,SWIR2.imgs)
  
  # Senescence season ----
  
  setwd(L2TC.folder)
  season2.imgs.list <- as.list(read.table(paste0(y,"_season2.imgs.list.txt"), sep = ","))
  season2.imgs <- c()
  blue.imgs <- c() 
  green.imgs <- c()
  red.imgs <- c()
  NIR.imgs <- c()
  SWIR1.imgs <- c()
  SWIR2.imgs <- c()
  
  for (i in 1:length(season2.imgs.list)){
    img <- season2.imgs.list[[i]]
    satellite <- substr(img,3,4)
    setwd(L2TC.folder)
    img.score <- rast(paste0(img,"_SCORE.tif"))
    season2.imgs <- c(season2.imgs, img.score)
    
    setwd(file.path(L2TC.folder,img))
    if (satellite == "04" | satellite == "05" | satellite == "07") {
      blue <- rast(list.files(pattern=glob2rx("*B1.*"))) %>% crop(aoi) %>% resample(img.score)
      green <- rast(list.files(pattern=glob2rx("*B2.*"))) %>% crop(aoi) %>% resample(img.score)
      red <- rast(list.files(pattern=glob2rx("*B3.*"))) %>% crop(aoi) %>% resample(img.score)
      NIR <- rast(list.files(pattern=glob2rx("*B4.*"))) %>% crop(aoi) %>% resample(img.score)
      SWIR1 <- rast(list.files(pattern=glob2rx("*B5.*"))) %>% crop(aoi) %>% resample(img.score)
      SWIR2 <- rast(list.files(pattern=glob2rx("*B7.*"))) %>% crop(aoi) %>% resample(img.score)
    } else if (satellite == "08" | satellite == "09") {
      blue <- rast(list.files(pattern=glob2rx('*B2.*'))) %>% crop(aoi) %>% resample(img.score)
      green <- rast(list.files(pattern=glob2rx('*B3.*'))) %>% crop(aoi) %>% resample(img.score)
      red <- rast(list.files(pattern=glob2rx('*B4.*'))) %>% crop(aoi) %>% resample(img.score)
      NIR <- rast(list.files(pattern=glob2rx('*B5.*'))) %>% crop(aoi) %>% resample(img.score)
      SWIR1 <- rast(list.files(pattern=glob2rx('*B6.*'))) %>% crop(aoi) %>% resample(img.score)
      SWIR2 <- rast(list.files(pattern=glob2rx('*B7.*'))) %>% crop(aoi) %>% resample(img.score)
    }
    
    blue.imgs <- c(blue.imgs,blue) 
    green.imgs <- c(green.imgs,green)
    red.imgs <- c(red.imgs,red)
    NIR.imgs <- c(NIR.imgs,NIR)
    SWIR1.imgs <- c(SWIR1.imgs,SWIR1)
    SWIR2.imgs <- c(SWIR2.imgs,SWIR2)
    
  }
  
  season2.imgs <- rast(season2.imgs)
  blue.imgs <- rast(blue.imgs) 
  green.imgs <- rast(green.imgs)
  red.imgs <- rast(red.imgs)
  NIR.imgs <- rast(NIR.imgs)
  SWIR1.imgs <- rast(SWIR1.imgs)
  SWIR2.imgs <- rast(SWIR2.imgs)
  
  score.matrix <- as.matrix(season2.imgs) %>% 
    as.data.frame() 
  columns <- unlist(season2.imgs.list)
  colnames(score.matrix) = columns
  score.matrix[is.na(score.matrix)] <- 0
  score.matrix <- score.matrix %>% 
    mutate_if(is.integer, as.numeric)
  score.matrix$BAP_Column <- apply(score.matrix,1,which.max)
  #score.matrix$BAP <- substr(as.character(score.matrix$BAP_Column),45,nchar(as.character(score.matrix$BAP_Column))-1)
  score.matrix$BAP <- as.numeric(score.matrix$BAP_Column)
  
  # Create score img
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- score.matrix$BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_SCORING.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  # blue 
  blue.imgs.matrix <- as.matrix(blue.imgs) %>% 
    as.data.frame() %>% 
    cbind(score.matrix$BAP)
  blue.BAP <- cbind(c3 = blue.imgs.matrix[cbind(seq_along(blue.imgs.matrix$`score.matrix$BAP`), blue.imgs.matrix$`score.matrix$BAP`)])
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- blue.BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_blue.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  # green
  green.imgs.matrix <- as.matrix(green.imgs) %>% 
    as.data.frame() %>% 
    cbind(score.matrix$BAP)
  green.BAP <- cbind(c3 = green.imgs.matrix[cbind(seq_along(green.imgs.matrix$`score.matrix$BAP`), green.imgs.matrix$`score.matrix$BAP`)])
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- green.BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_green.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  # red
  red.imgs.matrix <- as.matrix(red.imgs) %>% 
    as.data.frame() %>% 
    cbind(score.matrix$BAP)
  red.BAP <- cbind(c3 = red.imgs.matrix[cbind(seq_along(red.imgs.matrix$`score.matrix$BAP`), red.imgs.matrix$`score.matrix$BAP`)])
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- red.BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_red.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  # NIR
  NIR.imgs.matrix <- as.matrix(NIR.imgs) %>% 
    as.data.frame() %>% 
    cbind(score.matrix$BAP)
  NIR.BAP <- cbind(c3 = NIR.imgs.matrix[cbind(seq_along(NIR.imgs.matrix$`score.matrix$BAP`), NIR.imgs.matrix$`score.matrix$BAP`)])
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- NIR.BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_NIR.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  # SWIR1
  SWIR1.imgs.matrix <- as.matrix(SWIR1.imgs) %>% 
    as.data.frame() %>% 
    cbind(score.matrix$BAP)
  SWIR1.BAP <- cbind(c3 = SWIR1.imgs.matrix[cbind(seq_along(SWIR1.imgs.matrix$`score.matrix$BAP`), SWIR1.imgs.matrix$`score.matrix$BAP`)])
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- SWIR1.BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_SWIR1.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  # SWIR2
  SWIR2.imgs.matrix <- as.matrix(SWIR2.imgs) %>% 
    as.data.frame() %>% 
    cbind(score.matrix$BAP)
  SWIR2.BAP <- cbind(c3 = SWIR2.imgs.matrix[cbind(seq_along(SWIR2.imgs.matrix$`score.matrix$BAP`), SWIR2.imgs.matrix$`score.matrix$BAP`)])
  img <- rast(ncol=ncol(SWIR1), nrow=nrow(SWIR1), ext=ext(SWIR1), crs = crs(SWIR1))  #create an empty raster
  values(img) <- SWIR2.BAP
  setwd(bap.folder)
  writeRaster(img,filename=paste0(y,"_SS_195028_02_T1_SR_BAP_SWIR2.tif"),overwrite=TRUE,datatype="INT4S",filetype="GTiff")
  rm(img)
  
  rm(season2.imgs,season2.imgs.list,blue.imgs,green.imgs,red.imgs,NIR.imgs,SWIR1.imgs,SWIR2.imgs)
  
  setwd(L2TC.folder)
  print(paste("year",y,"done"))
}

#### Inspect results ----

setwd(bap.folder)
for (y in yy) {
  
  tci_gs <- rast(c(paste0(y,"_GS_195028_02_T1_SR_BAP_red.tif"),
                   paste0(y,"_GS_195028_02_T1_SR_BAP_green.tif"),
                   paste0(y,"_GS_195028_02_T1_SR_BAP_blue.tif")))
  
  tci_ss <- rast(c(paste0(y,"_SS_195028_02_T1_SR_BAP_red.tif"),
                   paste0(y,"_SS_195028_02_T1_SR_BAP_green.tif"),
                   paste0(y,"_SS_195028_02_T1_SR_BAP_blue.tif")))
  
  png(paste0(y,"_TCI.png"), width=10, height=5, units="in", res=300)
  par(mfrow=c(1,2),adj=0.3)
  plotRGB(tci_gs, r = 1, g = 2, b = 3, stretch = "lin", mar=1, main = "Growing season")
  plotRGB(tci_ss, r = 1, g = 2, b = 3, stretch = "lin", mar=1, main = "Senescence season")
  mtext(y, side = 3, line = - 2, outer = TRUE)
  dev.off() #only 129kb in size
  
  rm(tci_gs,tci_ss)
  print(paste("year",y,"done"))
  
}
