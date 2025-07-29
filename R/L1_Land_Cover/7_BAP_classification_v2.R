
rm(list = ls())

# Packages ----
list.of.packages <- c(
  "devtools",
  "terra",
  #"rgdal",
  "data.table",
  "ff",
  "lubridate",
  "dplyr",
  "ggplot2",
  "ggspatial",
  "tidyterra",
  "RStoolbox",
  "grid",
  "parallel",
  "ranger",
  "smoothr"
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
L2folder <- file.path(root,"Landsat/level2/AC")
ancillary.folder <- file.path(root,"Ancillary/5m")
rfFolder <- file.path(root,"Landsat/RF")

numCores <- detectCores()-2

#### Retrieve available images and scores  ----
imgs.list <- list.files(path = L2folder,pattern="T1")
f <- function(x){as.Date(substr(x,18,25),format="%Y%m%d")}
days.list <- lapply(imgs.list,FUN = f)
days.list <- sort(as.Date(unlist(days.list)))

first.date <- min(days.list)
last.date <- max(days.list)

yy <- (strftime(first.date, format = "%Y"):strftime(last.date, format = "%Y"))

aoi <- vect("D:/ABRESO/PNGP/Shapefile/PNGP.shp")


sel.pred <- c("L3","GS_BSI","GS_CI","GS_GLI","GS_mARI","GS_NBR2","GS_NDBI","GS_NDWI",
              "GS_SBI","SS_ARVI","SS_BSI","SS_CI","SS_CVI","GS_EVI","SS_EVI","SS_GLI",
              "SS_GVI","SS_GVMI","SS_mARI","SS_MCARI1","SS_NBR2","SS_NDBI",
              "SS_WET","GS_NIR","GS_SWIR1","SS_blue","SS_NIR","SS_SWIR1",
              "elev","roughness","slope","twi","Geologia")

sel.pred.names <- c("LULC","GS_BSI","GS_CI","GS_GLI","GS_mARI","GS_NBR2","GS_NDBI","GS_NDWI",
                    "GS_SBI","SS_ARVI","SS_BSI","SS_CI","SS_CVI","GS_EVI","SS_EVI","SS_GLI",
                    "SS_GVI","SS_GVMI","SS_mARI","SS_MCARI1","SS_NBR2","SS_NDBI",
                    "SS_WET","GS_NIR","GS_SWIR1","SS_blue","SS_NIR","SS_SWIR1",
                    "Elevation","Roughness","Slope","TWI","Geology")

# FUNCTIONS
# Function to train Random Forest models for the ensemble ----
train_rf_ensemble <- function(train_data, n_models = 5) {
  models <- list()
  set.seed(42)  # For reproducibility
  
  for (i in 1:n_models) {
    # Create a bootstrap sample of the training data
    boot_data <- train_data[sample(1:nrow(train_data), replace = TRUE), ]
    
    # Train the Random Forest model
    rf_model <- ranger(
      formula = L3 ~ .,  # Adjust formula to match your dataset
      data = boot_data,
      probability = TRUE,
      num.trees = 100,
      mtry = 16,
      sample.fraction = 0.51,
      min.node.size = 4
    )
    
    # Store the trained model
    models[[i]] <- rf_model
  }
  
  return(models)
}

# Custom predict function to use each ranger model with terra::predict() ----
rf_predict_function <- function(model, data) {
  # Convert the data (a matrix) into a data frame for the model
  data <- as.data.frame(data)
  
  # Predict class probabilities using the Random Forest model
  pred <- terra::predict(model, data = data, type = "response")$predictions
  
  # Return the probabilities (or class labels if needed)
  return(pred)
}



###################### Grow classifiers #####################################

for (y in yy) {
  
  # Combine all the pixels selectd for each class to build the training dataset 
  setwd("D:/ABRESO/PNGP/Landsat/z_statistics/Training")
  trainings <- rast(list.files(pattern=glob2rx(paste0(y,"*training*30m*8cl*")))) %>% mask(aoi) 
  all.trainings <- app(trainings, fun=sum, na.rm=TRUE) # Get training variables
  rm(trainings)
  
  #### Import spectral predictors ####
  setwd(l3.folder)
  indlist <- list.files(getwd(), pattern=glob2rx(paste0(y,'*.tif$')), all.files=TRUE, full.names=TRUE)
  ind.names <- lapply(indlist,function(x) {paste0(substr(x,43,44),"_",substr(x,66,nchar(x)-4))})
  indices<- rast(indlist) 
  names(indices) <- unlist(ind.names)
  
  setwd(bap.folder)
  bandslist <- list.files(getwd(), pattern=glob2rx(paste0(y,'*.tif$')), all.files=TRUE, full.names=TRUE)
  bands.names <- lapply(bandslist,function(x) {paste0(substr(x,43,44),"_",substr(x,66,nchar(x)-11))})
  bands<- rast(bandslist) 
  
  # Import static predictors with geology
  setwd("D:/ABRESO/PNGP/Ancillary/5m")
  predlist <- list.files(getwd(), pattern=glob2rx(paste0('*.tif$')), all.files=TRUE, full.names=TRUE)
  geo <- rast(predlist[6]) %>% terra::resample(indices, method="near")  
  predlist <- predlist[-6]
  pred.names <- lapply(predlist,function(x) {substr(x,29,nchar(x)-4)})
  predictors <- rast(predlist) %>% terra::resample(indices) 
  
  df.training <- c(all.trainings,geo,indices,bands,predictors) %>% mask(aoi) 
  all.names <- c("L3","Geologia",ind.names,bands.names,pred.names)
  names(df.training) <- unlist(all.names)
  
  df.train <- as.data.frame(df.training) %>% 
    na.omit()
  
  rm(df.training)
  gc()
  
  df.train$L3 <- as.factor(df.train$L3)
  df.train$Geologia <- as.factor(df.train$Geologia)
  
  # mantain only selected predictors and discard the others
  df.train <- select(df.train, all_of(sel.pred))
  dim(df.train)
  t <- df.train %>% group_by(L3) %>% slice_sample(prop = 0.25) %>% tally() %>% 
    as.data.frame()
  print(t)
  
  # Subset 59%
  #df <- df.train %>% group_by(L3) %>% slice_sample(prop = 0.59) %>% 
    #as.data.frame()
  #tt <- df %>% group_by(L3) %>%tally() %>% 
    #as.data.frame()
  #print(tt)
  
  df.train <- df.train %>% group_by(L3) %>% slice_sample(prop = 0.25) %>% 
    as.data.frame()
  
  #fwrite(df.train, file = paste0("D:/ABRESO/PNGP/LULC/",y,"_training_dataset_25proportional_30m_8cl.csv"), sep = ",")
  
  # Train 5 Random Forest models for the ensemble
  rf_models <- train_rf_ensemble(train_data = df.train, n_models = 5)
  
  #### Import spectral predictors ####
  setwd(l3.folder)
  indlist <- list.files(getwd(), pattern=glob2rx(paste0(y,'*.tif$')), all.files=TRUE, full.names=TRUE)
  ind.names <- lapply(indlist,function(x) {paste0(substr(x,43,44),"_",substr(x,66,nchar(x)-4))})
  indices<- rast(indlist) 
  names(indices) <- unlist(ind.names)
  
  setwd(bap.folder)
  bandslist <- list.files(getwd(), pattern=glob2rx(paste0(y,'*.tif$')), all.files=TRUE, full.names=TRUE)
  bands.names <- lapply(bandslist,function(x) {paste0(substr(x,43,44),"_",substr(x,66,nchar(x)-11))})
  bands<- rast(bandslist) 
  
  # Import static predictors with geology
  setwd("D:/ABRESO/PNGP/Ancillary/5m")
  predlist <- list.files(getwd(), pattern=glob2rx(paste0('*.tif$')), all.files=TRUE, full.names=TRUE)
  geo <- rast(predlist[6]) %>% terra::resample(indices, method="near")  
  predlist <- predlist[-6]
  pred.names <- lapply(predlist,function(x) {substr(x,29,nchar(x)-4)})
  predictors <- rast(predlist) %>% terra::resample(indices) 
  
  #load(paste0(rfFolder,"/",y,"_PNGP_LULC_classificator_8cl_SG.Rda"))
  
  stack <- c(geo,indices,bands,predictors)
  all.names <- c("Geologia",ind.names,bands.names,pred.names)
  names(stack) <- unlist(all.names)
  pngp.stack <- subset(stack, sel.pred[2:length(sel.pred)]) %>% mask(aoi)
  pngp.stack$Geologia <- as.factor(pngp.stack$Geologia)
  
  
  rf_predict_function <- function(model, data) {
    data_df <- as.data.frame(data)  # Convert raster matrix to a data frame
    pred <- predict(model, data = data_df, type = "response")
    return(pred$predictions)
  }
  
  # Apply terra::predict() to each model ----
  
  all_predictions <- rast(pngp.stack[[1]])
  for (i in 1:length(rf_models)) {
    # Create a wrapper for the current model
    model <- rf_models[[i]]
    
    # Predict using terra::predict with custom function
    prediction <- terra::predict(pngp.stack, model, fun=predict,const=NULL, 
                                 na.rm=TRUE,   cpkgs="ranger") #index=c,
    
    # Store the prediction raster
    all_predictions <- c(all_predictions,prediction)
  }
  
  # Aggregate the predictions from all models (averaging) ----
  for (c in 1:8) {
    layer_indices <- which(names(all_predictions) == paste0("X",c))
    
    # Subset by those indices
    subset_by_name <- all_predictions[[layer_indices]]
    
    # Calculate the average across all ensemble models
    ensemble_predictions <- mean(subset_by_name)
    
    setwd("D:/ABRESO/PNGP/Landsat/RF/ensemble")
    # Save the averaged probability map
    writeRaster(ensemble_predictions, paste0(y,"_ensemble_class_",c,"_probability.tif"), overwrite = TRUE)
  }
  
  
  setwd("D:/ABRESO/PNGP/Landsat/RF/ensemble")
  ensemble_predictions <- rast(list.files(pattern=glob2rx(paste0(y,"_ensemble_class_*_probability.tif$"))))
  
  # Final classification based on the highest probability ----
  final_class <- which.max(ensemble_predictions)
  
  # Save the final classification map
  writeRaster(final_class, paste0(y,"_PNGP_LULC_ensemble.tif"), overwrite = TRUE)
  
}
