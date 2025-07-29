################################ Libraries ################################
library(terra)
library(dplyr)
library(caret)
library(ggplot2)  
library(tidyr)

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
cond <- sapply(yy, function(x) { substr(tail(x,1),nchar(tail(x,1)),nchar(tail(x,1))) == 0 | 
    substr(tail(x,1),nchar(tail(x,1)),nchar(tail(x,1))) == 5})
yy <- yy[cond]

# Import data 
lulc <- rast("D:/ABRESO/PNGP/LULC/Training/LULC_L3_5m.tif")
setwd("D:/ABRESO/PNGP/Landsat/Validation/ensemble/postprocessing")

# Function to calculate the total area for each land cover class
calculate_landcover_area <- function(r) {
  
  # Get the unique classes (values) in the raster
  classes <- unique(values(r, na.rm=TRUE))
  
  # Calculate the area of each cell in square meters (adjust units based on your CRS)
  cell_area <- cellSize(r, unit="m")
  
  # Calculate the frequency (count) of each class in the raster
  freq_table <- freq(r)
  
  # Calculate the total area for each class
  total_area_by_class <- freq_table
  total_area_by_class$area_m2 <- total_area_by_class$count * as.numeric(global(cell_area, fun=mean))
  total_area_by_class$area_ha <- round(total_area_by_class$area_m2 / 10000,2)
  
  return(total_area_by_class)
}

# Initialize a dataframe to store validation metrics
validation_results <- data.frame()
validation_results_2 <- data.frame()
validation_results_3 <- data.frame()

yy <- c(1990,2000,2010,2020)
setwd("D:/ABRESO/PNGP/Landsat/RF/ensemble/postprocessing_NEW/validation")

for (y in yy) {
  
  ############################## Validation ################################
  #### Import ground references ####
  gt <- vect(list.files(path="D:/ABRESO/PNGP/LULC/Validation", pattern=glob2rx(paste0("val",y,".shp$")),full.names=TRUE)) %>% 
    rasterize(lulc, field="LULC",background=NA)
  
  pixel.area <- xres(gt)*yres(gt)
  
  # Get statistics on ground truth
  #tt <- gt %>% as.data.frame() %>% 
    #group_by(LULC) %>% 
    #tally() %>% 
    #as.data.frame()
  #tt$area_ha <- round((tt$n*pixel.area)/10000,0)
  #tt <- rbind(tt, colSums(tt))
  
  ### Import classified 
  #pred <- rast(list.files(path="D:/ABRESO/PNGP/Landsat/RF/ensemble/postprocessing", pattern=glob2rx(paste0(y,"_PNGP_LULC*.tif$")),full.names=TRUE)) %>% 
    #terra::resample(gt, method="near")
  if (y==1990){
    pred <- rast(list.files(path=file.path(root,"Landsat/RF/ensemble/postprocessing_NEW"), 
                            pattern=glob2rx(paste0("LC_1985_1995_LC_reclassified.tif$")),full.names=TRUE)) %>% 
      terra::resample(gt, method="near")
  } else if (y == 2000){
    pred <- rast(list.files(path=file.path(root,"Landsat/RF/ensemble/postprocessing_NEW"), 
                           pattern=glob2rx(paste0("LC_1995_2005_LC_reclassified.tif$")),full.names=TRUE)) %>% 
      terra::resample(gt, method="near")
  } else if (y == 2000){
    pred <- rast(list.files(path=file.path(root,"Landsat/RF/ensemble/postprocessing_NEW"), 
                           pattern=glob2rx(paste0("LC_2005_2015_LC_reclassified.tif$")),full.names=TRUE)) %>% 
      terra::resample(gt, method="near")
  } else {
    pred <- rast(list.files(path=file.path(root,"Landsat/RF/ensemble/postprocessing_NEW"), 
                            pattern=glob2rx(paste0("LC_2015_2023_LC_reclassified.tif$")),full.names=TRUE)) %>% 
      terra::resample(gt, method="near")
  }
  
  names(pred) <- "prediction"
  f <- function(x){x<=8}
  mask <- app(pred, fun=f)
  mask[mask==0] <- NA
  pred <- mask(pred,mask)
  
  df <- c(gt,pred) %>% 
    as.data.frame() %>% 
    na.omit() %>% 
    mutate_if(is.numeric, as.factor)
  
  cm <- confusionMatrix(df$prediction, df$LULC)
  
  write.csv(as.matrix(cm,what="overall"),file=paste0("overall_scores_",y,".csv")) # save overall metrics
  write.csv(as.matrix(cm, what = "classes"),file=paste0("class_scores_",y,".csv")) # save class metrics
  write.csv(as.matrix(cm),file=paste0("cm_table_",y,".csv")) # save table
  capture.output(cm, file = paste0("confusionMatrix_results_",y,".txt"))
  
  ##########################################################################
  
  # Extract values from both rasters, removing NA values
  #vals <- c(pred, gt)
  #vals <- na.omit(as.data.frame(vals))
  
  # Overall accuracy
  overall_accuracy <- cm$overall['Accuracy']
  
  # Kappa coefficient
  kappa <- cm$overall['Kappa']
  
  # User's accuracy (Precision), Producer's accuracy (Recall), F1-score
  user_accuracy <- cm$byClass[, 'Pos Pred Value']  # User's accuracy (Precision)
  producer_accuracy <- cm$byClass[, 'Sensitivity'] # Producer's accuracy (Recall)
  f1_score <- 2 * (user_accuracy * producer_accuracy) / (user_accuracy + producer_accuracy)
  
  # Area-based metrics
  # Total mapped area for each class
  mapped_area <- calculate_landcover_area(pred)
  
  # Weighted Accuracy
  class_frequency <- table(df$LULC)  # Ground truth class frequencies
  weighted_accuracy <- sum((class_frequency / sum(class_frequency)) * user_accuracy, na.rm=TRUE)
  
  # G-Mean
  g_mean <- sqrt(prod(producer_accuracy, na.rm=TRUE))
  
  # Intersection over Union (IoU) for each class
  iou <- cm$byClass[, 'Sensitivity'] / (cm$byClass[, 'Sensitivity'] + cm$byClass[, 'Pos Pred Value'] + cm$byClass[, 'Neg Pred Value'] - cm$byClass[, 'Sensitivity'])
  
  # True Skill Statistic (TSS)
  specificity <- cm$byClass[, 'Specificity']
  tss <- producer_accuracy + specificity - 1
  
  
  # Compile metrics
  metrics <- list(
    overall_accuracy = overall_accuracy,
    kappa = kappa,
    user_accuracy = user_accuracy,
    producer_accuracy = producer_accuracy,
    f1_score = f1_score,
    weighted_accuracy = weighted_accuracy,
    g_mean = g_mean,
    iou = iou,
    tss = tss,
    cmrix = cm$table
  )
  
  output_file <- paste0("validation_results_",y)
  # Save the confusion matrix
  write.csv(as.data.frame(metrics$cmrix), file = paste0(output_file, "_cmrix.csv"))
  
  # Save the overall accuracy, kappa, and MAE
  write.csv(data.frame(Overall_Accuracy = metrics$overall_accuracy, Kappa = metrics$kappa),
            file = paste0(output_file, "_overall_kappa.csv"))
  
  # Save class-wise metrics (User's accuracy, Producer's accuracy, F1-score, IoU, TSS)
  class_metrics <- data.frame(
    User_Accuracy = metrics$user_accuracy,
    Producer_Accuracy = metrics$producer_accuracy,
    F1_Score = metrics$f1_score,
    IoU = metrics$iou,
    TSS = metrics$tss
  )
  write.csv(class_metrics, file = paste0(output_file, "_class_metrics.csv"))
  
  # Save area-based metrics (Mapped area, spatial autocorrelation)
  write.csv(mapped_area, file = paste0(output_file, "_area_metrics.csv"))
  
  # Save additional metrics (Weighted Accuracy, G-Mean)
  additional_metrics <- data.frame(
    Weighted_Accuracy = metrics$weighted_accuracy,
    G_Mean = metrics$g_mean
  )
  
  #output_file <- paste0("validation_results_",y)
  write.csv(additional_metrics, file = paste0(output_file, "_additional_metrics.csv"))
  
  print(paste("Metrics saved to:", output_file))
  
  # Store results in the validation_results dataframe
  validation_results_2 <- rbind(validation_results_2, 
                              data.frame(Year = y,
                                         LULC = rownames(cm$table),
                                         PA = producer_accuracy,
                                         #PA_CI = result[, "PA±"],
                                         UA = user_accuracy,
                                         #UA_CI = result[, "UA±"],
                                         F1 =  f1_score
                              ))
  
  validation_results_3 <- rbind(validation_results_3, 
                                data.frame(Year = y,
                                           OA = overall_accuracy,
                                           K = kappa,
                                           WA = weighted_accuracy
                                ))
  
  
  ################### Adjusted area estimation #############################
  tb <- as.data.frame(pred) %>%
    group_by(prediction) %>% 
    summarise(pixels = n()) %>% 
    mutate(`area m^2` = pixels * pixel.area, #multiply by pixel area
           `area ha` = round(`area m^2`/10000,2)) 
  
  #print(tb)
  ###
  
  # import classification image and validation shapefile
  img.class <- pred
  shp.valid <- vect(list.files(path="D:/ABRESO/PNGP/LULC/Validation", pattern=glob2rx(paste0("val",y,".shp$")),full.names=TRUE))
  
  # create regular accuracy matrix 
  confmat <- as.matrix(cm)
  # get number of pixels per class and convert in ha
  #imgVal <- as.factor(values(img.class))
  nclass <- length(unique(shp.valid$LULC))
  #maparea <- sapply(1:nclass, function(x) sum(imgVal == x))
  #maparea <- maparea * res(img.class)[1] ^ 2 / 10000
  maparea <- tb$`area ha`
  
  # set confidence interval
  conf <- 1.96
  
  # total  map area
  A <- sum(maparea)
  # proportion of area mapped as class i
  W_i <- maparea / A
  # number of reference points per class
  n_i <- rowSums(confmat) 
  # population error matrix (Eq.4)
  p <- W_i * confmat / n_i
  p[is.na(p)] <- 0
  
  # area estimation
  p_area <- colSums(p) * A
  # area estimation confidence interval (Eq.10)
  p_area_CI <- conf * A * sqrt(colSums((W_i * p - p ^ 2) / (n_i - 1))) 
  
  # overall accuracy (Eq.1)
  OA <- sum(diag(p))
  # producers accuracy (Eq.2)
  PA <- diag(p) / colSums(p)
  # users accuracy (Eq.3)
  UA <- diag(p) / rowSums(p)
  # F1 score
  F1 <- 2 * (PA * UA) / (PA + UA)
  
  # overall accuracy confidence interval (Eq.5)
  OA_CI <- conf * sqrt(sum(W_i ^ 2 * UA * (1 - UA) / (n_i - 1)))
  # user accuracy confidence interval (Eq.6)
  UA_CI <- conf * sqrt(UA * (1 - UA) / (n_i - 1)) 
  # producer accuracy confidence interval (Eq.7)
  N_j <- sapply(1:nclass, function(x) sum(maparea / n_i * confmat[ , x]) )
  tmp <- sapply(1:nclass, function(x) sum(maparea[-x] ^ 2 * confmat[-x, x] / n_i[-x] * ( 1 - confmat[-x, x] / n_i[-x]) / (n_i[-x] - 1)) )
  PA_CI <- conf * sqrt(1 / N_j ^ 2 * (maparea ^ 2 * ( 1 - PA ) ^ 2 * UA * (1 - UA) / (n_i - 1) + PA ^ 2 * tmp))
  
  # gather results
  result <- matrix(c(p_area, p_area_CI, PA * 100, PA_CI * 100, UA * 100, UA_CI * 100, F1 * 100, c(OA * 100, rep(NA, nclass-1)), c(OA_CI * 100, rep(NA, nclass-1))), nrow = nclass)
  result <- round(result, digits = 2) 
  rownames(result) <- levels(as.factor(shp.valid$LULC))
  colnames(result) <- c("ha", "ha±", "PA", "PA±", "UA", "UA±", "F1","OA", "OA±")
  class(result) <- "table"
  #result
  
  write.csv(as.matrix(result),file=paste0("Adjusted_area_",y,".csv")) # save table
  capture.output(result, file = paste0("Adjusted_area_results_",y,".txt"))
  
  # Store results in the validation_results dataframe
  validation_results <- rbind(validation_results, 
                              data.frame(Year = y,
                                         LULC = rownames(result),
                                         PA = result[, "PA"],
                                         #PA_CI = result[, "PA±"],
                                         UA = result[, "UA"],
                                         #UA_CI = result[, "UA±"],
                                         F1 = result[, "F1"],
                                         OA = result[, "OA"][1]#,
                                         #OA_CI = result[, "OA±"][1]
                                         ))
  
  ##############################################################################
  
  rm(cm,df,pred,gt,tb,result,A,W_i,n_i,p,p_area,p_area_CI,confmat,
     OA,PA,UA,OA_CI,UA_CI,PA_CI,N_j,tmp,shp.valid,img.class)
  print(paste("Year",y,"done"))
}

# Visualization of validation results
validation_results <- validation_results %>% 
  gather(key = "Metric", value = "Value", PA,  UA,  F1, OA) # PA_CI,UA_CI,, OA_CI

validation_results_2 <- validation_results_2 %>% 
  gather(key = "Metric", value = "Value", PA,  UA,  F1)

validation_results_3 <- validation_results_3 %>% 
  gather(key = "Metric", value = "Value", OA,  K,  WA) 

# Custom labels and colors
c.labels <- c("Rocks, screes, debris",
              "Snow and glaciers",
              "Water",
              "Broadleaved",
              "Coniferous",
              "Grassland",
              "Shrubs",
              "Wetland")

#c.values <- c("#bdbdbd", "#ccece6", "#084594","#319812", "#024618", "#97e800", "#afbea2","#b98753")

c.values <- c(
  "#969696",  # grey - più scuro per distinguersi meglio
  "#99d8c9",  # teal chiaro - più saturo
  "#084594",  # blu-verdastro molto scuro (ex blu scuro, ma distinguibile)
  "#31a354",  # verde più saturo ma distinguibile
  "#006d00",  # verde scuro intenso (ex #024618)
  "#c7ea46",  # verde chiaro/giallastro più acceso
  "#8c9e7a",  # verde-grigio più contrastato
  "#b07d45"   # marrone più saturo
)


p1 <- ggplot(validation_results, aes(x = Year, y = Value, color = LULC, group = LULC)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_minimal() +
  labs( x = "Year", y = "Value", color = "LULC Class") + #title = "Validation Metrics Over Years",
  scale_color_manual(values = c.values, labels=c.labels)


p2 <- ggplot(validation_results_2, aes(x = Year, y = Value, color = LULC, group = LULC)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Metric, scales = "free_y") +
  ylim(c(0,1)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs( x = "Year", y = "Value", color = " ") + #title = "Validation Metrics Over Years",
  scale_color_manual(values = c.values, labels=c.labels)

p3 <- ggplot(validation_results_3, aes(x = Year, y = Value, color = Metric, group = Metric)) +
  geom_line() +           
  geom_point() +               
  labs(x = "Year", 
       y = "Value", 
       color = "Metric") + 
  ylim(min(validation_results_3$Value), 1) +
  theme_minimal() 

