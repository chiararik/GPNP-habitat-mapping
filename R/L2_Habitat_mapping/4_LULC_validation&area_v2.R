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
lulc <- rast("D:/ABRESO/PNGP/LULC/Training/LULC_L4_5m.tif")
setwd("D:/ABRESO/PNGP/Landsat/RF/ensemble/habitat_mapping/postprocessing")

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


############################## Validation ################################
#### Import ground references ####
#mask <- rast(list.files("D:/ABRESO/PNGP/Landsat/RF/ensemble/postprocessing_NEW", 
                        #pattern=glob2rx(paste0('2009_LC_reclassified.tif$')), all.files=TRUE, full.names=TRUE))

### Import classified 
#pred <- rast(list.files(path="D:/ABRESO/PNGP/Landsat/RF/ensemble/habitat_mapping/postprocessing", 
                        #pattern=glob2rx("2009_habitat_map_reclassified.tif$"),full.names=TRUE)) 
pred <- rast("D:/ABRESO/PNGP/Landsat/RF/ensemble/habitat_mapping/postprocessing/5Y/LC_2005_2015_LC_reclassified.tif")

gt <- rast("D:/ABRESO/PNGP/LULC/Training/LULC_L4_5m.tif") %>% terra::resample(pred, method="near")

gt <- gt %>% terra::resample(pred, method="near")

gt[gt==0]<- NA
gt[gt==11|gt==12|gt==13|gt==14] <- 1
gt[gt==21|gt==23]<-2
gt[gt==31|gt==35]<-3
gt[gt==81|gt==82|gt==83|gt==84|gt==85]<-8

pixel.area <- xres(gt)*yres(gt)

# Get statistics on ground truth
#tt <- gt %>% as.data.frame() %>% 
#group_by(LULC) %>% 
#tally() %>% 
#as.data.frame()
#tt$area_ha <- round((tt$n*pixel.area)/10000,0)
#tt <- rbind(tt, colSums(tt))

names(gt) <- "LULC"
names(pred) <- "prediction"
#f <- function(x){x<=8}
#mask <- app(pred, fun=f)
#mask[mask==0] <- NA
#pred <- mask(pred,mask)

df <- c(gt,pred) %>% 
  as.data.frame() %>% 
  na.omit() %>% 
  mutate_if(is.numeric, as.factor)

cm <- confusionMatrix(df$prediction, df$LULC)

setwd("validation")
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
write.csv(as.data.frame(metrics$cmrix), file = paste0(output_file,"_",y, "_cmrix.csv"))

# Save the overall accuracy, kappa, and MAE
write.csv(data.frame(Overall_Accuracy = metrics$overall_accuracy, Kappa = metrics$kappa),
          file = paste0(output_file,"_",y, "_overall_kappa.csv"))

# Save class-wise metrics (User's accuracy, Producer's accuracy, F1-score, IoU, TSS)
class_metrics <- data.frame(
  User_Accuracy = metrics$user_accuracy,
  Producer_Accuracy = metrics$producer_accuracy,
  F1_Score = metrics$f1_score,
  IoU = metrics$iou,
  TSS = metrics$tss
)
write.csv(class_metrics, file = paste0(output_file, "_",y,"_class_metrics.csv"))

# Save area-based metrics (Mapped area, spatial autocorrelation)
write.csv(mapped_area, file = paste0(output_file,"_",y, "_area_metrics.csv"))

# Save additional metrics (Weighted Accuracy, G-Mean)
additional_metrics <- data.frame(
  Weighted_Accuracy = metrics$weighted_accuracy,
  G_Mean = metrics$g_mean
)

#output_file <- paste0("validation_results_",y)
write.csv(additional_metrics, file = paste0(output_file,"_",y, "_additional_metrics.csv"))

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

##############################################################################

# Visualization of validation results ----
validation_results <- validation_results %>% 
  gather(key = "Metric", value = "Value", PA,  UA,  F1, OA) # PA_CI,UA_CI,, OA_CI

validation_results_2 <- validation_results_2 %>% 
  gather(key = "Metric", value = "Value", PA,  UA,  F1)

validation_results_3 <- validation_results_3 %>% 
  gather(key = "Metric", value = "Value", OA,  K,  WA) 

p1 <- ggplot(validation_results, aes(x = Year, y = Value, group = LULC)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_minimal() +
  labs( x = "Year", y = "Value", color = "LULC Class") + #title = "Validation Metrics Over Years",
  scale_color_manual(values = c.values, labels=c.labels)


p2 <- ggplot(validation_results_2, aes(x = Year, y = Value,  group = LULC)) + #color = LULC,
  geom_line() +
  geom_point() +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_minimal() +
  labs( x = "Year", y = "Value", color = "LULC Class") #+ #title = "Validation Metrics Over Years",
  #scale_color_manual(values = c.values, labels=c.labels)

p3 <- ggplot(validation_results_3, aes(x = Year, y = Value, color = Metric, group = Metric)) +
  geom_line() +           
  geom_point() +               
  labs(x = "Year", 
       y = "Value", 
       color = "Metric") + 
  ylim(min(validation_results_3$Value), 1) +
  theme_minimal() 

