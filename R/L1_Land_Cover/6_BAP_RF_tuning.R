
rm(list = ls())

# Packages ----
list.of.packages <- c(
  "devtools",
  "terra",
  "data.table",
  "ff",
  "dplyr",
  "ggplot2",
  "ggspatial",
  "tidyterra",
  "parallel",
  "ranger",
  "tuneRanger",
  "mlr",
  "spm",
  "caret"
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

y <- 2009

####################################################################

sel.pred <- c("L3","GS_BSI","GS_CI","GS_GLI","GS_mARI","GS_NBR2","GS_NDBI","GS_NDWI",
              "GS_SBI","SS_ARVI","SS_BSI","SS_CI","SS_CVI","SS_EVI","SS_GLI",
              "SS_GVI","SS_GVMI","SS_mARI","SS_MCARI1","SS_NBR2","SS_NDBI",
              "SS_WET","GS_NIR","GS_SWIR1","SS_blue","SS_NIR","SS_SWIR1",
              "elev","roughness","slope","twi","Geologia")

sel.pred.names <- c("LULC","GS_BSI","GS_CI","GS_GLI","GS_mARI","GS_NBR2","GS_NDBI","GS_NDWI",
                    "GS_SBI","SS_ARVI","SS_BSI","SS_CI","SS_CVI","SS_EVI","SS_GLI",
                    "SS_GVI","SS_GVMI","SS_mARI","SS_MCARI1","SS_NBR2","SS_NDBI",
                    "SS_WET","GS_NIR","GS_SWIR1","SS_blue","SS_NIR","SS_SWIR1",
                    "Elevation","Roughness","Slope","TWI","Geology")

#####################################################################
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

##################### RF tuning ####################################

df.train$L3 <- as.factor(df.train$L3)
df.train$Geologia <- as.factor(df.train$Geologia)


pngp.task = makeClassifTask(data = df.train, target = "L3")

resAcc = tuneRanger(
  pngp.task,
  measure = list(acc), #NULL,
  iters = 70,
  iters.warmup = 30,
  time.budget = NULL,
  num.threads = numCores,
  num.trees = 100,
  parameters = list(replace = FALSE, respect.unordered.factors = "order"),
  tune.parameters = c("mtry", "min.node.size", "sample.fraction"),
  save.file.path = paste0(rfFolder,"/RFtuning_Acc.RData"),
  build.final.model = TRUE,
  show.info = getOption("mlrMBO.show.info", TRUE))









###########################################################################
# Split 70% / 30%
hre_seed<- 50
set.seed(hre_seed)

intraining <- createDataPartition(df.train$L3, p = .70, list = FALSE)

training <- df.train[intraining,]
validation <- df.train[-intraining,]

dim(training)
dim(validation)
#

hyper_grid <- expand.grid(
  mtry       = 1:as.integer(sqrt((ncol(df.train)-1))),
  node_size  = 1:10,
  num.trees = seq(50,500,50),
  OOB_RMSE   = 0
)

system.time(
  for(i in 1:nrow(hyper_grid)) {
    # train model
    rf <- ranger(
      formula        = L3 ~ .,
      data           = df.train,
      num.trees      = hyper_grid$num.trees[i],
      mtry           = hyper_grid$mtry[i],
      min.node.size  = hyper_grid$node_size[i],
      importance = 'impurity')
    # add OOB error to grid
    hyper_grid$OOB_RMSE[i] <- sqrt(rf$prediction.error)
  })

nrow(hyper_grid) # number of models
position = which.min(hyper_grid$OOB_RMSE)
head(hyper_grid[order(hyper_grid$OOB_RMSE),],5)

# fit best model
rf.model <- ranger(L3 ~ .,data = df.train, num.trees = hyper_grid$num.trees[position], 
                   importance = 'impurity', probability = FALSE, 
                   min.node.size = hyper_grid$node_size[position], mtry = hyper_grid$mtry[position])




################ Cross-validation ####################################

trainx <- validation[,2:ncol(validation)]
trainy <- validation[,1]


n <- 20 # number of iterations, 60 to 100 is recommended.
measures <- NULL
for (i in 1:n) {
  rgcv1 <- rgcv(trainx,
                trainy, 
                mtry = 16, #if (!is.null(trainy) && !is.factor(trainy)) max(floor(ncol(trainx)/3), 1) else
                #floor(sqrt(ncol(trainx))), 
                num.trees = 100, 
                min.node.size = 2, 
                num.threads = numCores,
                verbose = TRUE,
                predacc = "ALL")
  measures <- rbind(measures, rgcv1$kappa) # for kappa, replace ccr with kappa
}
plot(measures ~ c(1:n), xlab = "Iteration for RF", ylab = "Correct classification rate (%)")
points(cumsum(measures) / c(1:n) ~ c(1:n), col = 2)
abline(h = mean(measures), col = 'blue', lwd = 2)
