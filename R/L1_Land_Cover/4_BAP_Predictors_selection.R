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

y <- 2009

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
aoi <- vect("D:/ABRESO/PNGP/Shapefile/PNGP.shp")
#reference<- shapefile("D:/ABRESO/PNGP/LULC/Training/tipologie_habitat_pngp_LULC_filtered.shp")
#reference$L3 <- as.factor(reference$L3)
lulc <- rast("D:/ABRESO/PNGP/LULC/Training/LULC_L3_5m.tif") %>% 
  mask(aoi) %>% resample(indices, method="near") 
lulc[lulc==0] <- NA

### Import static predictors ####
setwd("D:/ABRESO/PNGP/Ancillary/5m")
predlist <- list.files(getwd(), pattern=glob2rx(paste0('*.tif$')), all.files=TRUE, full.names=TRUE)
#predlist <- predlist[-6] # To remove Geology
pred.names <- lapply(predlist,function(x) {substr(x,29,nchar(x)-4)})
predictors <- rast(predlist) %>% resample(indices,method="near") 

all.names <- c("L3",ind.names,bands.names,pred.names)
all.pred <- c(lulc,indices,bands,predictors) %>% mask(aoi)
names(all.pred) <- unlist(all.names)

######################### Subset training #################################

# Compute z-statistics ----
# https://doi.org/10.1016/j.rse.2015.12.031

# To select only based on spectral bands

#all.names <- c("L3",bands.names)
#all.pred <- c(lulc,bands) %>% mask(aoi)
#names(all.pred) <- unlist(all.names)

df <- as.data.frame(all.pred)
df <- na.omit(df)

df$L3 <- as.factor(df$L3)
gc()

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
    writeRaster(z.stat.part,paste0(root,"/Landsat/z_statistics/",y,"_class_",f,"_variable_", names(l),"_zstat_30m_9cl.tif"),
                overwrite=TRUE,datatype="FLT4S",filetype="GTiff")
    n <- n+1
    r <- r+1
    rm(l,z.stat.part)
  }
  
  # compute final z-stat
  setwd(paste0(root,"/Landsat/z_statistics"))
  z.stats <- rast(list.files(pattern=glob2rx(paste0(y,"_class_",f,"_variable_*_zstat_30m_9cl.tif"))))
  f.z.stat <- function(x){return(sqrt(sum(x,na.rm=TRUE)))}
  z.stats <- app(z.stats, fun=f.z.stat, cores=8)
  z.stats <- mask(z.stats,lulc.f)
  z.stats[is.infinite(z.stats)] <- NA
  writeRaster(z.stats, filename=paste0(root,"/Landsat/z_statistics/",y,"_class_",f,"_TOT_zstat_30m_9cl.tif"), 
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
    writeRaster(training, paste0(root,"/Landsat/z_statistics/Training/",y,"_class_",f,"_training_30m_9cl.tif"), datatype='INT2S', filetype = "GTiff", overwrite=TRUE)
    
  } else {
    print("Sample size too small")
    training <- app(z.stats, fun = function(x) {return(!is.na(x))})
    training[training==1] <- as.numeric(f)
    training[training==0] <- NA
    writeRaster(training, paste0(root,"/Landsat/z_statistics/Training/",y,"_class_",f,"_training_30m_9cl.tif"), datatype='INT2S', filetype = "GTiff", overwrite=TRUE)
  }
  rm(z.stats,df.sub,data_long,values,lulc.f,my_plot,stat.sum,training,μ,σ,t1,t2)
  gc()
}

gc()


################################ Select predictors #################################

# Combine all the pixels selectd for each class to build the training dataset 
setwd("D:/ABRESO/PNGP/Landsat/z_statistics/Training")
trainings <- rast(list.files(pattern=glob2rx(paste0(y,"*training*30m*9cl*")))) %>% mask(aoi) 
all.trainings <- app(trainings, fun=sum, na.rm=TRUE) # Get training variables

rm(trainings)
gc()

# Import static predictors with geology
#setwd("D:/ABRESO/PNGP/Ancillary/5m")
#predlist <- list.files(getwd(), pattern=glob2rx(paste0('*.tif$')), all.files=TRUE, full.names=TRUE)
#pred.names <- lapply(predlist,function(x) {substr(x,29,nchar(x)-4)})
#predictors <- rast(predlist) %>% resample(indices) 

df.training <- c(all.trainings,indices,bands,predictors) %>% mask(aoi) 
all.names <- c("L3",ind.names,bands.names,pred.names)
names(df.training) <- unlist(all.names)

df.train <- as.data.frame(df.training) %>% 
  na.omit()

rm(df.training)
gc()

df.train$L3 <- as.factor(df.train$L3)
df.train$Geologia_5m <- as.factor(df.train$Geologia_5m)


# Naming variables
var.names <- c("LULC",
               unlist(ind.names),
               unlist(bands.names),
               "Aspect","CHM","Curvature","DAHI","Elevation","Geology","Roughness","Slope","TPI","TRI","TWI"   #unlist(pred.names)
               )
  
colnames(df.train) <- var.names

# Select a balanced sample size for each class
t <- df.train %>% group_by(LULC) %>%tally() %>% 
  as.data.frame()

l <- as.data.frame(lulc) %>% 
  na.omit()
tot <- nrow(l)

t$size <- round(as.numeric(t$n)*100/tot,1)

#sample.size <- min(t$n)
#reference <- df.train %>% group_by(L3) %>% slice_sample(n=sample.size) %>%  #balanced sampling
#as.data.frame()

reference <- df.train %>% group_by(LULC) %>% slice_sample(prop = 0.25) %>% 
  as.data.frame()

tt <- reference %>% group_by(LULC) %>%tally() %>% 
  as.data.frame()

tt$size.sample <- round(as.numeric(tt$n)*100/tot,1)

#fwrite(reference, file = paste0("D:/ABRESO/PNGP/LULC/predictors_selection_dataset_proportional_30m_9cl.csv"), sep = ",")

#### Divide ground references in training dataset (70%) and validation dataset (30%) ####
#reference <- as.data.frame(all.pred)
#reference$L3[reference$L3==0] <- NA
#reference <- na.omit(reference)
#reference$L3 <- as.factor(reference$L3)

hre_seed<- 50
set.seed(hre_seed)

intraining <- createDataPartition(reference$LULC, p = .70, list = FALSE)

training <- reference[intraining,]
validation <- reference[-intraining,]

dim(training)
dim(validation)

#### Mean and standard deviation 
mean_table <-aggregate(reference[, 2:ncol(reference)], list(reference$LULC), FUN=mean)
write.table(mean_table,paste("D:/ABRESO/PNGP/LULC/Mean_",y,"_30m_9cl.txt",sep=""),
            col.names = T,row.names = T,sep = "\t")

sd_table <-aggregate(reference[, 2:ncol(reference)], list(reference$LULC), sd)
write.table(sd_table,paste("D:/ABRESO/PNGP/LULC/Standard_deviation_",y,"_30m_9cl.txt",sep=""),
            col.names = T,row.names = T,sep = "\t")


################################ Predictor selection ################################
#### Area Under the receiver operator characteristic Curve (AUC) ####

d2<-reference[,-1]
d.topo <- d2[,67:76]
d.bands <- d2[,55:66]
d.ind <- d2[,1:55]
d<-reference[,1]
d3 <- reference[,-c(1,73)]

ROC<-filterVarImp(d3,d, nonpara = TRUE)
ROC

write.table(ROC,paste("D:/ABRESO/PNGP/LULC/AUC_",y,"_30m.txt",sep=""),col.names = T,row.names = T,sep = "\t")

roc.df <- as.data.frame(ROC) 
roc.df$name <- row.names(roc.df)

p <- ggplot(roc.df, aes(reorder(name, X1), y=X1)) + 
  geom_bar(stat = "identity", width=0.5, fill="#44546A") +
  coord_flip() + # Horizontal bar plot
  theme_minimal() +
  ylim(0, 1) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 6, hjust = 1, family = "TT Times New Roman"),
    axis.text.y = element_text(size = 6, hjust = 1, family = "TT Times New Roman")#,
    #plot.margin = margin(rep(15, 4))
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank()
  )  

p


#### Boruta ####
boruta_output <- Boruta(LULC ~., data=na.omit(reference), doTrace=0) 
boruta_signif <- names(boruta_output$finalDecision[boruta_output$finalDecision %in% c("Confirmed", "Tentative")])
print(boruta_signif)
plot <-plot(boruta_output, cex.axis=.6, las=2, xlab="", main="Variable Importance")  
plot

getSelectedAttributes(boruta_output, withTentative = F)
bot2<-attStats(boruta_output) 
print(bot2)

write.table(bot2,paste("D:/ABRESO/PNGP/LULC/Boruta_",y,"_30m.txt",sep=""),col.names = T,row.names = T,sep = "\t")

# Adjust the margins: (bottom, left, top, right)
par(mar = c(6, 5, 2, 0.5) + 0.1)

# Now, generate the plot again
plot(boruta_output, 
     cex.axis = .8, 
     las = 2, 
     xlab = "",  # Keep this empty
     ylab = "")#,  # Keep this empty
#main = "Boruta Feature Importance")

# Add x and y axis labels manually
#mtext("Features", side = 1, line = 6)  # x-axis label
mtext("Importance", side = 2, line = 4)  # y-axis label


#### Correlation matrix ####
correlation_matrix <- cor(d3)
print(correlation_matrix)

highly_correlated <- findCorrelation(correlation_matrix, cutoff=1)
print(highly_correlated)

corrplot(cor(correlation_matrix), method="color", type = "lower", tl.cex=0.5, tl.col="black")
corrplot(cor(correlation_matrix), method="number", type = "lower", tl.cex=0.5, tl.col="black", number.cex=0.5)

write.table(correlation_matrix,paste("D:/ABRESO/PNGP/LULC/correlation_",y,"_30m.txt",sep=""),col.names = T,row.names = T,sep = "\t")


#### Test the separability of the classes with Jeffries-Matusita
#install.packages("varSel")
library(varSel)

d2<-validation[,-1]
d.topo <- d2[,67:76]
d.bands <- d2[,55:66]
d.ind <- d2[,1:55]
d<-validation[,1]

#jm1 <- spectral.separability(d, d2, jeffries.matusita = TRUE)
jm2 <- JMdist(d, solve(d.bands))
