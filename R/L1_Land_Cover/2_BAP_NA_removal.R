
# Load necessary libraries
library(terra)

# Folders ----
sito <- "PNGP"
root <- file.path(paste0("D:/ABRESO/",sito))
bap.folder <- file.path(root,"Landsat/level2/BAP")
l3.folder <- file.path(root,"Landsat/level3/BAP")
L2TC.folder <- file.path(root,"Landsat/level2/ATC/improved_cosine")

# Set working directory
setwd(bap.folder)

#### Retrieve available images and scores  ----
imgs.list <- list.files(path = L2TC.folder,pattern="T1")
f <- function(x){as.Date(substr(x,18,25),format="%Y%m%d")}
days.list <- lapply(imgs.list,FUN = f)
days.list <- sort(as.Date(unlist(days.list)))

first.date <- min(days.list)
last.date <- max(days.list)

yy <- (strftime(first.date, format = "%Y"):strftime(last.date, format = "%Y"))

bands <- c("blue","green","red","NIR","SWIR1","SWIR2")

# DEFINE FUNCIONS

# Function to fill NA values with the mean of the previous and following 2 layers
fill_na_with_mean_terra <- function(raster_stack) {
  # Number of layers in the stack
  nlayers_stack <- nlyr(raster_stack)
  
  # Create a copy of the stack to store the results
  filled_stack <- raster_stack
  
  # Iterate over each layer in the stack
  for (i in 1:nlayers_stack) {
    # Determine the range of layers to average (prev 2, current, next 2)
    start_layer <- max(1, i - 2)
    end_layer <- min(nlayers_stack, i + 2)
    
    # Extract the relevant subset of layers
    subset_stack <- raster_stack[[start_layer:end_layer]]
    
    # Calculate the mean across these layers, ignoring NA values
    mean_raster <- app(subset_stack, fun = function(x) mean(x, na.rm = TRUE))
    
    # Replace NA values in the current layer with the calculated mean
    current_layer <- raster_stack[[i]]
    lname <- names(current_layer)
    filled_layer <- ifel(is.na(current_layer), mean_raster, current_layer)
    names(filled_layer) <- lname
    
    # Store the filled layer back into the filled_stack
    filled_stack[[i]] <- filled_layer
  }
  
  return(filled_stack)
}

# Save each layer separately
save_filled_layers <- function(raster_stack, output_directory) {
  # Ensure the output directory exists, create it if not
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  
  # Number of layers in the stack
  nlayers_stack <- nlyr(raster_stack)
  
  # Iterate over each layer in the stack
  for (i in 1:nlayers_stack) {
    # Get the current layer
    current_layer <- raster_stack[[i]]
    namefile <- names(current_layer)
    # Define the filename for the output
    layer_name <- paste0(namefile, "_filled.tif")
    output_path <- file.path(output_directory, layer_name)
    
    # Save the current layer to a file
    writeRaster(current_layer, filename = output_path, filetype = "GTiff", overwrite = TRUE)
  }
}


# Growing season ----
for (b in bands) {
  
  ### Import BAP images
  setwd(bap.folder)
  landstack <- rast(list.files(pattern=glob2rx(paste0("*GS*",b,".tif$")))) 
  names(landstack) <- unlist(paste0(substr(list.files(pattern=glob2rx(paste0("*GS*",b,".tif$"))),1,28),b))
  # Fill NA and save
  filled_raster_stack <- fill_na_with_mean_terra(landstack)
  save_filled_layers(filled_raster_stack, "D:/ABRESO/PNGP/Landsat/level2/BAP/GF")
  
  print(paste0("NA filling and saving of band ",b," completed successfully."))
}

# Senescence season ----

for (b in bands) {
  ### Import BAP images
  setwd(bap.folder)
  landstack <- rast(list.files(pattern=glob2rx(paste0("*SS*",b,".tif$"))))
  names(landstack) <- unlist(paste0(substr(list.files(pattern=glob2rx(paste0("*SS*",b,".tif$"))),1,28),b))
  # Fill NA and save
  filled_raster_stack <- fill_na_with_mean_terra(landstack)
  save_filled_layers(filled_raster_stack, "D:/ABRESO/PNGP/Landsat/level2/BAP/GF")
  print(paste0("NA filling and saving of band ",b," completed successfully."))
}
      
      
# Inspect results
setwd("D:/ABRESO/PNGP/Landsat/level2/BAP/GF")
for (y in yy) {
  
  tci_gs <- rast(c(paste0(y,"_GS_195028_02_T1_SR_BAP_red_filled.tif"),
                   paste0(y,"_GS_195028_02_T1_SR_BAP_green_filled.tif"),
                   paste0(y,"_GS_195028_02_T1_SR_BAP_blue_filled.tif")))
  
  tci_ss <- rast(c(paste0(y,"_SS_195028_02_T1_SR_BAP_red_filled.tif"),
                   paste0(y,"_SS_195028_02_T1_SR_BAP_green_filled.tif"),
                   paste0(y,"_SS_195028_02_T1_SR_BAP_blue_filled.tif")))
  
  png(paste0(y,"_TCI.png"), width=10, height=5, units="in", res=300)
  par(mfrow=c(1,2),adj=0.3)
  plotRGB(tci_gs, r = 1, g = 2, b = 3, stretch = "lin", mar=1, main = "Growing season")
  plotRGB(tci_ss, r = 1, g = 2, b = 3, stretch = "lin", mar=1, main = "Senescence season")
  mtext(y, side = 3, line = - 2, outer = TRUE)
  dev.off() #only 129kb in size
  
  rm(tci_gs,tci_ss)
  print(paste("year",y,"done"))
  
}
