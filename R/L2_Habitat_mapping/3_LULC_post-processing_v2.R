# Post-processing

# Rules:
# - If the pixel is assigned only once to a different class in the whole time 
# series, it is an error, replace it with the most frequent value;
# - if in 1 year the pixel belongs to one class, but in the 3 or more previous 
# and subsequent years it belongs to another class, reclassify it in the latter; 
# - if a pixel belongs more than half the time to class 3, reclassify it in all layers as 3;
# - if a pixel is NA, fill it with the most frequent value in the time series;
# - all other pixel values remain unchanged. 

library(terra)
library(dplyr)


# Import the raster stack
r <- rast(list.files("D:/ABRESO/PNGP/Landsat/RF/ensemble/habitat_mapping", 
                     pattern = "._PNGP_habitats.tif$", 
                     full.names = TRUE))

# Define a function that reclassifies a single pixel's time series
reclassify_pixel <- function(x) {
  # x: numeric vector of length 39 (one pixel across all years)
  
  ## 1) Fill NA with the most frequent non-NA value
  
  if (any(is.na(x))) {
    tbl_nonNA <- table(x[!is.na(x)])
    if (length(tbl_nonNA) > 0) {
      fill_class <- as.numeric(names(tbl_nonNA)[which.max(tbl_nonNA)])
      x[is.na(x)] <- fill_class
    }
  }
  
  ## 2) If a pixel belongs more than half the time to class 3, 
  ##    reclassify it in all layers as 3
  
  if (sum(x == 3, na.rm = TRUE) > (length(x) / 2)) {
    x[] <- 3
  }
  
  ## 3) If in 1 year the pixel is a different class than in the
  ##    3 previous & 3 subsequent years (which are the same), 
  ##    reclassify to that same neighbor class
  
  for (i in 4:(length(x) - 3)) {
    neighbors <- x[c(i-3, i-2, i-1, i+1, i+2, i+3)]
    # Skip if any neighbor is NA
    if (!any(is.na(neighbors))) {
      # Also skip if x[i] is NA (can't safely compare x[i] != neighbors[1])
      if (!is.na(x[i])) {
        # If all neighbors are the same
        if (length(unique(neighbors)) == 1) {
          # Reclassify if x[i] differs from the neighbor class
          if (x[i] != neighbors[1]) {
            x[i] <- neighbors[1]
          }
        }
      }
    }
  }
  
  ## 4) Replace single-occurrences with the most frequent value
  
  tbl <- table(x)
  if (length(tbl) > 0) {
    majority_class <- as.numeric(names(tbl)[which.max(tbl)])
    single_occurrences <- names(tbl)[tbl == 1]
    for (cls in single_occurrences) {
      cls_num <- as.numeric(cls)
      x[x == cls_num] <- majority_class
    }
  }
  
  
  return(x)
}


# Now apply the above function pixel-by-pixel over the entire 39-layer raster
r_reclassified <- app(r, reclassify_pixel)
names(r_reclassified) <- 1985:2023

# Finally, write out each of the 39 layers separately
# Modify the file names to suit your output preferences
out_dir <- "D:/ABRESO/PNGP/Landsat/RF/ensemble/habitat_mapping/postprocessing"
for (i in 1:nlyr(r_reclassified)) {
  yy <- names(r_reclassified[[i]])
  out_file <- file.path(out_dir, paste0(yy, "_habitat_map_reclassified.tif"))
  writeRaster(r_reclassified[[i]], filename = out_file, overwrite = TRUE)
}

