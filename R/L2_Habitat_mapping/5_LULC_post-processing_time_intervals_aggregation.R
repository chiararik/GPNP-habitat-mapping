# ------------------------------------------------------------------
# Load required library
library(terra)
library(dplyr)
library(ggplot2)
library(ggalluvial)

# ------------------------------------------------------------------
# 1. Load land cover time series stack
root <- file.path(paste0("D:/ABRESO/",sito))
rfFolder <- file.path(root,"Landsat/RF/ensemble/habitat_mapping/postprocessing")

# List all the land cover files in the directory
lc_stack <- rast(list.files(path=rfFolder, pattern=glob2rx("*.tif$"), full.names=TRUE))

# ------------------------------------------------------------------
# 2. Define the year vector and time intervals
years <- 1985:2023  # This should match the layers of your raster stack

# Find the indices of the layers belonging to each time interval
interval1_idx <- which(years >= 1985 & years <= 1995)  # 1985-1995
interval2_idx <- which(years >= 1995 & years <= 2005)  # 1995-2005
interval3_idx <- which(years >= 2005 & years <= 2015)  # 2005-2015
interval4_idx <- which(years >= 2015 & years <= 2023)  # 2015-2023

# ------------------------------------------------------------------
# 3. For each interval, calculate the mode (most frequent value) per pixel
#    'modal()' returns a SpatRaster where each pixel is the most common value
#    among the stacked layers supplied.

lc_1985_1995 <- modal(lc_stack[[interval1_idx]])
lc_1995_2005 <- modal(lc_stack[[interval2_idx]])
lc_2005_2015 <- modal(lc_stack[[interval3_idx]])
lc_2015_2023 <- modal(lc_stack[[interval4_idx]])

# ------------------------------------------------------------------
# 4. (Optional) Combine the results for the four intervals into a single stack
lc_intervals <- c(lc_1985_1995, lc_1995_2005, lc_2005_2015, lc_2015_2023)
names(lc_intervals) <- c("LC_1985_1995", "LC_1995_2005", 
                         "LC_2005_2015", "LC_2015_2023")

out_dir <- "D:/ABRESO/PNGP/Landsat/RF/ensemble/habitat_mapping/postprocessing/5Y"
for (i in 1:nlyr(lc_intervals)) {
  yy <- names(lc_intervals[[i]])
  out_file <- file.path(out_dir, paste0(yy, "_LC_reclassified.tif"))
  writeRaster(lc_intervals[[i]], filename = out_file, overwrite = TRUE)
}

# ------------------------------------------------------------------
# 5. (Optional) Write out to file
writeRaster(lc_intervals, "landcover_intervals_mode.tif", overwrite=TRUE)

classifyPixel <- function(x) {
  # x is a numeric vector of length 39 for a single pixel (time series)
  
  # Remove NAs if needed (uncomment if you want to handle NA specifically)
  # x <- x[!is.na(x)]
  # if (length(x) == 0) return(NA)  # Or a specific code for all-NA
  
  # 1) Check how many unique values
  uv <- unique(x)
  n_unique <- length(uv)
  
  # --- CASE 1: Stable ---
  if (n_unique == 1) {
    # All values in time series are the same
    return(1)
  }
  
  # --- CASE 2 or 4: Exactly 2 unique values ---
  if (n_unique == 2) {
    # Count how many times each value occurs
    freq <- table(x)
    
    # **Error (2)**: 38 of one value, 1 of the other
    if (any(freq == 38) && any(freq == 1)) {
      return(2)
    } else {
      # Count how many times the value changes between consecutive years
      transitions <- sum(x[-1] != x[-length(x)])
      # **Transition (4)**: Only one change event
      if (transitions == 1) {
        return(4)
      } else {
        # **Mixed (3)**: More than one change
        return(3)
      }
    }
  }
  
  # --- CASE 3: > 2 unique values => Mixed ---
  return(3)
}

# Apply classifyPixel to each pixel:
classified_rast <- app(lc_intervals, classifyPixel) %>% 
  mask(lc_intervals[[1]])

setwd(file.path(root,"Landsat/RF/ensemble/habitat_mapping/postprocessing/5Y"))
writeRaster(classified_rast, paste0("Pixel_types.tif"), datatype='INT2S', filetype = "GTiff", overwrite=TRUE)

# ------------------------------------------------------------------
# Result:
# 'lc_intervals' is a SpatRaster with four layers, each representing
# the most frequent (mode) land cover class per pixel for the respective
# time interval.

# -------------------------------------------------------------------
# -----------------------------------------------------------
# -----------------------------------------------------------
# 2. Suppose 'lc_intervals' is a SpatRaster with four layers:
#     1) LC_1985_1995
#     2) LC_1995_2005
#     3) LC_2005_2015
#     4) LC_2015_2023
#
#    Each layer has integer codes 1..8 (land-cover classes).
#    We convert these to a data frame of pixel values.

classified_rast[classified_rast!=4]<- NA
lc_intervals <- mask(lc_intervals, classified_rast)

df_raw <- as.data.frame(lc_intervals, xy=FALSE)
names(df_raw) <- c("t1", "t2", "t3", "t4")

# Remove rows with NA (optional) if you don't want incomplete pixels
df_raw <- df_raw[complete.cases(df_raw), ]

# -----------------------------------------------------------
# 3. Count how many pixels follow each combination of classes
transition_counts <- df_raw %>%
  group_by(t1, t2, t3, t4) %>%
  summarize(count = n(), .groups = "drop")

# -----------------------------------------------------------
# 4. Convert pixel counts to area (hectares) 
#    (assuming your raster is in meters, e.g. a projected CRS)
px_res    <- res(lc_intervals)             # pixel size in X and Y (m)
px_area_m <- abs(px_res[1] * px_res[2])    # pixel area (m^2)
px_area_ha <- px_area_m / 10000            # convert m^2 to hectares

transition_counts <- transition_counts %>%
  mutate(area_ha = count * px_area_ha)

# -----------------------------------------------------------
# 5. Prepare for ggalluvial:
#    - Rename the time-slice columns to descriptive strings
#    - Convert the numeric classes (1..8) into factor with c.labels

df_alluvial <- transition_counts %>%
  rename(
    "1985-1995" = t1,
    "1995-2005" = t2,
    "2005-2015" = t3,
    "2015-2023" = t4
  ) %>%
  filter(area_ha > 0)

# Convert each time-slice to a factor with the same levels & labels
df_alluvial$`1985-1995` <- factor(df_alluvial$`1985-1995`, 
                                  #levels = 1:8, 
                                  #labels = c.labels
                                  )
df_alluvial$`1995-2005` <- factor(df_alluvial$`1995-2005`, 
                                  #levels = 1:8, 
                                  #labels = c.labels
                                  )
df_alluvial$`2005-2015` <- factor(df_alluvial$`2005-2015`, 
                                  #levels = 1:8, 
                                  #labels = c.labels
                                  )
df_alluvial$`2015-2023` <- factor(df_alluvial$`2015-2023`, 
                                  #levels = 1:8, 
                                  #labels = c.labels
                                  )

# -----------------------------------------------------------
# 6. Plot the alluvial diagram
ggplot(
  df_alluvial,
  aes(
    axis1 = `1985-1995`,
    axis2 = `1995-2005`,
    axis3 = `2005-2015`,
    axis4 = `2015-2023`,
    y = area_ha
  )
) +
  # The flows (alluvia) are colored by the class in the first axis 
  geom_alluvium(aes(fill = `1985-1995`), width = 1/12, alpha = 0.8) +
  
  # The 'strata' are the blocks at each axis; color them by class as well
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/12, color = "grey") +
  
  # Label each stratum block with its factor label
  #geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            #size = 3, color = "black") +
  
  # Control the order of discrete x-axis categories
  scale_x_discrete(limits = c("1985-1995", "1995-2005", "2005-2015", "2015-2023")) +
  
  # Map our factor levels to the provided color palette
  #scale_fill_manual(values = c.values, name = "Land cover") +
  
  # Axis labels, theme, etc.
  ylab("Area (ha)") +
  xlab("Time intervals") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title.x    = element_text(size=12),
    axis.title.y    = element_text(size=12)
  )



