# --------------------------
# 1) Install / load packages
# --------------------------
required <- c("terra","sf","dplyr","caret","randomForest","e1071")
for (pkg in required) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)

library(terra)
library(sf)
library(dplyr)
library(caret)
library(randomForest)

# --------------------------
# 2) User inputs
# --------------------------
raster_path   <- "Chillaw_RGB.tif"
shapefile_path<- "SamplePoints.shp"
output_classified <- "rf_classified.tif"

# --------------------------
# 3) Read data
# --------------------------
rgb <- rast(raster_path)
print(rgb)

if (nlyr(rgb) < 3) stop("The input raster must have at least 3 bands (R,G,B).")

band_count <- nlyr(rgb)
if (is.null(names(rgb)) || any(grepl("^lyr|^band|^X", names(rgb)))) {
  names(rgb) <- paste0("B", 1:band_count)
}

samples <- vect(shapefile_path)
print(samples)

# --------------------------
# 4) Use 'Point_Type' column for classes (1=Vegetation, 2=Builtup, 3=Water)
# --------------------------
if (!"Point_Type" %in% names(samples)) {
  stop("The shapefile must contain a column named 'Point_Type'.")
}

samples$Point_Type <- factor(samples$Point_Type,
                             levels = c(1,2,3),
                             labels = c("Vegetation","Builtup","Water"))
names(samples)[names(samples) == "Point_Type"] <- "class"

# --------------------------
# 5) Reproject shapefile if needed
# --------------------------
if (crs(samples) != crs(rgb)) {
  samples <- project(samples, crs(rgb))
  message("Reprojected samples to raster CRS.")
}

# --------------------------
# 6) Extract raster values at sample locations
# --------------------------
vals <- extract(rgb, samples, ID=TRUE)
samples_df <- as.data.frame(samples)
training_df <- cbind(samples_df, vals[ , -1, drop = FALSE])
training_df <- na.omit(training_df)

training_df$class <- as.factor(training_df$class)
band_cols <- setdiff(names(vals), "ID")

# --------------------------
# 7) Train/test split
# --------------------------
set.seed(42)
trainIndex <- createDataPartition(training_df$class, p = 0.7, list = FALSE)
train_data <- training_df[trainIndex , c("class", band_cols)]
test_data  <- training_df[-trainIndex, c("class", band_cols)]

# --------------------------
# 8) Train Random Forest model
# --------------------------
rf_model <- randomForest(class ~ ., data = train_data, ntree = 500, importance = TRUE)
print(rf_model)

# --------------------------
# 9) Accuracy assessment
# --------------------------
pred_test <- predict(rf_model, newdata = test_data)
cm <- confusionMatrix(pred_test, test_data$class)
print(cm)

# --------------------------
# 10) Predict raster
# --------------------------
if (file.exists(output_classified)) file.remove(output_classified)

pred_raster <- terra::predict(
  rgb, rf_model,
  filename = output_classified,
  overwrite = TRUE,
  na.rm = TRUE
)

# --------------------------
# 11) Plot with custom colors
# --------------------------
class_labels <- c("Vegetation", "Builtup", "Water")
class_colors <- c("green", "pink", "blue")

plot(pred_raster, col = class_colors, legend = FALSE,
     main = "Random Forest Classified Map")
#legend("topright", legend = class_labels, fill = class_colors, border = "black", bg = "white")



# --------------------------
# Notes
# --------------------------
# - This script uses randomForest::randomForest for classification.
# - ntree=500 is a good default; increase for stability.
# - The output raster is written as 'rf_classified.tif'.
# - Colors: Vegetation=green, Builtup=pink, Water=blue.

# End of script
