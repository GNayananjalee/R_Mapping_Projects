# SVM classification of an RGB GeoTIFF using point samples in R (terra + caret)
# Save this script and run in RStudio.


required <- c("terra","sf","dplyr","caret","kernlab","e1071")
for (pkg in required) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)

library(terra)
library(sf)
library(dplyr)
library(caret)
library(kernlab)
library(e1071)

# 2) User inputs
raster_path   <- "Chillaw_RGB.tif"       # RGB tif
shapefile_path<- "SamplePoints.shp"     # sample point shapefile



# 3) Read data
rgb <- rast(raster_path)
print(rgb)

if (nlyr(rgb) < 3) stop("The input raster must have at least 3 bands (R,G,B).")

band_count <- nlyr(rgb)
if (is.null(names(rgb)) || any(grepl("^lyr|^band|^X", names(rgb)))) {
  names(rgb) <- paste0("B", 1:band_count)
}

samples <- vect(shapefile_path)
print(samples)


# 4) Use the 'Point_Type' column for classes (1=Vegetation, 2=Builtup, 3=Water)
if (!"Point_Type" %in% names(samples)) {
  stop("The shapefile must contain a column named 'Point_Type'.")
}

# Convert numeric codes to factor with labels
samples$Point_Type <- factor(samples$Point_Type,
                             levels = c(1,2,3),
                             labels = c("Vegetation","Builtup","Water"))

# Rename to 'class' for convenience
names(samples)[names(samples) == "Point_Type"] <- "class"


# 5) Reproject raster and sample points to a common CRS
# Choose a standard CRS (e.g., UTM zone or WGS84). Here we use the raster CRS as the target.
if (is.na(crs(rgb))) {
  stop("The raster has no defined CRS. Please define it before continuing.")
}

if (crs(samples) != crs(rgb)) {
  samples <- project(samples, crs(rgb))
  message("Reprojected samples to raster CRS.")
}


# 6) Extract raster values at sample locations
vals <- extract(rgb, samples, ID=TRUE)
samples_df <- as.data.frame(samples)
training_df <- cbind(samples_df, vals[ , -1, drop = FALSE])
training_df <- na.omit(training_df)

# Ensure factor
training_df$class <- as.factor(training_df$class)

band_cols <- setdiff(names(vals), "ID")


# 7) Train/test split
set.seed(42)
trainIndex <- createDataPartition(training_df$class, p = 0.7, list = FALSE)
train_data <- training_df[trainIndex , c("class", band_cols)]
test_data  <- training_df[-trainIndex, c("class", band_cols)]


# 8) Train SVM model
svm_formula <- as.formula(paste("class ~ ", paste(band_cols, collapse = " + ")))

trctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 3)
svm_model <- train(svm_formula,
                   data = train_data,
                   method = "svmRadial",
                   trControl = trctrl,
                   preProcess = c("center", "scale"),
                   tuneLength = 6)

print(svm_model)


# 9) Accuracy assessment
pred_test <- predict(svm_model, newdata = test_data)
cm <- confusionMatrix(pred_test, test_data$class)
print(cm)




# 10) Predict raster

output_classified <- "svm_classified_v2.tif"
# Train the SVM using e1071
svm_model2 <- e1071::svm(class ~ ., data = train_data, kernel = "radial", scale = TRUE)

# If output file already exists, remove it
if (file.exists(output_classified)) file.remove(output_classified)

# Predict across the raster
pred_raster <- terra::predict(
  rgb, svm_model2,
  filename = output_classified,
  overwrite = TRUE,
  na.rm = TRUE
)

# Plot result
plot(pred_raster, main = "SVM classified: Vegetation, Builtup, Water")





# IF WANT TO DEFINE THE COLORS LIKE THIS, NEED TO USE BELOW CODE
# Define class labels and colors
class_labels <- c("Vegetation", "Builtup", "Water")
class_colors <- c("green", "pink", "blue")

# Plot with colors
plot(pred_raster,
     col = class_colors,
     legend = FALSE,
     main = "SVM Classified Map")

# Add legend
legend("topright",
       legend = class_labels,
       fill = class_colors,
       border = "black",
       bg = "white")




# --------------------------
# Notes
# --------------------------
# - The script now explicitly uses 'Point_Type' from SamplePoints.shp.
# - Values 1/2/3 are relabeled to Vegetation/Builtup/Water.
# - Shapefile and raster are reprojected to a common CRS.
# - For very large rasters, terra::predict writes to disk in chunks.
# - You can replace caret with e1071::svm if you want a simpler model.

# End of script
