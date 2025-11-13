# --------------------------
# 1) Install / load packages
# --------------------------
required <- c("terra","sf","dplyr","caret","nnet")
for (pkg in required) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)

library(terra)
library(sf)
library(dplyr)
library(caret)
library(nnet)

# --------------------------
# 2) User inputs
# --------------------------
raster_path   <- "Chilaw_RGBCloudfree.tif"
shapefile_path<- "SamplePoints.shp"
output_classified <- "nn_classified.tif"

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
# 4) Use 'Point_Type' column for classes
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
# 8) Train Neural Network model
# --------------------------
ctrl <- trainControl(method = "cv", number = 5)  # 5-fold CV
nn_model <- train(
  class ~ .,
  data = train_data,
  method = "nnet",
  trControl = ctrl,
  tuneGrid = expand.grid(size = c(3,5,7), decay = c(0.01, 0.1)),
  trace = FALSE,
  maxit = 500
)

print(nn_model)

# --------------------------
# 9) Accuracy assessment
# --------------------------
pred_test <- predict(nn_model, newdata = test_data)
cm <- confusionMatrix(pred_test, test_data$class)
print(cm)

# --------------------------
# 10) Predict raster
# --------------------------
if (file.exists(output_classified)) file.remove(output_classified)

pred_raster <- terra::predict(
  rgb, nn_model,
  filename = output_classified,
  overwrite = TRUE,
  na.rm = TRUE,
  type = "raw"   # caret model requires raw prediction
)

# --------------------------
# 11) Plot with custom colors
# --------------------------
class_labels <- c("Vegetation", "Builtup", "Water")
class_colors <- c("green", "pink", "blue")

plot(pred_raster, col = class_colors, legend = TRUE,
     main = "Neural Network Classified Map")

# default placement (topright)
terra::north()  

terra::north(xy = "topright")

legend("topright", legend = class_labels, fill = class_colors, cex = 0.8)
