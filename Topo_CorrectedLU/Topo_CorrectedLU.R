
# Install / load packages
req <- c("terra","sf","randomForest","caret")
inst <- req %in% installed.packages()[, "Package"]
if (any(!inst)) install.packages(req[!inst], repos="https://cloud.r-project.org")
library(terra); library(sf); library(randomForest); library(caret)




# Filenames (change to your actual filenames)
srtm_file  <- "output_SRTMGL1.tif"    # SRTM DEM
study_shp  <- "CH1.shp"               # study area polygon (CH1)
train_shp  <- "landuse_training.shp"
print (train_shp)
plot(train_shp)# training polygons/points (field 'class')
# Landsat band files (adjust to your set and correct order)
band_files <- c("LC08_L2SP_141055_20250204_20250208_02_T1_SR_B2.TIF","LC08_L2SP_141055_20250204_20250208_02_T1_SR_B3.TIF","LC08_L2SP_141055_20250204_20250208_02_T1_SR_B4.TIF","LC08_L2SP_141055_20250204_20250208_02_T1_SR_B5.TIF","LC08_L2SP_141055_20250204_20250208_02_T1_SR_B6.TIF","LC08_L2SP_141055_20250204_20250208_02_T1_SR_B7.TIF")

# Output filenames
out_corrected_stack <- "landsat_corrected.tif"
out_classified       <- "landsat_classified.tif"

# Solar geometry (fallback). If you have MTL, code below will try to read it.
solar_zenith_deg  <- 30   # fallback (degrees)
solar_azimuth_deg  <- 140 # fallback (degrees)

# training class field name in shapefile
class_field <- "Class"    # change if different
# -----------------------------------------------------------------------------

# 1) Read data
cat("Reading data...\n")
srtm  <- rast(srtm_file)
bands <- rast(band_files)           # multi-layer SpatRaster
study <- vect(study_shp)            # terra vector for cropping/masking
train_sf <- st_read(train_shp, quiet=TRUE)  # sf to access attributes

# 2) Ensure all have same CRS: reproject SRTM and study/train to bands' CRS
bands_crs <- crs(bands)
if (is.na(bands_crs)) stop("Bands have no CRS — set CRS and retry.")

if (!compareGeom(srtm, bands, stopOnError = FALSE)) {
  cat("Reprojecting SRTM to match Landsat CRS...\n")
  srtm <- project(srtm, crs(bands))
}

if (crs(study) != crs(bands)) {
  cat("Reprojecting study polygon to match Landsat CRS...\n")
  study <- project(study, crs(bands))
}
# Reproject training shapefile (sf -> terra vect)
train_sf <- st_transform(train_sf, crs(bands))
train_vect <- vect(train_sf)

# 3) Crop & mask rasters to study area (safe checks)
# Make sure they overlap
if (!relate(bands, study, relation="intersects")) {
  stop("ERROR: bands and study polygon do not overlap after reprojection. Check extents/CRS.")
}

bands <- crop(bands, study); bands <- mask(bands, study)
srtm  <- crop(srtm, study);  srtm  <- mask(srtm, study)

# 4) Compute slope & aspect in radians from SRTM
cat("Computing slope & aspect from SRTM...\n")
slope_deg  <- terrain(srtm, v="slope", unit="degrees")
aspect_deg <- terrain(srtm, v="aspect", unit="degrees")
slope_rad  <- slope_deg * pi / 180
aspect_rad <- aspect_deg * pi / 180

# 5) Attempt to read solar angles from MTL if present (optional)
# If there is an MTL file in the workspace, try to parse SUN_AZIMUTH / SUN_ELEVATION
mtl_files <- list.files(workspace, pattern = "_MTL.txt$", full.names = TRUE)
if (length(mtl_files) >= 1) {
  cat("MTL file detected. Attempting to read solar angles from MTL...\n")
  mtl_text <- readLines(mtl_files[1])
  get_mtl_val <- function(tag) {
    line <- mtl_text[grepl(tag, mtl_text)]
    if (length(line)==0) return(NA)
    # extract value after =, remove quotes/spaces
    v <- sub(".*=\\s*", "", line[1])
    v <- gsub('"', '', v); v <- trimws(v)
    as.numeric(v)
  }
  sun_el <- get_mtl_val("SUN_ELEVATION")
  sun_az <- get_mtl_val("SUN_AZIMUTH")
  if (!is.na(sun_el)) {
    solar_zenith_deg <- 90 - sun_el
    cat(sprintf("Using solar zenith from MTL: %.3f deg (SUN_ELEVATION %.3f)\n", solar_zenith_deg, sun_el))
  }
  if (!is.na(sun_az)) {
    solar_azimuth_deg <- sun_az
    cat(sprintf("Using solar azimuth from MTL: %.3f deg\n", solar_azimuth_deg))
  }
} else {
  cat("No MTL file found — using fallback solar angles.\n")
}

solar_zenith_rad <- solar_zenith_deg * pi / 180
solar_azimuth_rad <- solar_azimuth_deg * pi / 180

# 6) Compute cos(i) using lapp() — then align cos_i to bands grid
cat("Computing illumination cos(i) and aligning to bands...\n")
cos_i_raw <- lapp(c(slope_rad, aspect_rad), fun = function(s, a) {
  cos(s) * cos(solar_zenith_rad) + sin(s) * sin(solar_zenith_rad) * cos(solar_azimuth_rad - a)
})

# Align cos_i to band grid exactly (resolution, extent & origin)
cos_i <- resample(cos_i_raw, bands[[1]], method = "bilinear")
cos_i <- crop(cos_i, bands[[1]])
cos_i <- mask(cos_i, bands[[1]])

# Verify alignment
if (!compareGeom(bands[[1]], cos_i, stopOnError = FALSE)) {
  stop("ERROR: cos_i and bands are not aligned after resampling. Inspect inputs.")
}

# 7) C-correction per band (safe handling of empty samples)
cat("Running per-band C-correction (C-correction)...\n")
corrected_layers <- list()
band_names <- names(bands)

for (i in seq_along(band_names)) {
  b <- bands[[i]]
  cat(sprintf("Processing band %d: %s\n", i, band_names[i]))
  
  # Combine band and cos_i into two-layer raster, convert to df removing NA rows
  combo <- c(b, cos_i)
  names(combo) <- c("L", "cosi")
  df <- as.data.frame(combo, na.rm=TRUE)
  
  if (nrow(df) == 0) {
    warning(sprintf("Band %s: no valid (non-NA) overlap with cos_i — skipping this band.", band_names[i]))
    next
  }
  
  # Subsample to a reasonable size for regression if huge
  if (nrow(df) > 200000) df <- df[sample(nrow(df), 200000), ]
  
  # Fit linear model L = a + b*cosi (C-correction uses C = a/b)
  lmfit <- tryCatch(lm(L ~ cosi, data=df), error=function(e) NULL)
  if (is.null(lmfit)) {
    warning(sprintf("Band %s: lm failed — skipping.", band_names[i])); next
  }
  a_coef <- coef(lmfit)[1]
  b_coef <- coef(lmfit)[2]
  if (is.na(a_coef) || is.na(b_coef) || abs(b_coef) < 1e-12) {
    warning(sprintf("Band %s: invalid regression coefficients — skipping correction (a=%.3e b=%.3e).", band_names[i], a_coef, b_coef))
    next
  }
  Cval <- a_coef / b_coef
  message(sprintf("Band %s: intercept=%.4e slope=%.4e C=%.4f", band_names[i], a_coef, b_coef, Cval))
  
  # Apply C-correction using lapp; protect divisions
  Lcorr <- lapp(c(b, cos_i), fun = function(Li, ci) {
    num <- cos(solar_zenith_rad) + Cval
    den <- ci + Cval
    out <- Li * (num / den)
    out[is.infinite(out) | is.nan(out)] <- NA
    out
  })
  
  names(Lcorr) <- paste0(band_names[i], "_c")
  corrected_layers[[length(corrected_layers) + 1]] <- Lcorr
}

# Combine corrected layers (ensure at least one)
if (length(corrected_layers) == 0) stop("No corrected bands produced — cannot continue.")
corrected_stack <- rast(corrected_layers)
cat("Writing corrected band stack to:", out_corrected_stack, "\n")
writeRaster(corrected_stack, out_corrected_stack, overwrite=TRUE)

names(train_sf)

# 8) Prepare training samples: pixel-level samples from training polygons
cat("Extracting training samples from corrected stack...\n")
# Ensure class field exists
if (!(class_field %in% names(train_sf))) stop("Training shapefile missing class field: ", class_field)
# Use terra::extract to get pixel-level samples; df=TRUE returns rows per pixel with ID linking to polygon
pix_samples <- terra::extract(corrected_stack, train_vect, df=TRUE)

# pix_samples has column ID (polygon index), and band columns; attach class label
poly_df <- data.frame(ID = seq_len(nrow(train_sf)), class = train_sf[[class_field]])
pix_samples2 <- merge(pix_samples, poly_df, by.x="ID", by.y="ID", all.x=TRUE)

# Identify corrected band columns
bandcols <- grep("_c$", names(pix_samples2), value=TRUE)
# remove rows with NA spectral data or NA class
pix_samples2 <- pix_samples2[complete.cases(pix_samples2[, bandcols]) & !is.na(pix_samples2$class), ]

if (nrow(pix_samples2) == 0) stop("No valid training pixels found. Check training shapefile and corrected stack overlap.")

# Convert class to factor
pix_samples2$class <- as.factor(pix_samples2$class)

# Balance dataset by downsampling majority classes (optional but recommended)
class_counts <- table(pix_samples2$class)
min_n <- min(class_counts)
set.seed(42)
balanced_samples <- do.call(rbind, lapply(levels(pix_samples2$class), function(cl) {
  rows <- which(pix_samples2$class == cl)
  if (length(rows) > min_n) rows <- sample(rows, min_n)
  pix_samples2[rows, ]
}))

train_df <- balanced_samples[, c(bandcols, "class")]
train_df <- train_df[complete.cases(train_df), ]



# Recode your training data classes explicitly
train_sf[[class_field]] <- factor(train_sf[[class_field]],
                                  levels = c(1, 2, 3, 4),
                                  labels = c("Vegetation", "WaterBodies", "BuiltUp", "PlantationArea"))



# FIX: Recode class levels BEFORE training and prediction


# Recode class labels to ensure correct mapping BEFORE sampling
train_sf[[class_field]] <- factor(train_sf[[class_field]],
                                  levels = c(1, 2, 3, 4),
                                  labels = c("Vegetation", "WaterBodies", "BuiltUp", "PlantationArea"))

# Update the merged pixel-level samples (so that factor order is consistent)
pix_samples2 <- merge(pix_samples, 
                      data.frame(ID = seq_len(nrow(train_sf)), class = train_sf[[class_field]]),
                      by.x="ID", by.y="ID", all.x=TRUE)

bandcols <- grep("_c$", names(pix_samples2), value=TRUE)
pix_samples2 <- pix_samples2[complete.cases(pix_samples2[, bandcols]) & !is.na(pix_samples2$class), ]

# Balance dataset
class_counts <- table(pix_samples2$class)
min_n <- min(class_counts)
set.seed(42)
balanced_samples <- do.call(rbind, lapply(levels(pix_samples2$class), function(cl) {
  rows <- which(pix_samples2$class == cl)
  if (length(rows) > min_n) rows <- sample(rows, min_n)
  pix_samples2[rows, ]
}))
train_df <- balanced_samples[, c(bandcols, "class")]
train_df <- train_df[complete.cases(train_df), ]

# rain Random Forest

cat("Training Random Forest classifier...\n")
set.seed(100)
rf_model <- randomForest(x = train_df[, bandcols], y = train_df$class, ntree = 500, importance = TRUE)
print(rf_model)
cat("OOB confusion:\n"); print(rf_model$confusion)


# Predict & generate final classified raster

cat("Predicting land use classes...\n")
pred <- terra::predict(corrected_stack, rf_model, type="response", na.rm=TRUE,
                       filename=out_classified, overwrite=TRUE)

classified <- rast(out_classified)

# Apply class names and colors correctly
class_names <- c("Vegetation", "WaterBodies", "BuiltUp", "PlantationArea")
cols <- c("forestgreen", "steelblue", "grey40", "darkolivegreen3")

# Apply proper legend info
levels(classified) <- data.frame(value = 1:4, class = class_names)

# Plot with correct colors and labeled legend
plot(classified, col = cols, legend = FALSE,
     main = "Topographically Corrected Land Use Classification in Central Highland")
legend("topright",
       legend = class_names,
       fill = cols,
       border = "black",
       cex = 0.7,       # make text and box 30% smaller
       box.lwd = 0.8,   # thinner box border
       bg = "white")  
legend("topright",
       legend = class_names,
       fill = cols,
       border = "black",
       cex = 1.2,        # Adjusts text and box size
       text.col = "black",
       bg = "white",     # Legend background
       box.lwd = 0.5,    # Thin border
       x.intersp = 0.7,  # Adjust space between box and text
       y.intersp = 0.7,  # Adjust vertical spacing
       bty = "o")  


