install.packages("terra")
library(terra)

install.packages("sf")
library(sf)

#1)Download and Add bands 

#b2
b2_raster<-rast("T44NLP_20250305T045659_B02_20m.jp2")
print(b2_raster)
plot(b2_raster)

#b3
b3_raster<-rast("T44NLP_20250305T045659_B03_20m.jp2")
print(b3_raster)
plot(b3_raster)

#b4
b4_raster<-rast("T44NLP_20250305T045659_B04_20m.jp2")
print(b4_raster)
plot(b4_raster)

r<-c(b2_raster,b3_raster, b4_raster)
plot(r)

# 2)Composite bands

rgb_composite <- c(b4_raster, b3_raster, b2_raster)

plotRGB(rgb_composite, r = 1, g = 2, b = 3, 
        stretch = "lin", 
        main = "Sentinel-2 True-Color Composite (B4, B3, B2)")


writeRaster(rgb_composite, "sentinel_true_color_compressed.tif", 
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))


setwd("")
list.files()



#3) Add shapefile

Chillaw_shp<-vect("chillaw.shp")
print(Chillaw_shp)
plot(Chillaw_shp,col=NA, border="orange")

#check coordi. referenc. system
print(crs(rgb_composite))
print(crs(Chillaw_shp))

#reproject the shapefile
Chillaw_shp <- project(Chillaw_shp, rgb_composite)
print(crs(Chillaw_shp))

#check again the extent== maximum and minimum coordi.
print(ext(rgb_composite))
print(ext(Chillaw_shp))


#4) Clipped shapefile with composite

#crop & clip
composite_crop<- crop(rgb_composite,Chillaw_shp)
composite_clip<- mask(composite_crop,Chillaw_shp)

#plot raster
plotRGB(composite_clip, r=1, g=2, b=3, stretch="lin",
        main="Sentinel-2 Clipped RGB Composite")

#save
writeRaster(composite_clip, "Chillaw_RGB.tif",
            overwrite = TRUE,
            gdal = c("COMPRESS=LZW"))





#5)cloud mask

scl <- rast("T44NLP_20250305T045659_SCL_20m.jp2")
Chillaw_tiff<-rast("Chillaw_RGB.tif")

# resample scl to match composite resolution (10m)
scl_resampled <- resample(scl, Chillaw_tiff, method = "near")

# define cloud classes (7, 8, 9, 10, 11 = clouds)
cloud_mask <- scl_resampled %in% c(3, 7, 8, 9, 10, 11)

# apply mask
rgb_cloudfree <- mask(Chillaw_tiff, cloud_mask, maskvalues=1, updatevalue=NA)

# plot result
plotRGB(rgb_cloudfree, r=1, g=2, b=3, stretch="lin",
        main="Sentinel-2 RGB Cloud-Masked")

writeRaster(rgb_cloudfree, "Chilaw_RGBCloudfree.tif",
            overwrite = TRUE,
            gdal = c("COMPRESS=LZW"))



#6) Develop an Index

#NDVI
b4 <- rast("T44NLP_20250305T045659_B04_20m.jp2") 
b8 <- rast("T44NLP_20250305T045659_B8A_20m.jp2") 


# Reproject / resample to match your Chillaw_shp clipped raster
b4 <- crop(b4, Chillaw_shp) |> mask(Chillaw_shp)
b8 <- crop(b8, Chillaw_shp) |> mask(Chillaw_shp)


# NDVI formula: (NIR - Red) / (NIR + Red)
ndvi <- (b8 - b4) / (b8 + b4)

# Plot NDVI
plot(ndvi, col = rev(terrain.colors(100)), main = "NDVI - Chillaw")

# Save NDVI raster
writeRaster(ndvi, "Chillaw_NDVI.tif",
            overwrite = TRUE,
            gdal = c("COMPRESS=LZW"))

#MSI
b11<-rast("T44NLP_20250305T045659_B11_20m.jp2")
b8 <- rast("T44NLP_20250305T045659_B8A_20m.jp2")

b11 <- crop(b11, Chillaw_shp) |> mask(Chillaw_shp)
b8 <- crop(b8, Chillaw_shp) |> mask(Chillaw_shp)

msi<-b11/b8  # to check potential water stress

plot(msi, col = rev(heat.colors(100)), main = "MSI - Chillaw")

#VWSI-Vegetation Water Stress Index

msi[msi==0]<-NA
vwsi<- (1-ndvi)/msi
plot(vwsi, col = rev(terrain.colors(100)), 
     main = "Vegetation Water Stress Index(VWSI) in Chillaw")

rcl <- matrix(c(-Inf, 0.5, 1,   # Low stress
                0.5, 1.0, 2,    # Moderate stress
                1.0, Inf, 3),   # High stress
              ncol = 3, byrow = TRUE)


# Apply classification
vwsi_class <- classify(vwsi, rcl)

# Define labels and colors
class_labels <- c("Low Stress", "Moderate Stress", "High Stress")
class_colors <- c("green", "yellow", "red")

# Plot classified raster
plot(vwsi_class, col = class_colors, main = "Classified VWSI",
     axes = FALSE, legend = FALSE)

# Add legend
legend("topright", legend = class_labels,
       fill = class_colors, border = "black", title = "VWSI Class")
