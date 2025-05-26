library(terra)
library(geometry)

# Set directory structure
spPath <- "E:/Research/DroughtForecasts/Data/Occurrences/BySpecies/Filtered/Over10/"
spList <- "E:/Research/PhylogeneticSignal/Data/RangeSizes/rangeSizes.csv"
predPath <- "E:/Research/DroughtForecasts/Data/Predictors/"
outPath <- "E:/Research/PhylogeneticSignal/Data/NicheSizes/"
PCAPAth <- "E:/Research/PhylogeneticSignal/Data/PCA/"
densityPath <- "E:/Research/PhylogeneticSignal/Data/DensityRasters/"

# List species
spList <- read.csv(spList, header = T)
spList <- unique(spList$species)

# List predictor files
predFiles <- list.files(predPath, "current")

# Run spatial PCA
set.seed(2025)
spatSample1 <- spatSample(rast(paste0(predPath, predFiles[1])), 100000, "random", na.rm = T)
spatSample2 <- spatSample(rast(paste0(predPath, predFiles[2])), 100000, "random", na.rm = T)
spatSample3 <- spatSample(rast(paste0(predPath, predFiles[3])), 100000, "random", na.rm = T)
spatSample1 <- spatSample1[complete.cases(spatSample1),]
spatSample2 <- spatSample2[complete.cases(spatSample2),]
spatSample3 <- spatSample3[complete.cases(spatSample3),]
spatPCA1 <- prcomp(spatSample1, center = T, scale. = T)
spatPCA2 <- prcomp(spatSample2, center = T, scale. = T)
spatPCA3 <- prcomp(spatSample3, center = T, scale. = T)

# Export PCA
save(spatSample1, spatSample2, spatSample3, spatPCA1, spatPCA2, spatPCA3, file = paste0(PCAPAth, "PCA.rda"))
# load(paste0(PCAPAth, "PCA.rda"))

# Read in predictors
preds1 <- rast(paste0(predPath, predFiles[1]))
preds2 <- rast(paste0(predPath, predFiles[2]))
preds3 <- rast(paste0(predPath, predFiles[3]))

# Transform predictors
predsAggr1 <- aggregate(preds1, fact = 2, cores = 4, na.rm = T)
predsTrans1 <- predict(spatPCA1, newdata = values(predsAggr1, na.rm = T))
predsTrans1 <- predsTrans1[, 1:2]
predsTrans1 <- as.data.frame(predsTrans1)
predsAggr2 <- aggregate(preds2, fact = 2, cores = 4, na.rm = T)
predsTrans2 <- predict(spatPCA2, newdata = values(predsAggr2, na.rm = T))
predsTrans2 <- predsTrans2[, 1:2]
predsTrans2 <- as.data.frame(predsTrans2)
predsAggr3 <- aggregate(preds3, fact = 2, cores = 4, na.rm = T)
predsTrans3 <- predict(spatPCA3, newdata = values(predsAggr3, na.rm = T))
predsTrans3 <- predsTrans3[, 1:2]
predsTrans3 <- as.data.frame(predsTrans3)

# Create density rasters
rast1 <- rast(nrows = 100, ncols = 100, xmin = min(predsTrans1[,1]), xmax = max(predsTrans1[,1]), ymin = min(predsTrans1[,2]), ymax = max(predsTrans1[,2]))
crs(rast1) <- ""
resx <- res(rast1)[1]
resy <- res(rast1)[2]
for (i in 1:ncell(rast1)) {
  PC1min <- xFromCell(rast1, i) - (resx / 2)
  PC1max <- xFromCell(rast1, i) + (resx / 2)
  PC2min <- yFromCell(rast1, i) - (resy / 2)
  PC2max <- yFromCell(rast1, i) + (resy / 2)
  count <- subset(predsTrans1, PC1 > PC1min & PC1 < PC1max & PC2 > PC2min & PC2 < PC2max)
  values(rast1)[i] <- nrow(count)
  print(nrow(count))
  plot(rast1)
}

rast2 <- rast(nrows = 100, ncols = 100, xmin = min(predsTrans2[,1]), xmax = max(predsTrans2[,1]), ymin = min(predsTrans2[,2]), ymax = max(predsTrans2[,2]))
crs(rast2) <- ""
resx <- res(rast2)[1]
resy <- res(rast2)[2]
for (i in 1:ncell(rast2)) {
  PC1min <- xFromCell(rast2, i) - (resx / 2)
  PC1max <- xFromCell(rast2, i) + (resx / 2)
  PC2min <- yFromCell(rast2, i) - (resy / 2)
  PC2max <- yFromCell(rast2, i) + (resy / 2)
  count <- subset(predsTrans2, PC1 > PC1min & PC1 < PC1max & PC2 > PC2min & PC2 < PC2max)
  values(rast2)[i] <- nrow(count)
  print(nrow(count))
  plot(rast2)
}

rast3 <- rast(nrows = 100, ncols = 100, xmin = min(predsTrans3[,1]), xmax = max(predsTrans3[,1]), ymin = min(predsTrans3[,2]), ymax = max(predsTrans3[,2]))
crs(rast3) <- ""
resx <- res(rast3)[1]
resy <- res(rast3)[2]
for (i in 1:ncell(rast3)) {
  PC1min <- xFromCell(rast3, i) - (resx / 2)
  PC1max <- xFromCell(rast3, i) + (resx / 2)
  PC2min <- yFromCell(rast3, i) - (resy / 2)
  PC2max <- yFromCell(rast3, i) + (resy / 2)
  count <- subset(predsTrans3, PC1 > PC1min & PC1 < PC1max & PC2 > PC2min & PC2 < PC2max)
  values(rast3)[i] <- nrow(count)
  plot(rast3)
}

# Export density rasters
writeRaster(rast1, paste0(densityPath, "rast1.tif"))
writeRaster(rast2, paste0(densityPath, "rast2.tif"))
writeRaster(rast3, paste0(densityPath, "rast3.tif"))

nicheSizes <- data.frame(species = spList,
                         nicheSize_GFDL = numeric(length(spList)),
                         nicheSize_MPI = numeric(length(spList)),
                         nicheSize_UKESM = numeric(length(spList)),
                         nicheDensity_GFDL = numeric(length(spList)),
                         nicheDensity_MPI = numeric(length(spList)),
                         nicheDensity_UKESM = numeric(length(spList)))
# Perform PCA and calculate niche size/density
for (i in 1:length(spList)) {
  print(spList[i])
  # Read occurrences
  spOcc <- vect(paste0(spPath, spList[i], ".shp"))
  
  # Extract environmental variables
  extr1 <- terra::extract(preds1, spOcc, cells = F, ID = F)
  extr1 <- extr1[complete.cases(extr1),]
  # Transform predictors
  trans1 <- predict(spatPCA1, newdata = extr1)
  # Calculate niche size
  trans1 <- trans1[, 1:2]
  convex1 <- convhulln(trans1, options = "FA")
  size1 <- convex1$vol
  # Calculate niche density
  hull_indices <- unique(as.vector(convex1$hull))
  hull_coords <- trans1[hull_indices, ]
  center <- colMeans(hull_coords)
  angles <- atan2(hull_coords[, 2] - center[2], hull_coords[, 1] - center[1])
  ordered_coords <- hull_coords[order(angles), ]
  if (!all(ordered_coords[1, ] == ordered_coords[nrow(ordered_coords), ])) {
    ordered_coords <- rbind(ordered_coords, ordered_coords[1, ])
  }
  polygon1 <- vect(ordered_coords, type = "polygons")
  density1 <- extract(rast1, polygon1, fun = sum, na.rm = TRUE)

  # Extract environmental variables
  extr2 <- terra::extract(preds1, spOcc, cells = F, ID = F)
  extr2 <- extr2[complete.cases(extr2),]
  # Transform predictors
  trans2 <- predict(spatPCA2, newdata = extr2)
  # Calculate niche size
  trans2 <- trans2[, 1:2]
  convex2 <- convhulln(trans2, options = "FA")
  size2 <- convex2$vol
  # Calculate niche density
  hull_indices <- unique(as.vector(convex2$hull))
  hull_coords <- trans1[hull_indices, ]
  center <- colMeans(hull_coords)
  angles <- atan2(hull_coords[, 2] - center[2], hull_coords[, 1] - center[1])
  ordered_coords <- hull_coords[order(angles), ]
  if (!all(ordered_coords[1, ] == ordered_coords[nrow(ordered_coords), ])) {
    ordered_coords <- rbind(ordered_coords, ordered_coords[1, ])
  }
  polygon2 <- vect(ordered_coords, type = "polygons")
  density2 <- extract(rast2, polygon2, fun = sum, na.rm = TRUE)
  
  # Extract environmental variables
  extr3 <- terra::extract(preds3, spOcc, cells = F, ID = F)  
  extr3 <- extr3[complete.cases(extr3),]
  # Transform predictors
  trans3 <- predict(spatPCA3, newdata = extr3)
  # Calculate niche size
  trans3 <- trans3[, 1:2]
  convex3 <- convhulln(trans3, options = "FA")
  size3 <- convex3$vol
  # Calculate niche density
  hull_indices <- unique(as.vector(convex3$hull))
  hull_coords <- trans3[hull_indices, ]
  center <- colMeans(hull_coords)
  angles <- atan2(hull_coords[, 2] - center[2], hull_coords[, 1] - center[1])
  ordered_coords <- hull_coords[order(angles), ]
  if (!all(ordered_coords[1, ] == ordered_coords[nrow(ordered_coords), ])) {
    ordered_coords <- rbind(ordered_coords, ordered_coords[1, ])
  }
  polygon3 <- vect(ordered_coords, type = "polygons")
  density3 <- extract(rast3, polygon3, fun = sum, na.rm = TRUE)
  
  # Export niche sizes
  nicheSizes[i, "species"] <- spList[i]
  nicheSizes[i, "nicheSize_GFDL"] <- size1
  nicheSizes[i, "nicheSize_MPI"] <- size2
  nicheSizes[i, "nicheSize_UKESM"] <- size3
  nicheSizes[i, "nicheDensity_GFDL"] <- density1[2]
  nicheSizes[i, "nicheDensity_MPI"] <- density2[2]
  nicheSizes[i, "nicheDensity_UKESM"] <- density3[2]
}

# Export results
write.csv(nicheSizes, paste0(outPath, "nicheSizes.csv"))