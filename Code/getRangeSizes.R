# Set path
speciesPath <- "E:/Research/DroughtForecasts/Outputs/"
outPath <- "E:/Research/PhylogeneticSignal/Data/RangeSizes/rangeSizes.csv"

# List directories 
speciesDirs <- list.dirs(speciesPath, full.names = F, recursive = F)

# Extract range sizes
rangeSizes <- data.frame(species = speciesDirs)
for (i in 1:nrow(rangeSizes)) {
  if (file.exists(paste0(speciesPath, speciesDirs[i], "/RangeSizesBestModel/rangeSizes.csv"))) {
    tmp <- read.csv(paste0(speciesPath, speciesDirs[i], "/RangeSizesBestModel/rangeSizes.csv"))
    rangeSizes[i, "rangeSize100"] <- tmp[1, "RangeSize"]
    tmp <- read.csv(paste0(speciesPath, speciesDirs[i], "/RangeSizesBestModel/rangeSizesNoDisp.csv"))
    rangeSizes[i, "rangeSize0"] <- tmp[1, "RangeSize"]
    rangeSizes[i, "GCM"] <- tmp[1, "GeneratingModel"]
  }
}
rangeSizes <- rangeSizes[complete.cases(rangeSizes),]
write.csv(rangeSizes, outPath, row.names = F)