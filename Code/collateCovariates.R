# Set path structure
nicheData <- "E:/Research/PhylogeneticSignal/Data/NicheSizes/nicheSizes.csv"
geoData <- "E:/Research/PhylogeneticSignal/Data/RangeSizes/rangeSizes.csv"
speciesPath <- "E:/Research/DroughtForecasts/Manuscript/Tables/SupplementaryTable2.csv"
rangeSizePath <- "Data/RangeSizes/rangeSizes.csv"
outPath <- "E:/Research/PhylogeneticSignal/Data/Covariates/finalData.csv"

# Read data
niches <- read.csv(nicheData)[, 2:8]
rangeSizes <- read.csv(geoData)
rangeChanges <- read.csv(speciesPath)

# Remove Leucostele skottsbergii and Lobivia cinnabarina
niches <- subset(niches, species != "Leucostele skottsbergii" & species != "Lobivia cinnabarina")
rangeSizes <- subset(rangeSizes, species != "Leucostele skottsbergii" & species != "Lobivia cinnabarina")

# Compile data
tmp <- rangeSizes
for (i in 1:nrow(tmp)) {
  print(i)
  GCM <- tmp[i, "GCM"]
  if (GCM == "MPI-ESM1-2-HR") {
    tmp[i, "nicheSize"] <- niches[i, "nicheSize_MPI"]
    tmp[i, "nicheDensity"] <- niches[i, "nicheDensity_MPI"]
  }
  if (GCM == "UKESM1-0-LL") {
    tmp[i, "nicheSize"] <- niches[i, "nicheSize_UKESM"]
    tmp[i, "nicheDensity"] <- niches[i, "nicheDensity_UKESM"]
  }
  if (GCM == "GFDL-ESM4") {
    tmp[i, "nicheSize"] <- niches[i, "nicheSize_GFDL"]
    tmp[i, "nicheDensity"] <- niches[i, "nicheDensity_GFDL"]
  }
  tmp[i, "rangeChange0_SSP1"] <- subset(rangeChanges, species == tmp[i, "species"] & bestModel == GCM & Dispersal == "Dispersal: none" & SSP == "SSP1-2.6")$RangeChange
  tmp[i, "rangeChange0_SSP3"] <- subset(rangeChanges, species == tmp[i, "species"] & bestModel == GCM & Dispersal == "Dispersal: none" & SSP == "SSP3-7.0")$RangeChange
  tmp[i, "rangeChange0_SSP5"] <- subset(rangeChanges, species == tmp[i, "species"] & bestModel == GCM & Dispersal == "Dispersal: none" & SSP == "SSP5-8.5")$RangeChange
  tmp[i, "rangeChange100_SSP1"] <- subset(rangeChanges, species == tmp[i, "species"] & bestModel == GCM & Dispersal == "Dispersal: 100 km" & SSP == "SSP1-2.6")$RangeChange
  tmp[i, "rangeChange100_SSP3"] <- subset(rangeChanges, species == tmp[i, "species"] & bestModel == GCM & Dispersal == "Dispersal: 100 km" & SSP == "SSP3-7.0")$RangeChange
  tmp[i, "rangeChange100_SSP5"] <- subset(rangeChanges, species == tmp[i, "species"] & bestModel == GCM & Dispersal == "Dispersal: 100 km" & SSP == "SSP5-8.5")$RangeChange
}

# Export data
write.csv(tmp, outPath, row.names = F)