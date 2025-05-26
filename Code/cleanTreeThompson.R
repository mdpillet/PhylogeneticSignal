library(phylobase)
library(TNRS)
library(ape)

# Set path
relPath <- "E:/Research/PhylogeneticSignal/"
treePath <- "Data/Trees/"
namePath <- "Data/Taxonomy/"

# Read in tree
tree <- readNewick(paste0(relPath, treePath, "Raw/Thompson.tre"))

# Standardize species names
species <- tree@label
occNames <- TNRS(taxonomic_names = species,
                 sources = "cact",
                 url = "http://vegbiendev.nceas.ucsb.edu:9975/tnrs_api.php")
write.csv(occNames, paste0(relPath, namePath, "rawNamesThompson.csv"), row.names = F)
taxonNames <- read.csv(paste0(relPath, namePath, "finalNamesThompson.csv"), header = T)
for (i in 1:length(species)) {
  tmp <- subset(taxonNames, Name_submitted == species[i])
  if (tmp$Name_override != "") tree@label[i] <- tmp$Name_override
  else tree@label[i] <- tmp$Accepted_species
}

# Drop tips for which no taxonomic opinion can be rendered (Opuntia pusilla and Mammillaria lloydii)
nTips(tree)
prunedTree <- prune(tree, tips.exclude = c("Opuntia pusilla", "Mammillaria lloydii"))
nTips(prunedTree)

# Drop duplicate tips (first tip index is retained)
dupTips <- duplicated(tipLabels(prunedTree))
dupNames <- tipLabels(prunedTree)[dupTips]
droppedTips <- NULL
for (i in 1:length(dupNames)) {
 tmpTips <- as.integer(names(tipLabels(prunedTree) == dupNames[i])[which(tipLabels(prunedTree) == dupNames[i])])
 tmpTips <- tmpTips[2:length(tmpTips)]
 droppedTips <- c(droppedTips, tmpTips)
}
droppedTips <- unique(droppedTips)
prunedTree <- prune(prunedTree, tips.exclude = droppedTips)
nTips(prunedTree)

# Export cleaned tree
tmpTree <- as(prunedTree, "phylo")
write.tree(tmpTree, paste0(relPath, treePath, "Processed/ThompsonPruned.tre"))