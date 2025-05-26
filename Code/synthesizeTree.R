library(rtrees)
library(ape)

# Set path
relPath <- "E:/Research/PhylogeneticSignal/"
treePath <- "Data/Trees/"
speciesPath <- "E:/Research/DroughtForecasts/Manuscript/Tables/SupplementaryTable2.csv"

# Read tree
megaTree <- read.tree(paste0(relPath, treePath, "Processed/ThompsonPruned.tre"))

# Construct tree with grafted tips
species <- unique(read.csv(speciesPath, header = T)$species)
speciesList <- sp_list_df(species, "plant")
speciesList <- subset(speciesList, genus != "Cactus") # Get rid of Cactus coquimbanus
speciesList <- speciesList[!grepl("_x_", speciesList$species, fixed = T),] # Get rid of hybrids
speciesList$close_genus <- NA
speciesList[speciesList$genus == "Borzicactus", "close_genus"] <- "Matucana"
speciesList[speciesList$species == "Borzicactus_acanthurus", "close_genus"] <- "Haageocereus"

graftedTree <- get_tree(sp_list = speciesList,
                        tree = megaTree,
                        taxon = "plant",
                        scenario = "at_basal_node",
                        show_grafted = F,
                        tree_by_user = T)
# plot(graftedTree, type = "fan", cex = 0.35)
write.tree(graftedTree, paste0(relPath, treePath, "Processed/ThompsonFinalGrafted.tre"))

# Construct tree without grafted tips
graftedTips <- subset(graftedTree$graft_status, status == "exisiting species in the megatree")$species
speciesList <- speciesList[speciesList$species %in% graftedTips, ]
ungraftedTree <- get_tree(sp_list = speciesList,
                        tree = megaTree,
                        taxon = "plant",
                        scenario = "at_basal_node",
                        show_grafted = F,
                        tree_by_user = T)
plot(ungraftedTree, type = "fan", cex = 0.5)
write.tree(ungraftedTree, paste0(relPath, treePath, "Processed/ThompsonFinalNotGrafted.tre"))