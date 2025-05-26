library(phylobase)
library(adephylo)
library(phytools)
library(geomorph)
library(geiger)
library(ggplot2)
library(reshape2)

# Set seed
set.seed(2025)

# Set path
relPath <- "E:/Research/PhylogeneticSignal/"
treePath <- "Data/Trees/"
speciesPath <- "E:/Research/PhylogeneticSignal/Data/Covariates/finalData.csv"
outPath <- "E:/Research/PhylogeneticSignal/Data/Signals/"
figPath <- "E:/Research/PhylogeneticSignal/Figures/"
outName <- "signalsThompson.csv"

# Read in tree
chosenTree <- readNewick(paste0(relPath, treePath, "Processed/ThompsonFinalNotGrafted.tre"))

# Attach traits to tree and prune
covariates <- read.csv(speciesPath, header = T)
covariates <- covariates[, c(1, 5, 6, 4, 3, 2, 7:12)]
covariates$species <- gsub(" ", "_", covariates$species)
row.names(covariates) <- covariates$species
tree <- as(chosenTree, "phylo")
remove <- name.check(tree, covariates)
covariates <- subset(covariates, !(species %in% remove$data_not_tree))
tree <- drop.tip(tree, tip = remove$tree_not_data)
phylotraits <- phylo4d(tree, covariates, label.type = "column", label.column = "species")

# Calculate phylogenetic signal (Moran and Cmean)
moran <- abouheif.moran(phylotraits, method = "Abouheif", nrepet = 9999)
Cmean <- abouheif.moran(phylotraits, method = "oriAbouheif", nrepet = 9999)
signals <- data.frame(trait = moran$names,
                      signalM = moran$obs,
                      pM = moran$pvalue,
                      signalC = Cmean$obs,
                      pC = Cmean$pvalue)

# Calculate phylogenetic signal (Pagel and Blomberg) and standardized effect sizes
effectSizes <- list()
counter <- 1
for (i in c(2:12)) {
  print(i)
  if (i == 4) next
  tmp <- as.matrix(covariates[, i])
  row.names(tmp) <- covariates$species
  Z <- physignal.z(A = tmp, phy = tree, lambda = "all", verbose = F, print.progress = F, iter = 9999)
  lambda <- phylosig(tree, tmp, method = "lambda", test = T, nsim = 9999)
  K <- phylosig(tree, tmp, method = "K", test = T, nsim = 9999)
  signals[i - 1, "signalL"] <- lambda$lambda
  signals[i - 1, "pL"] <- lambda$P
  signals[i - 1, "signalK"] <- K$K
  signals[i - 1, "pK"] <- K$P
  signals[i - 1, "signalZ"] <- Z$Z
  signals[i - 1, "pZ"] <- Z$pvalue
  effectSizes[[counter]] <- Z
  names(effectSizes)[[counter]] <- colnames(covariates)[i]
  counter <- counter + 1
}

# Adjust p-values for multiple comparisons and export signals
signals <- signals[c(1:2, 4:11),]
signals$pMadj <- p.adjust(signals$pM, method = "BH")
signals$pCadj <- p.adjust(signals$pC, method = "BH")
signals$pLadj <- p.adjust(signals$pL, method = "BH")
signals$pKadj <- p.adjust(signals$pK, method = "BH")
signals$pZadj <- p.adjust(signals$pZ, method = "BH")
signals <- signals[, c(1:3, 12, 4:5, 13, 6:7, 14, 8:9, 15, 10:11, 16)]
write.csv(signals, paste0(outPath, outName), row.names = F)

# Compare effect sizes and adjust for multiple comparisons
zcomps <- compare.physignal.z(effectSizes, two.tailed = F)
zpvals <- zcomps$pairwise.P
zpvals[row(zpvals) <= col(zpvals)] <- NA
lower_vals <- zpvals[lower.tri(zpvals)]
adjusted_vals <- p.adjust(lower_vals, method = "BH")
zpvals[lower.tri(zpvals)] <- adjusted_vals
zvals <- zcomps$pairwise.z
zvals[row(zvals) <= col(zvals)] <- NA

# Plot effect size comparisons

# Reshape and prepare data
zp_df <- melt(zpvals, varnames = c("Row", "Col"), value.name = "Adjusted_p")
zp_df$significant <- !is.na(zp_df$Adjusted_p) & zp_df$Adjusted_p < 0.05
zp_df$z <- melt(zvals)$value 
# Format labels
zp_df$label <- ifelse(
  is.na(zp_df$Adjusted_p), "",
  ifelse(zp_df$Adjusted_p < 0.001,
         paste0("P < 0.001\nZ = ", formatC(zp_df$z, format = "f", digits = 2)),
         paste0("P = ", formatC(zp_df$Adjusted_p, format = "f", digits = 3),
                "\nZ = ", formatC(zp_df$z, format = "f", digits = 2))))
xlabels <- c("Niche size",
             "Niche density",
             "Range size (0 km)",
             "Range size (100 km)",
             "Range change (0 km, SSP1)",
             "Range change (0 km, SSP3)",
             "Range change (0 km, SSP5)",
             "Range change (100 km, SSP1)",
             "Range change (100 km, SSP3)",
             "")
ylabels <- c("",
             "Niche density",
             "Range size (0 km)",
             "Range size (100 km)",
             "Range change (0 km, SSP1)",
             "Range change (0 km, SSP3)",
             "Range change (0 km, SSP5)",
             "Range change (100 km, SSP1)",
             "Range change (100 km, SSP3)",
             "Range change (100 km, SSP5)")
# Plot
plot <- ggplot(zp_df, aes(x = Col, y = Row)) +
  geom_tile(aes(fill = Adjusted_p), size = 0.8) +
  geom_text(aes(label = label), size = 2.5) +
  geom_text(data = subset(zp_df, significant),
            aes(label = "*"),
            vjust = -0.8, size = 4) +
  scale_fill_gradientn(colors = c("green", "yellow", "red"), na.value = "white", name = expression("Adjusted "*italic(P)*" value")) +
  scale_x_discrete(position = "top", labels = xlabels) +
  scale_y_discrete(labels = ylabels) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 37.5, hjust = 0, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "inside",
    legend.title = element_text(size = 8),
    legend.position.inside = c(0.8, 0.25),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(color = "white"),
    plot.margin = margin(l = 10, r = 40)
  )
ggsave(paste0(figPath, "Fig1.png"), plot, width = 180, height = 180, units = "mm")