library(ape)
library(geiger)
library(nlme)
library(visreg)
library(ggplot2)
library(ggpubr)
library(HDtweedie)

# Set path
treePath <- "E:/Research/PhylogeneticConservatism/Data/Trees/Processed/ThompsonFinalNotGrafted.tre"
figPath <- "E:/Research/PhylogeneticConservatism/Figures/"
modelPath <- "E:/Research/PhylogeneticConservatism/Data/ModelFits/modelFits.csv"

# Read in data
covariates <- read.csv("E:/Research/PhylogeneticConservatism/Data/Covariates/finalData.csv", header = T)

# Transform names
covariates[, "species"] <- gsub(" ", "_", covariates$species)
rownames(covariates) <- covariates[, "species"]

# Read phylogenetic tree and remove taxa not in tree
tree <- read.tree(treePath)
remove <- name.check(tree, covariates)
covariates <- subset(covariates, !(species %in% remove$data_not_tree))

# Calculate correlations
cor.test(log(covariates$nicheSize), log(covariates$nicheDensity))
cor.test(log(covariates$nicheSize), log(covariates$rangeSize0))
cor.test(log(covariates$nicheSize), log(covariates$rangeSize100))
cor.test(log(covariates$nicheDensity), log(covariates$rangeSize0))
cor.test(log(covariates$nicheDensity), log(covariates$rangeSize100))

# Perform regressions
lm1 <- lm(log(rangeSize100) ~ log(nicheSize) + log(nicheDensity) + log(nicheDensity):log(nicheSize), covariates)
summary(lm1)
confint(lm1)
visreg(lm1, "nicheSize", by = "nicheDensity", xtrans = log, breaks = 10)
lm2 <- lm(log(rangeSize0) ~ log(nicheSize) + log(nicheDensity) + log(nicheDensity):log(nicheSize), covariates)
summary(lm2)
confint(lm2)
visreg(lm2, "nicheSize", by = "nicheDensity", xtrans = log, breaks = 10)

# Check evolutionary model
tmp <- covariates[, "rangeSize100"]
names(tmp) <- rownames(covariates)
BM100 <- fitContinuous(tree, tmp, model = "BM")
OU100 <- fitContinuous(tree, tmp, model = "OU") # Lowest AIC
EB100 <- fitContinuous(tree, tmp, model = "EB")
LA100 <- fitContinuous(tree, tmp, model = "lambda")
WN100 <- fitContinuous(tree, tmp, model = "white")
tmp <- covariates[, "rangeSize0"]
names(tmp) <- rownames(covariates)
BM0 <- fitContinuous(tree, tmp, model = "BM")
OU0 <- fitContinuous(tree, tmp, model = "OU") # Lowest AIC
EB0 <- fitContinuous(tree, tmp, model = "EB")
LA0 <- fitContinuous(tree, tmp, model = "lambda")
WN0 <- fitContinuous(tree, tmp, model = "white")
tmp <- covariates[, "nicheSize"]
names(tmp) <- rownames(covariates)
BMs <- fitContinuous(tree, tmp, model = "BM")
OUs <- fitContinuous(tree, tmp, model = "OU") # Lowest AIC
EBs <- fitContinuous(tree, tmp, model = "EB")
LAs <- fitContinuous(tree, tmp, model = "lambda")
WNs <- fitContinuous(tree, tmp, model = "white")
tmp <- covariates[, "nicheDensity"]
names(tmp) <- rownames(covariates)
BMd <- fitContinuous(tree, tmp, model = "BM")
OUd <- fitContinuous(tree, tmp, model = "OU") # Lowest AIC
EBd <- fitContinuous(tree, tmp, model = "EB")
LAd <- fitContinuous(tree, tmp, model = "lambda")
WNd <- fitContinuous(tree, tmp, model = "white")
tmp <- covariates[, "rangeChange100_SSP1"]
names(tmp) <- rownames(covariates)
BMc100s1 <- fitContinuous(tree, tmp, model = "BM")
OUc100s1 <- fitContinuous(tree, tmp, model = "OU") # Lowest AIC
EBc100s1 <- fitContinuous(tree, tmp, model = "EB")
LAc100s1 <- fitContinuous(tree, tmp, model = "lambda")
WNc100s1 <- fitContinuous(tree, tmp, model = "white")
tmp <- covariates[, "rangeChange0_SSP1"]
names(tmp) <- rownames(covariates)
BMc0s1 <- fitContinuous(tree, tmp, model = "BM")
OUc0s1 <- fitContinuous(tree, tmp, model = "OU") # Lowest AIC
EBc0s1 <- fitContinuous(tree, tmp, model = "EB")
LAc0s1 <- fitContinuous(tree, tmp, model = "lambda")
WNc0s1 <- fitContinuous(tree, tmp, model = "white")
tmp <- covariates[, "rangeChange100_SSP3"]
names(tmp) <- rownames(covariates)
BMc100s3 <- fitContinuous(tree, tmp, model = "BM")
OUc100s3 <- fitContinuous(tree, tmp, model = "OU") # Lowest AIC
EBc100s3 <- fitContinuous(tree, tmp, model = "EB")
LAc100s3 <- fitContinuous(tree, tmp, model = "lambda")
WNc100s3 <- fitContinuous(tree, tmp, model = "white")
tmp <- covariates[, "rangeChange0_SSP3"]
names(tmp) <- rownames(covariates)
BMc0s3 <- fitContinuous(tree, tmp, model = "BM")
OUc0s3 <- fitContinuous(tree, tmp, model = "OU") # Lowest AIC
EBc0s3 <- fitContinuous(tree, tmp, model = "EB")
LAc0s3 <- fitContinuous(tree, tmp, model = "lambda")
WNc0s3 <- fitContinuous(tree, tmp, model = "white")
tmp <- covariates[, "rangeChange100_SSP5"]
names(tmp) <- rownames(covariates)
BMc100s5 <- fitContinuous(tree, tmp, model = "BM")
OUc100s5 <- fitContinuous(tree, tmp, model = "OU") # Lowest AIC
EBc100s5 <- fitContinuous(tree, tmp, model = "EB")
LAc100s5 <- fitContinuous(tree, tmp, model = "lambda")
WNc100s5 <- fitContinuous(tree, tmp, model = "white")
tmp <- covariates[, "rangeChange0_SSP5"]
names(tmp) <- rownames(covariates)
BMc0s5 <- fitContinuous(tree, tmp, model = "BM")
OUc0s5 <- fitContinuous(tree, tmp, model = "OU") # Lowest AIC
EBc0s5 <- fitContinuous(tree, tmp, model = "EB")
LAc0s5 <- fitContinuous(tree, tmp, model = "lambda")
WNc0s5 <- fitContinuous(tree, tmp, model = "white")
modelFits <- data.frame(NicheSize = c(BMs$opt$aicc, OUs$opt$aicc, EBs$opt$aicc, LAs$opt$aicc, WNs$opt$aicc),
           NicheDensity = c(BMd$opt$aicc, OUd$opt$aicc, EBd$opt$aicc, LAd$opt$aicc, WNd$opt$aicc),
           RangeSize100 = c(BM100$opt$aicc, OU100$opt$aicc, EB100$opt$aicc, LA100$opt$aicc, WN100$opt$aicc),
           RangeSize0 = c(BM0$opt$aicc, OU0$opt$aicc, EB0$opt$aicc, LA0$opt$aicc, WN0$opt$aicc),
           RangeChange100S1 = c(BMc100s1$opt$aicc, OUc100s1$opt$aicc, EBc100s1$opt$aicc, LAc100s1$opt$aicc, WNc100s1$opt$aicc),
           RangeChange100S3 = c(BMc100s3$opt$aicc, OUc100s3$opt$aicc, EBc100s3$opt$aicc, LAc100s3$opt$aicc, WNc100s3$opt$aicc),
           RangeChange100S5 = c(BMc100s5$opt$aicc, OUc100s5$opt$aicc, EBc100s5$opt$aicc, LAc100s5$opt$aicc, WNc100s5$opt$aicc),
           RangeChange0S1 = c(BMc0s1$opt$aicc, OUc0s1$opt$aicc, EBc0s1$opt$aicc, LAc0s1$opt$aicc, WNc0s1$opt$aicc),
           RangeChange0S3 = c(BMc0s3$opt$aicc, OUc0s3$opt$aicc, EBc0s3$opt$aicc, LAc0s3$opt$aicc, WNc0s3$opt$aicc),
           RangeChange0S5 = c(BMc0s5$opt$aicc, OUc0s5$opt$aicc, EBc0s5$opt$aicc, LAc0s5$opt$aicc, WNc0s5$opt$aicc))
row.names(modelFits) <- c("BM", "OU", "EB", "LA", "WN")
modelFits["minRow+3", ] <- sapply(modelFits, function(col) {
  min_val <- min(col, na.rm = TRUE)
  close_rows <- rownames(modelFits)[abs(col - min_val) <= 3]
  paste(close_rows, collapse = ", ")
})
write.csv(modelFits, modelPath, row.names = T)

# Perform regressions for range size
pgls100 <- gls(log(rangeSize100) ~ log(nicheSize) + log(nicheDensity) + log(nicheDensity):log(nicheSize),
               correlation = corMartins(value = 1, phy = tree, form = ~species, fixed = F),
               data = covariates,
               method = "ML")
pgls0 <- gls(log(rangeSize0) ~ log(nicheSize) + log(nicheDensity) + log(nicheDensity):log(nicheSize),
             correlation = corMartins(value = 1, phy = tree, form = ~species, fixed = F),
             data = covariates,
             method = "ML")
summary(pgls100)
confint(pgls100)
summary(pgls0)
confint(pgls0)

# Plot regressions
df <- visreg(pgls100, "nicheSize", plot = FALSE)$fit
plotA <- ggplot(df, aes(x = log(nicheSize), y = visregFit)) +
  geom_line(color = "black", size = 1.5) +
  labs(x = "log(niche size)", y = "log(range size (100-km buffer))")
df <- visreg(pgls100, "nicheDensity", plot = FALSE)$fit
plotB <- ggplot(df, aes(x = log(nicheDensity), y = visregFit)) +
  geom_line(color = "black", size = 1.5) +
  labs(x = "log(niche density)", y = "log(range size (100-km buffer))")
covariates$fit <- fitted(pgls100)
plotC <- ggplot(covariates, aes(x = log(nicheSize), y = log(nicheDensity), z = fit)) +
  geom_point(aes(color = fit)) +
  labs(x = "log(niche size)", y = "log(niche density)") +
  scale_color_gradient(low = "red", high = "green") +
  theme(legend.position = "inside", 
        legend.position.inside = c(0.3,0.875),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  labs(color = "log (range size (100-km buffer))")
plotD <- ggplot(covariates, aes(x = log(rangeSize100), y = fit)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, size = 1.5) +
  labs(x = "log(range size (100-km buffer)) (observed)", y = "log(range size (100-km buffer)) (predicted)") +
  theme(axis.title = element_text(size = 12))
ggsave(paste0(figPath, "Fig2.jpg"), ggarrange(plotC, plotD, ncol = 2, labels = c("a","b")), width = 180, height = 180, units = "mm")

df <- visreg(pgls0, "nicheSize", plot = FALSE)$fit
plotA <- ggplot(df, aes(x = log(nicheSize), y = visregFit)) +
  geom_line(color = "black", size = 1.5) +
  labs(x = "log(niche size)", y = "log(range size (0-km buffer))")
df <- visreg(pgls0, "nicheDensity", plot = FALSE)$fit
plotB <- ggplot(df, aes(x = log(nicheDensity), y = visregFit)) +
  geom_line(color = "black", size = 1.5) +
  labs(x = "log(niche density)", y = "log(range size (0-km buffer))")
covariates$fit <- fitted(pgls0)
plotC <- ggplot(covariates, aes(x = log(nicheSize), y = log(nicheDensity), z = fit)) +
  geom_point(aes(color = fit)) +
  labs(x = "log(niche size)", y = "log(niche density)") +
  scale_color_gradient(low = "red", high = "green") +
  theme(legend.position = "inside", 
        legend.position.inside = c(0.3,0.875),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  labs(color = "log (range size (0-km buffer))")
plotD <- ggplot(covariates, aes(x = log(rangeSize0), y = fit)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, size = 1.5) +
  labs(x = "log(range size (0-km buffer)) (observed)", y = "log(range size (0-km buffer)) (predicted)") +
  theme(axis.title = element_text(size = 12))
ggsave(paste0(figPath, "ExtendedDataFig1.jpg"), ggarrange(plotC, plotD, ncol = 2, labels = c("a","b")), width = 180, height = 180, units = "mm")

# Regression for range size change
coef(cv.HDtweedie(y = covariates$rangeChange100_SSP1, x = as.matrix(covariates[, c("nicheSize", "nicheDensity", "rangeSize100")]),
                  lambda = exp(seq(12, 0.5, length.out = 500)), nfolds = 10, alpha = 0, p = 1.36))
coef(cv.HDtweedie(y = covariates$rangeChange100_SSP3, x = as.matrix(covariates[, c("nicheSize", "nicheDensity", "rangeSize100")]),
                  lambda = exp(seq(12, 0.5, length.out = 500)), nfolds = 10, alpha = 0, p = 1.36))
coef(cv.HDtweedie(y = covariates$rangeChange100_SSP5, x = as.matrix(covariates[, c("nicheSize", "nicheDensity", "rangeSize100")]),
                  lambda = exp(seq(12, 0.5, length.out = 500)), nfolds = 10, alpha = 0, p = 1.36))
coef(cv.HDtweedie(y = covariates$rangeChange0_SSP1, x = as.matrix(covariates[, c("nicheSize", "nicheDensity", "rangeSize0")]),
                  lambda = exp(seq(12, 0.5, length.out = 500)), nfolds = 10, alpha = 0, p = 1.36))
coef(cv.HDtweedie(y = covariates$rangeChange0_SSP3, x = as.matrix(covariates[, c("nicheSize", "nicheDensity", "rangeSize0")]),
                  lambda = exp(seq(12, 0.5, length.out = 500)), nfolds = 10, alpha = 0, p = 1.36))
coef(cv.HDtweedie(y = covariates$rangeChange0_SSP5, x = as.matrix(covariates[, c("nicheSize", "nicheDensity", "rangeSize0")]),
                  lambda = exp(seq(12, 0.5, length.out = 500)), nfolds = 10, alpha = 0, p = 1.36))

coef(cv.HDtweedie(y = covariates$rangeChange100_SSP1, x = as.matrix(covariates[, c("nicheSize", "nicheDensity", "rangeSize100")]),
                  lambda = exp(seq(12, 0.5, length.out = 500)), nfolds = 10, alpha = 1, p = 1.36))
coef(cv.HDtweedie(y = covariates$rangeChange100_SSP3, x = as.matrix(covariates[, c("nicheSize", "nicheDensity", "rangeSize100")]),
                  lambda = exp(seq(12, 0.5, length.out = 500)), nfolds = 10, alpha = 1, p = 1.36))
coef(cv.HDtweedie(y = covariates$rangeChange100_SSP5, x = as.matrix(covariates[, c("nicheSize", "nicheDensity", "rangeSize100")]),
                  lambda = exp(seq(12, 0.5, length.out = 500)), nfolds = 10, alpha = 1, p = 1.36))
coef(cv.HDtweedie(y = covariates$rangeChange0_SSP1, x = as.matrix(covariates[, c("nicheSize", "nicheDensity", "rangeSize0")]),
                  lambda = exp(seq(12, 0.5, length.out = 500)), nfolds = 10, alpha = 1, p = 1.36))
coef(cv.HDtweedie(y = covariates$rangeChange0_SSP3, x = as.matrix(covariates[, c("nicheSize", "nicheDensity", "rangeSize0")]),
                  lambda = exp(seq(12, 0.5, length.out = 500)), nfolds = 10, alpha = 1, p = 1.36))
coef(cv.HDtweedie(y = covariates$rangeChange0_SSP5, x = as.matrix(covariates[, c("nicheSize", "nicheDensity", "rangeSize0")]),
                  lambda = exp(seq(12, 0.5, length.out = 500)), nfolds = 10, alpha = 1, p = 1.36))