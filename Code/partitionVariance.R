library(lme4)

# Set directory structure
relPath <- "E:/Research/PhylogeneticConservatism/"
covariatePath <- "Data/Covariates/finalData.csv"

# Read covariate data
traits <- read.csv(paste0(relPath, covariatePath), header = T)

# Extract genus
traits$genus <- unlist(lapply(strsplit(traits$species, " ", fixed = T), "[", 1))

# Perform mixed models
rangeSize100 <- lmer(rangeSize100 ~ 1 + (1 | genus), data = traits)
rangeSize0 <- lmer(rangeSize0 ~ 1 + (1 | genus), data = traits)
rangeChange100_SSP3 <- lmer(rangeChange100_SSP3 ~ 1 + (1 | genus), data = traits)
rangeChange100_SSP5 <- lmer(rangeChange100_SSP5 ~ 1 + (1 | genus), data = traits)
rangeChange100_SSP1 <- lmer(rangeChange100_SSP1 ~ 1 + (1 | genus), data = traits)
rangeChange0_SSP3 <- lmer(rangeChange0_SSP3 ~ 1 + (1 | genus), data = traits)
rangeChange0_SSP5 <- lmer(rangeChange0_SSP5 ~ 1 + (1 | genus), data = traits)
rangeChange0_SSP1 <- lmer(rangeChange0_SSP1 ~ 1 + (1 | genus), data = traits)
nicheSize <- lmer(nicheSize ~ 1 + (1 | genus), data = traits)
nicheDensity <- lmer(nicheDensity ~ 1 + (1 | genus), data = traits)

# Calculate variance explained by genus
round(as.data.frame(VarCorr(nicheSize))[1, "vcov"] / (as.data.frame(VarCorr(nicheSize))[1, "vcov"] + as.data.frame(VarCorr(nicheSize))[2, "vcov"]) * 100, digits = 1)
round(as.data.frame(VarCorr(nicheDensity))[1, "vcov"] / (as.data.frame(VarCorr(nicheDensity))[1, "vcov"] + as.data.frame(VarCorr(nicheDensity))[2, "vcov"]) * 100, digits = 1)
round(as.data.frame(VarCorr(rangeSize100))[1, "vcov"] / (as.data.frame(VarCorr(rangeSize100))[1, "vcov"] + as.data.frame(VarCorr(rangeSize100))[2, "vcov"]) * 100, digits = 1)
round(as.data.frame(VarCorr(rangeSize0))[1, "vcov"] / (as.data.frame(VarCorr(rangeSize0))[1, "vcov"] + as.data.frame(VarCorr(rangeSize0))[2, "vcov"]) * 100, digits = 1)
round(as.data.frame(VarCorr(rangeChange100_SSP3))[1, "vcov"] / (as.data.frame(VarCorr(rangeChange100_SSP3))[1, "vcov"] + as.data.frame(VarCorr(rangeChange100_SSP3))[2, "vcov"]) * 100, digits = 1)
round(as.data.frame(VarCorr(rangeChange100_SSP5))[1, "vcov"] / (as.data.frame(VarCorr(rangeChange100_SSP5))[1, "vcov"] + as.data.frame(VarCorr(rangeChange100_SSP5))[2, "vcov"]) * 100, digits = 1)
round(as.data.frame(VarCorr(rangeChange100_SSP1))[1, "vcov"] / (as.data.frame(VarCorr(rangeChange100_SSP1))[1, "vcov"] + as.data.frame(VarCorr(rangeChange100_SSP1))[2, "vcov"]) * 100, digits = 1)
round(as.data.frame(VarCorr(rangeChange0_SSP3))[1, "vcov"] / (as.data.frame(VarCorr(rangeChange0_SSP3))[1, "vcov"] + as.data.frame(VarCorr(rangeChange0_SSP3))[2, "vcov"]) * 100, digits = 1)
round(as.data.frame(VarCorr(rangeChange0_SSP5))[1, "vcov"] / (as.data.frame(VarCorr(rangeChange0_SSP5))[1, "vcov"] + as.data.frame(VarCorr(rangeChange0_SSP5))[2, "vcov"]) * 100, digits = 1)
round(as.data.frame(VarCorr(rangeChange0_SSP1))[1, "vcov"] / (as.data.frame(VarCorr(rangeChange0_SSP1))[1, "vcov"] + as.data.frame(VarCorr(rangeChange0_SSP1))[2, "vcov"]) * 100, digits = 1)