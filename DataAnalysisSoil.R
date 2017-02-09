# Experiment setup
#-----------------
# Two genotypes (one from Portugal and one from Florida)
# Two types of soil (Portugal and Florida)
# Two soil concentrations
# Presence or absence of rhizobia

# Data collected
#---------------
# Final weight above-ground mass
# Final weight below-ground mass
# Presence or absence of nodules
# Weight of nodules, if present
# Plants establish or fail

# Questions to answer with the Soil Rhizobia data
#------------------------------------------------
# 1. Do invasive (Fl) plants grow better in co-evolved or un-coevolved soil?
# 2. Do native (PT) plants grow better in co-evolved or un-coevolved soil?
       # a) Total plant mass
       # b) Above ground bio mass
# 3. Relative to native plants, do invasive plants grow better in co-evolved      or un-coevolved soil?
       # a) Total plant mass
       # b) Above ground bio mass
# 4. Do invasive plants change their root:shoot allocation in un-coevolved         soil?
# 5. Do native plants change their root:shoot allocation in un-coevolved           soil?
# 6. Do invasive plants grow better with rhizobia?
       # a) Total plant mass
       # b) Above ground bio mass
# 7. Do native plants grow better with rhizobia?
       # a) Total plant mass
       # b) Above ground bio mass
# 8. Does the presence of rhizobia increase the probability of success in           invasive plants?
# 9. Does the presence of rhizobia increase the probability of success in          native plants?
# 10. Does the concentration of other microbes in soil affect:
       # a) plant growth
            # i) Total plant mass
            # ii) Above ground biomass
            # iii) Below ground biomass
      # b) Plant success
      # c) Nodule formation
# 11. Do plants form more nodules with coevolved or non-coevolved soil?
# 12. Do invasive plants form more nodules compared to invasive plants?
# 13. Does nodule formation depend on what else is present in the soil?
# 14. Putting it all together- are there any interactions?
      # a) Total plant mass
      # b) Above ground bio mass
      # c) Below ground bio mass
      # d) Success establishing
      # e) Nodule formation
      # f) Nodule number
# 15. Overall do invasive plants grow better than native plants?
      # a) Total biomass
      # b) Above ground
      # c) Below ground
#16. Does rhizobia affect plant growth?
      # a) Total biomass
      # b) Above ground
      # c) Below ground    


# Library Packages
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(car)
library(multcomp)

# Data Analysis
setwd("~/Documents/Friesen lab/Soil_Rhizobia")
#-------------------------------------

# Starting with the TissueWeight.csv file. Seeing if can do a missing cell model, as mentioned in Logan Ch.12.

TW <- read.csv("TissueWeight.csv")
TW_shoot <- subset(TW, TissueType == "Shoot")
boxplot(log(Weight) ~ Genotype * Inoculate * SoilConc * SoilLocation, TW_shoot)
replications(log(Weight) ~ Genotype * Inoculate * SoilConc * SoilLocation, TW_shoot)
# This shows that the data is unbalanced and missing cells
TW_shoot$SingleFactor <- as.factor(paste(TW_shoot$Genotype, TW_shoot$Inoculate, TW_shoot$SoilConc, TW_shoot$SoilLocation, sep = ""))
unique(TW_shoot$SingleFactor)
TW_shoot %>% group_by(SingleFactor) %>% summarise(no_rows = length(SingleFactor))

# First test: All possible outcomes- Are there any differences?
mod1 <- aov(log(Weight) ~ SingleFactor, data = TW_shoot)
summary(glht(mod1, linfct = mcp(SingleFactor = "Tukey"), test=adjusted("Shaffer")))
summary(mod1)

Mod1_cld <- cld(glht(mod1, linfct = mcp(SingleFactor = "Tukey")))

TW_shoot$CLD <- ifelse(TW_shoot$SingleFactor == "PI493292BufferlowFL", "a", ifelse(TW_shoot$SingleFactor == "PI493292BufferlowPT", "ac", ifelse(TW_shoot$SingleFactor == "PI493292Buffernonenone", "ab", ifelse(TW_shoot$SingleFactor == "StAug2BufferlowPT", "ac", ifelse(TW_shoot$SingleFactor == "StAug2BufferlowFL", "a", ifelse(TW_shoot$SingleFactor == "StAug2Buffernonenone", "a", ifelse(TW_shoot$SingleFactor == "StAug2WSMnonenone", "bc", "c")))))))

ggplot(TW_shoot, aes(x = Genotype, y = Weight, colour = Inoculate, label = CLD)) + stat_summary(fun.data = "mean_se", na.rm = T) + facet_grid(SoilConc ~ SoilLocation) + geom_text(position = "dodge")

TWS <- ddply(TW_shoot, "SingleFactor", summarise, MEAN = mean(Weight), se = sqrt(var(Weight) / length(Weight)))
PlantLetters <- cld(glht(mod1, linfct = mcp(SingleFactor = "Tukey")))$mcletters
P_lett <- data.frame(SingleFactor = names(PlantLetters$Letters), LeT =PlantLetters$Letters)
TWS <- merge(TWS, P_lett, by ="SingleFactor")

TWS$TextPlace <- ifelse(!is.na(TWS$se), TWS$MEAN + TWS$se + .03, TWS$MEAN + .03)
ggplot(TWS, aes(x = SingleFactor, y = MEAN)) + geom_errorbar(aes(ymin = MEAN -se, ymax = MEAN + se)) + geom_bar(stat = "identity", position = "dodge") + geom_text(aes(SingleFactor, TextPlace, label = LeT), vjust = 0) + theme(axis.text.x = element_text(angle = 90)) + labs(x = "Treatments combined", y = "Average shoot weight", title = "Significant differences between all treatments: Shoots 8 Feb 2017")
ggsave(file = "AllCompShoots8Feb2017.pdf")

#---------------------------------
# Setting up contrast matrix for the planned contrasts

PCcsv <- read.csv("PlannedContrasts.csv", header = F) # The contrasts were setup following this website: http://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html -However, because the matrix is not square, cannot use solve- had to use ginv instead. Then, because glht requires constrasts to be by the row instead of the column, had to extract each column of the matrix and bind as rows.
FlvsPt <- PCcsv$V1
BufvsWSM <- PCcsv$V2
Soil <- PCcsv$V3
mat.temp <- rbind(constant = 1/14, FlvsPt, BufvsWSM, Soil)
mat <- ginv(mat.temp)
mat <- mat[,-1]
FlvsPt <- mat[,1]
BufvsWSM <- mat[,2]
Soil <- mat[,3]
MyContMat <- rbind("FlvsPt" = FlvsPt, "BufvsWSM" = BufvsWSM, "Soil" = Soil)
# ---------------------------------------------
Mod_GLHT <- glht(mod1, linfct = mcp(SingleFactor = MyContMat))
summary(Mod_GLHT, test = Ftest())
plot(print(confint(Mod_GLHT)))
mat
contrasts(TW_shoot$SingleFactor) = mat
mod1 <- aov(log(Weight) ~ SingleFactor, data = TW_shoot)
summary.aov(mod1, split = list(SingleFactor = list("Fl vs Pt" = 1, "Buffer vs WSM" = 2, "Soil vs None" = 3)), type = "marginal")
p.adjust(c(5.19e-08, 0.00253,1.50e-08,2.86e-05), "bonferroni")

ckTb <- table(TW_shoot$SingleFactor)
dfmy <- sum(ckTb)-14

# This answers Questions 15b (No sig difference); 16a
