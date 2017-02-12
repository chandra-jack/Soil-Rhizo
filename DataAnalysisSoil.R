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
# 17. Does presence of other microbes effect plant growth?
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
summary(Mod_GLHT, test = adjusted("fdr"))
plot(print(confint(Mod_GLHT)))
mat
contrasts(TW_shoot$SingleFactor) = mat
TW_shoot$Weight2 <- log(TW_shoot$Weight)
mod1 <- aov(Weight2 ~ SingleFactor, data = TW_shoot)
summary(mod1, split = list(SingleFactor = list("Fl vs Pt" = 1, "Buffer vs WSM" = 2, "Soil vs None" = 3)), type = "III")

#----------
ckTb <- table(TW_shoot$SingleFactor)
MjK <-tapply(TW_shoot$Weight2, TW_shoot$SingleFactor, mean)
dfmy <- sum(ckTb)-14
SSE <- sum((TW_shoot$Weight2 - ave(TW_shoot$Weight2, TW_shoot$SingleFactor, FUN=mean))^2)
MSE     <- SSE / dfmy
(psiHat <- sum(MyContMat[1, ] * MjK))
lenSq <- sum(MyContMat[1, ]^2 / ckTb) 
(SE   <- sqrt(lenSq*MSE))
(tStat <- psiHat / SE)
(pVal <- 2 * (1-pt(abs(tStat), dfmy)))

Gen.Means <- summarise(group_by(TW_shoot, Inoculate), y.mean = mean(log(Weight)))

# This answers Questions 15b (No sig difference); 16b (Differs); 17b (differs)

TWS_2 <- ddply(TW_shoot, "Genotype", summarise, MEAN = mean(Weight), se = sqrt(var(Weight) / length(Weight)))
ggplot(TW_shoot, aes(x = Genotype, y = Weight)) + geom_boxplot() + labs(y = "Average shoot weight")+ theme(axis.text.x = element_text(size = 14))
ggsave(file = "ShootGenotypeContrast8Feb2017.pdf")

ggplot(TWS_2, aes(x = Genotype, y = MEAN))+ geom_bar(stat = "identity")+ geom_errorbar(aes(ymin = MEAN -se, ymax = MEAN + se), width = 0.5)


ggplot(TW_shoot, aes(x = Inoculate, y = Weight)) + geom_boxplot() + labs(y = "Average shoot weight")+ theme(axis.text.x = element_text(size = 14))
ggsave(file = "ShootInoculateContrast8Feb2017.pdf")

# m1 <- aov(Weight2 ~ Genotype*Inoculate*SoilConc, data = TW_shoot)
# m2 <- update(m1, ~. - Genotype:Inoculate:SoilLocation)
# m3 <- update(m2, ~. - Genotype:SoilLocation)
# m4 <- update(m3, ~. - Genotype:Inoculate)
# Anova(m1, type = "3")
# m1

# Setting up post hoc contrast for affect of soil
PostHocSoil <- read.csv("PostHocSoil.csv", header = F)
NoneVsLow <- PostHocSoil$V1
LowVsHigh <- PostHocSoil$V2
NoneVsHigh <- PostHocSoil$V3
mat.temp <- rbind(constant = 1/14, NoneVsLow, LowVsHigh, NoneVsHigh)
mat <- ginv(mat.temp)
mat <- mat[,-1]
NoneVsLow <- mat[,1]
LowVsHigh <- mat[,2]
NoneVsHigh <- mat[,3]
ContMatSoilC <- rbind("NoneVsLow" = NoneVsLow, "LowVsHigh" = LowVsHigh, "NoneVsHigh" = NoneVsHigh)

mod1 <- aov(Weight2 ~ SingleFactor, data = TW_shoot)
Mod_GLHT_Soil <- glht(mod1, linfct = mcp(SingleFactor = ContMatSoilC))
summary(Mod_GLHT_Soil,  test = adjusted("fdr"))

plot(print(confint(Mod_GLHT_Soil)))

SC <- ddply(TW_shoot, "SoilConc", summarise, MEAN = mean(Weight), se = sqrt(var(Weight) / length(Weight)))
ggplot(SC, aes(SoilConc, MEAN)) + geom_bar(stat = "identity")+ geom_errorbar(aes(ymin = MEAN -se, ymax = MEAN + se), width = 0.5)
ggsave(file = "SoilConc10Feb1017.pdf")


# Setting up posthoc contrast for soil location
PostHocLoc <- read.csv("PostHocSoilLoc.csv", header = F)
FlvsPT <- PostHocLoc$V1
MatchNoMatch <- PostHocLoc$V2
mat.temp <- rbind(constant = 1/14, FlvsPT, MatchNoMatch)
mat <- ginv(mat.temp)
mat <- mat[,-1]
FlvsPT <- mat[,1]
MatchNoMatch <- mat[,2]

ContMatSoilLoc <- rbind("FlvsPT" = FlvsPT, "MatchNoMatch" = MatchNoMatch)
mod1 <- aov(Weight2 ~ SingleFactor, data = TW_shoot)
Mod_GLHT_Soil_Loc <- glht(mod1, linfct = mcp(SingleFactor = ContMatSoilLoc))
summary(Mod_GLHT_Soil_Loc, test = adjusted("fdr"))

 SL <- ddply(TW_shoot, c("Genotype", "SoilLocation"), summarise, MEAN = mean(Weight), se = sqrt(var(Weight) / length(Weight)))
 SL <- subset(SL, SoilLocation != "none")
 
ggplot(SL, aes(SoilLocation, MEAN)) + geom_bar(stat = "identity")+ geom_errorbar(aes(ymin = MEAN -se, ymax = MEAN + se), width = 0.5) + facet_wrap(~Genotype)
ggsave(file = "SoilLocation11Feb2017.pdf")

# Setting up post hoc contrasts for interactions

IC <- read.csv("InteractionContrast.csv", header = F)
Fl_SoilLoc <- IC$V1
Fl_soil_ctrl <- IC$V2
PI_SoilLoc <- IC$V3
PI_soil_ctrl <- IC$V4
RhizY_St.Aug_PI <-IC$V5
RhizN_St.Aug_PI <-IC$V6
Rhiz_St.Aug_FL_PT <-IC$V7
Rhiz_PI_FL_PT <-IC$V8

mat.temp <- rbind(constant = 1/14, Fl_SoilLoc, Fl_soil_ctrl, PI_SoilLoc, PI_soil_ctrl, RhizY_St.Aug_PI, RhizN_St.Aug_PI, Rhiz_St.Aug_FL_PT, Rhiz_PI_FL_PT)
mat <- ginv(mat.temp)
mat <- mat[,-1]

Fl_SoilLoc <- mat[,1]
Fl_soil_ctrl <- mat[,2]
PI_SoilLoc <- mat[,3]
PI_soil_ctrl <- mat[,4]
RhizY_St.Aug_PI <-mat[,5]
RhizN_St.Aug_PI <-mat[,6]
Rhiz_St.Aug_FL_PT <-mat[,7]
Rhiz_PI_FL_PT <-mat[,8]

InteractionContrast <- rbind("Fl=SoilLoc" = Fl_SoilLoc, "Fl:soil=ctrl" = Fl_soil_ctrl, "Pt=SoilLoc" = PI_SoilLoc, "Pt:soil=ctrl" = PI_soil_ctrl, "RhizY::St.Aug=PI" = RhizY_St.Aug_PI, "RhizN::St.Aug=PI" = RhizN_St.Aug_PI, "Rhiz:StAug:Fl=PT" = Rhiz_St.Aug_FL_PT, "Rhiz:PI:Fl=PT" = Rhiz_PI_FL_PT)

mod1 <- aov(Weight2 ~ SingleFactor, data = TW_shoot)
Mod_GLHT_Inter <- glht(mod1, linfct = mcp(SingleFactor = InteractionContrast))
summary(Mod_GLHT_Inter, test = Ftest())
