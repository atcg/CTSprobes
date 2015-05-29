setwd("~/src/manderCap/")
mappingStats <- read.csv("mappingStats.csv", stringsAsFactors=FALSE)
mappingStatsNoDSN <- mappingStats[-c(18,21),]
mappingStatsNoDSN$post.50_uL_PCR_totalDNA <- mappingStatsNoDSN$post.50_uL_PCR_concentration_cleaned * 50
mappingStatsNoDSN$CTSonly6iterRawMapRatePERCENT <- mappingStatsNoDSN$CTSonly6iterRawMapRate * 100
mappingStatsNoDSN$CTSonly6iterDiscountMapRatePERCENT <- mappingStatsNoDSN$CTSonly6iterDiscountMapRate * 100
mappingStatsNoDSN$CTSonly6iterPercDupesPERCENT <- mappingStatsNoDSN$CTSonly6iterPercDupes * 100



png("~/Dropbox/UCLA/Research/Papers/ambystoma_optimization/PCRconcentration_to_rawMappingRate.png", width=1000, height=500, units="px")
par(mfrow=c(1,2))
plot(mappingStatsNoDSN$CTSonly6iterRawMapRatePERCENT ~ mappingStatsNoDSN$post.50_uL_PCR_totalDNA, col=mappingStatsNoDSN$individual, pch=20, xlab="Total DNA after amplifying 1/2 of the post-enrichment pool (ng)", ylab="Raw mapping rate (%)", main="All Libraries", cex=1.2)
abline(lm(mappingStatsNoDSN$CTSonly6iterRawMapRatePERCENT ~ mappingStatsNoDSN$post.50_uL_PCR_totalDNA), lwd=2, col="red")
plot(mappingStatsNoDSN$CTSonly6iterRawMapRatePERCENT[-18] ~ mappingStatsNoDSN$post.50_uL_PCR_totalDNA[-18], col=mappingStatsNoDSN$individual, pch=20, xlab="Total DNA after amplifying 1/2 of the post-enrichment pool (ng)", ylab="Raw mapping rate (%)", main="Removing Outlier", cex=1.2)
abline(lm(mappingStatsNoDSN$CTSonly6iterRawMapRatePERCENT[-18] ~ mappingStatsNoDSN$post.50_uL_PCR_totalDNA[-18]), col="red", lwd=2)
dev.off()

png("~/Dropbox/UCLA/Research/Papers/ambystoma_optimization/qPCRave_to_rawMappingRate.png", width=600, height=600, units="px")
plot(mappingStatsNoDSN$CTSonly6iterRawMapRatePERCENT ~ mappingStatsNoDSN$Ave_delta_Cp, col=mappingStatsNoDSN$individual, pch=20, xlab="Mean change in qPCR cycle number after enrichment", ylab="Raw mapping rate (%)", main="Predicting enrichment efficiency from qPCR validation", cex=1.2)
abline(lm(mappingStatsNoDSN$CTSonly6iterRawMapRatePERCENT ~ mappingStatsNoDSN$Ave_delta_Cp), lwd=2, col="red")
dev.off()


# Figure 5--4 panel showing reduction in PCR duplication rates with different variables and increase in unique read on target rate
png("~/Dropbox/UCLA/Research/Papers/ambystoma_optimization/four-panel.png", width=1200, height=1200, units="px")
par(mfrow=c(2,2), mar=c(5,5,5,5))

    # PCR duplication rates
plot(mappingStatsNoDSN$individual_dna_in_capture, mappingStatsNoDSN$CTSonly6iterPercDupesPERCENT, col=c(mappingStatsNoDSN$individual), cex=2, pch=20, xlab="Amount of individual DNA in capture reaction (ng)", ylab="PCR duplication rate (%)", main="Input DNA vs. PCR duplication rate", cex.lab=1.8, cex.main=1.8, cex.axis=1.8)
abline(lm(mappingStatsNoDSN$CTSonly6iterPercDupesPERCENT ~ mappingStatsNoDSN$individual_dna_in_capture), lwd=2, col="red")
plot(mappingStatsNoDSN$cot.multiplier, mappingStatsNoDSN$CTSonly6iterPercDupesPERCENT, col=mappingStatsNoDSN$individual, cex=2, pch=20, xlab="Amount of c0t-1 in capture reaction (multiplier)", ylab="PCR duplication rate (%)", main="Amount of c0t-1 vs. PCR duplication rate", cex.lab=1.8, cex.main=1.8, cex.axis=1.8)
abline(lm(mappingStatsNoDSN$CTSonly6iterPercDupesPERCENT ~ mappingStatsNoDSN$cot.multiplier), lwd=2, col="red")

      # Discounted map rates
plot(mappingStatsNoDSN$individual_dna_in_capture, mappingStatsNoDSN$CTSonly6iterDiscountMapRatePERCENT, col=mappingStatsNoDSN$individual, cex=2, pch=20, xlab="Amount of individual DNA in capture reaction (ng)", ylab="Unique reads on target (%)", main="Input DNA vs. unique reads on target", cex.lab=1.8, cex.main=1.8, cex.axis=1.8)
abline(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRatePERCENT ~ mappingStatsNoDSN$individual_dna_in_capture), lwd=2, col="red")
plot(mappingStatsNoDSN$cot.multiplier, mappingStatsNoDSN$CTSonly6iterDiscountMapRatePERCENT, col=mappingStatsNoDSN$individual, cex=2, pch=20, xlab="Amount of c0t-1 in capture reaction (multiplier)", ylab="Unique reads on target (%)", main="Amount of c0t-1 vs. unique reads on target", cex.lab=1.8, cex.main=1.8, cex.axis=1.8)
abline(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRatePERCENT ~ mappingStatsNoDSN$cot.multiplier), lwd=2, col="red")

dev.off()


# Figure showing average target depth


# Figure showing predicted vs. expected values for two-variable model
png("~/Dropbox/UCLA/Research/Papers/ambystoma_optimization/predictedVsRealValues.png", width=600, height=600, units="px")
par(mar=c(5,5,5,5))
twoVarModel <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRatePERCENT ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual_dna_in_capture)
plot(twoVarModel$fitted.values, mappingStatsNoDSN$CTSonly6iterDiscountMapRatePERCENT, cex=1.5, pch=20, xlab="Predicted unique reads on target from c0t-1 + individual input DNA model (%)", ylab="Actual unique reads on target (%)", main="Predicted vs. actual percentages of unique reads on target", cex.lab=1.4, cex.main=1.4, cex.axis=1.4, col=mappingStatsNoDSN$individual)
abline(coef=c(0,1), lwd=3, col="red")
dev.off()

summary(twoVarModel)

# Table 3: Model comparison for discounted mapping rate
  # Input DNA only
summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture))
AIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture))
BIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture))
  # cot1 only
summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier))
AIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier))
BIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier))
  # individual only
summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual))
AIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual))
BIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual))
  # Input DNA and cot-1
summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier))
AIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier))
BIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier))
  # cot-1 and individual ID
summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual))
AIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual))
BIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual))
  # Input DNA and individual ID
summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$individual))
AIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$individual))
BIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$individual))
  # Input DNA and cot-1 and individual ID
summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual))
AIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual))
BIC(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual))



# Table 4: Model comparison for depth across HSP
# Input DNA only
summary(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual_dna_in_capture))
AIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual_dna_in_capture))
BIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual_dna_in_capture))
# cot1 only
summary(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$cot.multiplier))
AIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$cot.multiplier))
BIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$cot.multiplier))
# individual only
summary(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual))
AIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual))
BIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual))
# Input DNA and cot-1
summary(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier))
AIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier))
BIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier))
# cot-1 and individual ID
summary(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual))
AIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual))
BIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual))
# Input DNA and individual ID
summary(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$individual))
AIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$individual))
BIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$individual))
# Input DNA and cot-1 and individual ID
summary(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual))
AIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual))
BIC(lm(mappingStatsNoDSN$AverageTargetDepthAcrossHSP ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual))
















# Generate the nested models
cot1Mod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier)
inputDNAmod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture)
individualMod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual)
cot1AndInputDNAmod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual_dna_in_capture)
cot1AndIndividualMod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual)
inputDNAandIndividualMod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$individual)
all3Mod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual)


cot1AndInputDNAmod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual_dna_in_capture + G(mappingStatsNoDSN$individual))



totalDNAmod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$total_dna_in_capture)
all4Mod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual + mappingStatsNoDSN$total_dna_in_capture)
all3v2Mod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$total_dna_in_capture) 



# Are CTS/BTS/F1 distinguishable?
plot(mappingStatsNoDSN$condition, mappingStatsNoDSN$CTSonly6iterDiscountMapRate, col=mappingStatsNoDSN$individual, pch=20, ylab="Discounted capture rate", xlab="Reaction condition number", main="CTS-only assembly after six ARC iterations--looking for differences between CTS/BTS/F1")


# Looks at predictive power of our variables
plot(mappingStatsNoDSN$individual_dna_in_capture, mappingStatsNoDSN$CTSonly6iterDiscountMapRate, col=mappingStatsNoDSN$individual, pch=20, cex=1, ylab="Discounted capture rate", xlab="Individual DNA in capture", main="CTS-only assembly after six ARC iterations")
plot(mappingStatsNoDSN$individual_dna_in_capture, mappingStatsNoDSN$CTSonly6iterRawMapRate, col=mappingStatsNoDSN$individual, pch=20, cex=1, ylab="Discounted capture rate", xlab="Individual DNA in capture", main="RAW MAP")
plot(mappingStatsNoDSN$individual_dna_in_capture, mappingStatsNoDSN$CTSonly6iterPercDupes, col=mappingStatsNoDSN$individual, pch=20, cex=1, ylab="Discounted capture rate", xlab="Individual DNA in capture", main="CTS-only assembly after six ARC iterations")


indDNAmod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture)
summary(indDNAmod)
abline(indDNAmod, col="red", lwd=2)
indDNAmodResiduals <- resid(indDNAmod)
indDNAmodInd<- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$individual)
summary(indDNAmodInd)


all3mod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual)
summary(all3mod)
anova(all3mod)

bothMod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier)
plot(bothMod)
plot(fitted(bothMod), mappingStatsNoDSN$CTSonly6iterDiscountMapRate)
abline(coef=c(0,1), col="red", lwd=2)

# Plot of the residuals from the individual DNA model
plot(mappingStatsNoDSN$individual_dna_in_capture, indDNAmodResiduals, col=mappingStatsNoDSN$individual, pch=20, cex=1, ylab="Residual", xlab="Individual DNA in capture", main="CTS-only assembly after six ARC iterations--residuals from individual input DNA model")




# cot-1 Analyses
plot(mappingStatsNoDSN$cot.multiplier, mappingStatsNoDSN$CTSonly6iterDiscountMapRate, col=mappingStatsNoDSN$individual, pch=20, cex=1, ylab="Discounted capture rate", xlab="cot-1 multiplier", main="CTS-only assembly after six ARC iterations")
cot1mod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier)
summary(cot1mod)
abline(cot1mod, col="red", lwd=2)
cot1modResiduals <- resid(cot1mod)

# Plot of the residuals from the cot-1 model
plot(mappingStatsNoDSN$individual_dna_in_capture, cot1modResiduals, col=mappingStatsNoDSN$individual, pch=20, cex=1, ylab="Residual", xlab="Individual DNA in capture", main="CTS-only assembly after six ARC iterations--residuals from cot-1 model")



# Predict individual DNA residuals based on cot-1
plot(mappingStatsNoDSN$cot.multiplier, indDNAmodResiduals, col=mappingStatsNoDSN$individual, pch=20)
summary(lm(indDNAmodResiduals ~ mappingStatsNoDSN$cot.multiplier))
abline(lm(indDNAmodResiduals ~ mappingStatsNoDSN$cot.multiplier), col="red", lwd=2)

# Predict cot-1 residuals based on input DNA
plot(mappingStatsNoDSN$individual_dna_in_capture, cot1modResiduals, col=mappingStatsNoDSN$individual, pch=20)
summary(lm(cot1modResiduals ~ mappingStatsNoDSN$individual_dna_in_capture))
abline(lm(cot1modResiduals ~ mappingStatsNoDSN$individual_dna_in_capture), col="red", lwd=2)




summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier))
summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier))
summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture))
