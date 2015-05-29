setwd("~/src/manderCap/")

# targetStats <- read.csv("targetMetrics.wholeHSPcounts.txt", sep="\t")
targetStatsAve <- read.csv("targetMetricsRENAMED.aveDepths.txt", sep="\t")
targetLengths <- read.csv("targetLengths.txt", sep="\t")
targetAveWithLength <- merge(targetStatsAve, targetLengths, "Target")

test <- read.csv("targetMetrics.aveDepths.sepCategories.txt", sep="\t")


SNPS_across_baits <- read.table("~/src/manderCap/SNPcountAcrossBaits.txt", quote="\"", header=TRUE)
#SNPpercentages_across_baits <- read.table("~/src/manderCap/percentagesOfSNPsAcrossBaits.txt", quote="\"")

targetAveWithLengthAndSNPcount <- merge(targetAveWithLength, SNPS_across_baits, "Target")

# Make the subsets based on depth cutoffs:
SNPtargetOver10 <- (targetAveWithLengthAndSNPcount[(targetAveWithLengthAndSNPcount$max100AverageTargetDepth >= 10),])
SNPtargetOver5 <- (targetAveWithLengthAndSNPcount[(targetAveWithLengthAndSNPcount$max100AverageTargetDepth >= 5),])
SNPtargetOver4 <- (targetAveWithLengthAndSNPcount[(targetAveWithLengthAndSNPcount$max100AverageTargetDepth >= 4),])
SNPtargetOver3 <- (targetAveWithLengthAndSNPcount[(targetAveWithLengthAndSNPcount$max100AverageTargetDepth >= 3),])






png("~/Box Sync/UCLA/Research/Papers/ambystoma_optimization/SNPS_per_target.png", width=600, height=600, units="px")
hist(SNPs_per_target$V1, breaks=100, xlab="Number of SNPs across entire target contig", ylab="Number of targets", main="# of SNPs for the different targets")
text(40, 400, "6,959 total targets with\nat least 1 qualifying SNP\n(homozygous in CTS, het in F1)")
dev.off()


# Percentage of sites within target regions that are potentially ancestry-informative
png("~/Box Sync/UCLA/Research/Papers/ambystoma_optimization/SNP_percentages_across_targets.png", width=600, height=600, units="px")
hist(SNPpercentages_across_baits$V1 * 100, breaks=100, xlab="Percentage of potentially-ancestry informative sites within target (%)", ylab="Number of targets", main="Percentage of sites within target regions that are\nhomozygous in CTS and heterozygous in F1")

text(4, 175, "min=0.2%, max=6.1%.\nmean=1.34%, sd=0.8%")
dev.off()


# How many SNPs within target regions?
png("~/Box Sync/UCLA/Research/Papers/ambystoma_optimization/SNPS_per_baitRegion.png", width=600, height=600, units="px")
hist(SNPS_across_baits$V1, breaks=100, xlab="Number of SNPs across bait region", ylab="Number of targets", main="Distribution of number of SNPs across designed\ntarget regions (no assembly overhang", axes="FALSE")
axis(1, at=seq(0, 25, by=1), las=2)
axis(2, at=seq(0, 2000, by=100) ,las=1)
text(10, 1300, "6,959 targets with at least\none qualifying SNP\n(homozygous in CTS and het in F1)")
dev.off()



png("~/Dropbox/UCLA/Research/Papers/ambystoma_optimization/targetAveDepths.png", width=600, height=600, units="px")
par(mar=c(5,5,5,5))
hist(targetStatsAve$AverageTargetDepthAcrossHSP, breaks=500, xlim=c(0,30), xlab="Average depth across target", ylab="Number of targets", main="Distribution of sequencing effort-corrected depth across targets", cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()


hist(targetStatsAve$max100AverageTargetDepth, breaks=500, xlim=c(0,30), xlab="Average depth across target", ylab="Number of targets", main="Distribution of sequencing effort-corrected depth across targets", cex.axis=1.5, cex.lab=1.5, cex.main=1.5)


median(targetStatsAve$max100AverageTargetDepth)
median(targetStatsAve$AverageTargetDepthAcrossHSP)


targetsBetterThan5AveCoverage <- (targetAveWithLength[(targetAveWithLength$AverageTargetDepthAcrossHSP > 5),])
targetsBetterThan10AveCoverage <- (targetAveWithLength[(targetAveWithLength$AverageTargetDepthAcrossHSP > 10),])
targetsBetterThan15AveCoverage <- (targetAveWithLength[(targetAveWithLength$AverageTargetDepthAcrossHSP > 15),])
targetsBetterThan20AveCoverage <- (targetAveWithLength[(targetAveWithLength$AverageTargetDepthAcrossHSP > 20),])
targetsBetterThan5MaxCoverage <- (targetAveWithLength[(targetAveWithLength$max100AverageTargetDepth > 5),])
targetsBetterThan10MaxCoverage <- (targetAveWithLength[(targetAveWithLength$max100AverageTargetDepth > 10),])
targetsBetterThan15MaxCoverage <- (targetAveWithLength[(targetAveWithLength$max100AverageTargetDepth > 15),])
targetsBetterThan20MaxCoverage <- (targetAveWithLength[(targetAveWithLength$max100AverageTargetDepth > 20),])
targetsBetterThan6MaxCoverage <- (targetAveWithLength[(targetAveWithLength$max100AverageTargetDepth > 6),])
targetsBetterThan7MaxCoverage <- (targetAveWithLength[(targetAveWithLength$max100AverageTargetDepth > 7),])
targetsBetterThan8MaxCoverage <- (targetAveWithLength[(targetAveWithLength$max100AverageTargetDepth > 8),])