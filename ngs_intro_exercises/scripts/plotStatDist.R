#!/usr/bin/env Rscript

# plotStatDist.R <file of vcf stats> <outfile prefix>

# parse arguments
args <- commandArgs(trailingOnly=TRUE)
df <- read.table(args[1],head=TRUE, na.strings=".")
suppressWarnings(df$EXCHET <- as.numeric(as.character(df$EXCHET))) # need this because R interprets this as a 'factor' due to many 'NA' at beginning
outprefix <- args[2]

# calculate percentiles
percdf=matrix(nrow=ncol(df)-2, ncol=13)
i=1
for (j in 3:ncol(df)) { percdf[i,] = quantile(df[,j], c(0.001, 0.01, seq(from=0.1, to=0.9, by=0.1), 0.99, 0.999), na.rm=TRUE); i=i+1 }
colnames(percdf) = names(quantile(df$DP, c(0.001, 0.01, seq(from=0.1, to=0.9, by=0.1), 0.99, 0.999), na.rm=TRUE))
rownames(percdf) = colnames(df[3:ncol(df)])

# print quantiles
cat("\n\nPercentiles\n\n")
print(percdf); cat("\n")

# plot stats
pdf(file=paste0(outprefix,".pdf"))
par(mfrow=c(2,2))
hist(df$DP, breaks=20, ylab="Number sites", xlab="Total site depth", main="DP", cex.lab=1.2)
hist(df$MQ, breaks=20, ylab="Number sites", xlab="Mapping quality", main="MQ", cex.lab=1.2)
hist(df$STRAND_BIAS, breaks=20, ylab="Number sites", xlab="Strand bias p-value", main="PV[0]", cex.lab=1.2)
hist(df$BASEQ_BIAS, breaks=20, ylab="Number sites", xlab="Base quality p-value", main="PV[1]", cex.lab=1.2)
hist(df$BASEQ_BIAS, breaks=10000, ylab="Number sites", xlab="Base quality p-value", main="PV[1] Zoom", cex.lab=1.2, xlim=c(0,0.001))
hist(df$MQ_BIAS, breaks=20, ylab="Number sites", xlab="Mapping quality bias p-value", main="PV[2]", cex.lab=1.2)
hist(df$POS_BIAS, breaks=20, ylab="Number sites", xlab="tail distance bias p-value", main="PV[3]", cex.lab=1.2)
hist(df$EXCHET, breaks=20, ylab="Number sites", xlab="Excess heterozygosity p-value", main="ExcHet", cex.lab=1.2)
invisible(dev.off())
