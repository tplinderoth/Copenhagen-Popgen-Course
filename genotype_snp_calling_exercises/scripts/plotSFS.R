#!/usr/bin/env Rscript

# plotSFS.R <SFS file> <output prefix>

# parse inputs
args <- commandArgs(trailingOnly=TRUE)

sfs <- scan(args[1])
outprefix <- args[2]
n <- (length(sfs)-1)/2

sfs = sfs[2:(n+1)] # restrict to folded categories and remove fixed class

# plot
pdf(file=paste0(outprefix,".pdf"),width=14,height=7)
barplot(sfs, xlab="MAF", ylab="Number SNPs", names=1:n, cex.names=0.8, cex.axis=1.2, cex.lab=1.2)
invisible(dev.off())
