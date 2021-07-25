#!/usr/bin/env Rscript

# plotAFDist.R <mafs.gz file 1> <mafz.gz file 2> <output prefix>

library(ggplot2)

# parse arguments
args <- commandArgs(trailingOnly=TRUE)

maf1 = read.table(gzfile(args[1],'rt'), head=TRUE)
maf2 = read.table(gzfile(args[2],'rt'), head=TRUE)
outprefix = args[3]

# combine the data so that allele frequency densities can be plotted together
maf1$SNP_pval <- "1e-6"
maf2$SNP_pval <- "1e-2"
maf.comb <- rbind(maf1, maf2)
maf.comb$SNP_pval <- factor(maf.comb$SNP_pval, levels=c("1e-6","1e-2"), order=TRUE)

# plot
pdf(file=paste0(outprefix,".pdf"))
ggplot(maf.comb, aes(knownEM, fill = SNP_pval)) + geom_density(alpha = 0.4, bw=0.015) + theme_classic(base_size=16) + theme(axis.line = element_line(size=0.5)) + xlab("MAF") + ylab("Density")
invisible(dev.off())
