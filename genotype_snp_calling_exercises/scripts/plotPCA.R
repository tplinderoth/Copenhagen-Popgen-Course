#!/usr/bin/env Rscript

# plotAFDist.R <covariance file> <metadata file> <output prefix>

# parse inputs
args <- commandArgs(trailingOnly=TRUE)

covmat <- as.matrix(read.table(args[1], head=FALSE))
meta <- read.table(args[2], head=TRUE)
outprefix <- args[3]

# set up ecomorph colors for plotting
morphcol <- replace(as.character(meta$ecomorph.phenotype), which(meta$ecomorph.phenotype == "Littoral"), "yellow4")
morphcol <- replace(morphcol, which(meta$ecomorph.phenotype == "Benthic"), "midnightblue")

# perform eigen decomposition on the covariance matrix
eig <- eigen(covmat)

# calculate the percent genetic variance explained by each PC
pcvar <- (eig$values/sum(eig$values))*100

# plot PC1 vs PC2

pdf(file=paste0(outprefix,".pdf"))

plot(x=eig$vectors[,1], y=eig$vectors[,2], xlab=paste0("PC1 (",sprintf("%.2f",pcvar[1]), "%)"), ylab=paste0("PC2 (",sprintf("%.2f",pcvar[2]),"%)"), 
main="Lake Masoko A. calliptera", col=morphcol, lwd=1.5, cex=1.3, cex.lab=1.3, cex.axis=1.3) # PC1 vs PC2

legend('bottomright', c("Littoral","Benthic"), fill=c("yellow4","midnightblue"), bty='n',title="Ecomorph", cex=1.2)

invisible(dev.off())
