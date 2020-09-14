# Please make sure to be using R-3.4 or later

args <- commandArgs(TRUE)
source("https://bioconductor.org/biocLite.R")
if (!require("Tnseq")) biocLite("Tnseq")
library(Tnseq)

dat <- read.table(args[1], header=T, sep='\t') 
loci <- dat$loci
counts <- dat$counts
counts.adjusted <- list()
counts.adjusted[[1]] <- loci
biasfactor <- BiasFactor(loci, counts, window=as.numeric(args[2]))
counts.adjusted[[2]] <- counts/biasfactor$bias.factor
counts.adjusted.df <- data.frame(counts.adjusted)
colnames(counts.adjusted.df) <- c("loci", "counts.adjusted")
write.table(counts.adjusted.df, file=args[3], sep='\t', row.names=F)