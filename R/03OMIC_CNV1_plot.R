rm(list=ls());

library(OmicCircos)
## input files
cfile <- dir("../data_sets/OMIC_CNV1", "txt")
## colors
cols  <- rainbow(10, alpha=0.5)

for (i in 1:length(cfile)){
  id.n   <- cfile[i]
  ## input CNV file
  infile <- paste0("../data_sets/OMIC_CNV1/", cfile[i])
  ## read CNV data
  cnv    <- read.table(infile, sep="\t", header=T)
  ## input fusion file
  inf2   <- paste0("../data_sets/OMIC_FUSION/", cfile[i]);
  ## read fusion data
  dat    <- read.table(inf2, sep="\t", header=T)
  ## calculate cutoff of CNV
  c.m    <- mean(cnv[,3])
  ## output file
  pdff <- paste0("../out/", id.n, ".pdf")
  pdf(pdff, 8, 8);
  par(mar=c(2, 2, 2, 2));
  plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="", main=id.n);
  ## plot chromosome
  circos(R=350, cir="hg19", W=10, type="chr", print.chr.lab=T, scale=T);
  ## plot CNV
  circos(R=200, cir="hg19", W=160,  mapping=cnv, col.v=3, type="b2", 
         B=F, lwd=1, cutoff=c.m, col=cols[c(2,7)]);
  ## plot fusion
  circos(R=200, cir="hg19", mapping=dat, type="link", lwd=2, col=cols[1]);
  dev.off()
}