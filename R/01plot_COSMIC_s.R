rm(list=ls());

library(OmicCircos)

## colors
cols  <- rainbow(10, alpha=0.8) 


cnv.f  <- "../data_sets/COSMIC_CNV/COSO29914830_lung.txt"
gexp.f <- "../data_sets/COSMIC_EXP/COSO29914830_lung.txt"
fus.f  <- "../data_sets/COSMIC_FUSION/COSO29914830_lung.txt"
cnv    <- read.table(cnv.f, header=T)
gexp   <- read.table(gexp.f, header=T)
fus    <- read.table(fus.f, header=T)
## pdff <- paste0("../out/OC_COSO29914830_lung.pdf")
## pdf(pdff, 8, 8);
## output file name
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="");
## plot chromosome
circos(R=400, cir="hg19", W=4, type="chr", print.chr.lab=T, scale=T);
## plot heatmap
circos(R=260, cir="hg19", W=140, mapping=gexp[,c(1:20)],  col.v=4,  type="heatmap2", 
       lwd=0.1);
## plot CNV
circos(R=110, cir="hg19", W=160,  mapping=cnv, col.v=4, type="b2", 
       B=F, lwd=1, cutoff=2, col=cols[c(2,7)]);
## plot fussion
circos(R=120, cir="hg19", mapping=fus, type="link", lwd=1, col=cols[1]);
#dev.off();
