rm(list=ls());

library(OmicCircos)
## data set names
COSMIC_PHENOTYPE_ID  <- c("COSO29914830", "COSO32054826", "COSO29914826")
PRIMARY_SITE <- c("lung", "prostate", "lung")
## colors
cols   <- rainbow(10, alpha=0.8)[c(1,2,9,7)] 
## read cancer gene data
gene   <- read.table("../data_sets/COSMIC_GENE/gene.txt", header=T)
row.i  <- 1:nrow(gene)
g.i    <- which(row.i%%2==1)
## output file name
pdff <- paste0("../out/02plot_COSMIC_gene.pdf")
pdf(pdff, 8, 8);
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="");
## plot chromosome
circos(R=290, cir="hg19", W=4, type="chr", print.chr.lab=T, scale=T);
## plot genes
circos(R=325, cir="hg19", W=20, mapping=gene[g.i,], type="label", 
       side="out", col=c("black"), cex=c(0.4));
circos(R=280, cir="hg19", W=20, mapping=gene[-g.i,], type="label", 
       side="in", col=c("black"), cex=c(0.4));
## plot fussions
for (p.i in 1:length(COSMIC_PHENOTYPE_ID)){
  tid    <- COSMIC_PHENOTYPE_ID[p.i]
  tname  <- PRIMARY_SITE[p.i]
  fus.f  <- paste0("../data_sets/COSMIC_FUSION/", tid, "_", tname, ".txt")
  fus  <- read.table(fus.f, header=T)
  circos(R=200, cir="hg19", mapping=fus, type="link", lwd=1, col=cols[p.i]);
}

dev.off();

