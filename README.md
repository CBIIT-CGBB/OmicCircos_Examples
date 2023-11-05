# OmicCircos_Examples

<img src="out/OMIC_SampleID1.txt.png" width="600" height="600"> 

[R codes](R/03OMIC_CNV1_plot.R)


```r
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
```

<img src="out/04OMIC_CNV2_plot.png" width="600" height="600"> 

[R codes](R/03OMIC_CNV2_plot.R)

```r
rm(list=ls());

library(OmicCircos)

## select samples
cfile <- dir("../data_sets/OMIC_CNV1", "txt")
## colors
cols  <- rainbow(10, alpha=0.5)[c(1,7,2,9,6)]
## output file
pdff <- paste0("../out/04OMIC_CNV2_plot.pdf")
pdf(pdff, 8, 8);
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="")
## plot chromosome
circos(R=400, cir="hg19", W=10, type="chr", print.chr.lab=T, scale=T);

R <- 340
for (i in 1:length(cfile)){
  id.n   <- cfile[i]
  ## input CNV file
  infile <- paste0("../data_sets/OMIC_CNV2/", cfile[i])
  ## read CNV data
  cnv    <- read.table(infile, sep="\t", header=T)
  ## calculate cutoff of CNV
  c.m    <- mean(cnv[,5])
  ## input fusion file
  inf2   <- paste0("../data_sets/OMIC_FUSION/", cfile[i]);
  ## read fusion data
  dat    <- read.table(inf2, sep="\t", header=T)
  b.i <- i%%2
  ## plot CNV
  if (b.i==1){
    circos(R=R, cir="hg19", W=40,  mapping=cnv, col.v=5, type="ml3", 
           B=F, lwd=2, cutoff=2, col=cols[i]);
  } else {
    circos(R=R, cir="hg19", W=40,  mapping=cnv, col.v=5, type="ml3", 
           B=T, lwd=2, cutoff=2, col=cols[i]);
  }
  ## plot fusions
  R <- R - 40
  circos(R=160, cir="hg19", mapping=dat, type="link", lwd=2, col=cols[i]);
}

legend("topleft", legend=paste0("SampleID1-5"), 
       bty="n", cex=0.8, title="Outside(CNV)")
legend("topright", legend=paste0("SampleID",1:5), lwd=2, col=cols, 
       bty="n", cex=0.8, title="Inside(Fusion)")

dev.off()

```

<img src="out/OC_COSO29914826_lung.png" width="600" height="600"> 

[R codes](R/01plot_COSMIC.R)

```r
rm(list=ls());

library(OmicCircos)
## data set names
COSMIC_PHENOTYPE_ID  <- c("COSO29914830", "COSO32054826", "COSO29914826")
PRIMARY_SITE <- c("lung", "prostate", "lung")
## colors
cols  <- rainbow(10, alpha=0.8) 

for (p.i in 1:length(COSMIC_PHENOTYPE_ID)){
  tid    <- COSMIC_PHENOTYPE_ID[p.i]
  tname  <- PRIMARY_SITE[p.i]
  ## read input data sets
  cnv.f  <- paste0("../data_sets/COSMIC_CNV/", tid, "_", tname, ".txt")
  gexp.f <- paste0("../data_sets/COSMIC_EXP/", tid, "_", tname, ".txt")
  fus.f  <- paste0("../data_sets/COSMIC_FUSION/", tid, "_", tname, ".txt")
  cnv    <- read.table(cnv.f, header=T)
  gexp   <- read.table(gexp.f, header=T)
  fus    <- read.table(fus.f, header=T)
  ## output file name
  pdff <- paste0("../out/OC_", tid, "_", tname, ".pdf")
  pdf(pdff, 8, 8);
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
  dev.off()
}

```


