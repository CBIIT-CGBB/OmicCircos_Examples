## Remove all objects from the current R session to start with a clean environment
rm(list=ls())
## Load the OmicCircos library, which is used for generating circular plots for omics data
library(OmicCircos);
## load the data sets from OmicCircos package
## fusion data
data("TCGA.BC.fus");
## CNV data
data("TCGA.BC.cnv.2k.60");
## gene expression data
data("TCGA.BC.gene.exp.2k.60");
## sample ID and cancer types
data("TCGA.BC.sample60");
## p values for the association between CNV and gene expression
data("TCGA.BC_Her2_cnv_exp");
## transform p values by -log10
pvalue <- -1 * log10(TCGA.BC_Her2_cnv_exp[,5]);
## p value mapping data
pvalue <- cbind(TCGA.BC_Her2_cnv_exp[,c(1:3)], pvalue);
## Her2 sub-type data
Her2.i <- which(TCGA.BC.sample60[,2] == "Her2");
Her2.n <- TCGA.BC.sample60[Her2.i,1];
## CNV mapping data for Her2 sub-type
Her2.j <- which(colnames(TCGA.BC.cnv.2k.60) %in% Her2.n);
cnv    <- TCGA.BC.cnv.2k.60[,c(1:3,Her2.j)]; 
cnv.m  <- cnv[,c(4:ncol(cnv))];
cnv.m[cnv.m >  2] <- 2;
cnv.m[cnv.m < -2] <- -2;
cnv <- cbind(cnv[,1:3], cnv.m);
## gene expression mapping data for Her2 sub-type
Her2.j   <- which(colnames(TCGA.BC.gene.exp.2k.60) %in% Her2.n);
gene.exp <- TCGA.BC.gene.exp.2k.60[,c(1:3,Her2.j)]; 
## initial colors
colors <- rainbow(10, alpha=0.5);
## Set the margins of the plot
par(mar=c(0, 0, 0, 0));
## Initialize a blank plot with custom dimensions
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
## Plot the chromosomes on the circular plot with specified radius, chromosome annotation, and scaling
circos(R=400, cir="hg18", W=4,   type="chr", print.chr.lab=TRUE, scale=TRUE);
## Plot heatmap
circos(R=300, cir="hg18", W=100, mapping=gene.exp,  col.v=8,  type="heatmap2", 
       cluster=FALSE, col.bar=TRUE, lwd=0.1, col="blue");
## Plot CNV
circos(R=220, cir="hg18", W=80,  mapping=cnv,   col.v=4,   type="ml3", B=FALSE, lwd=1, cutoff=0);
## Plot CNV
circos(R=140, cir="hg18", W=80,  mapping=pvalue,  col.v=4,    type="l",   B=TRUE, lwd=1, col=colors[1]);
## Add fusion data to the circular plot as links
cols        <- rep(colors[7], nrow(TCGA.BC.fus));
col.i       <- which(TCGA.BC.fus[,1]==TCGA.BC.fus[,4]);
cols[col.i] <- colors[1];
circos(R=130, cir="hg18", W=10,  mapping=TCGA.BC.fus, type="link2", lwd=2, col=cols);

