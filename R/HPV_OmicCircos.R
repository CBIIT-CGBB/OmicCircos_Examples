rm(list=ls());

set.seed(1234);
options(stringsAsFactors = FALSE);
library(OmicCircos);

## human data
data(UCSC.hg19.chr);
UCSC.hg19.chr[,1] <- gsub("chr", "Chr", UCSC.hg19.chr[,1])

## HPV data
HPV     <- read.table("HPV_anno.txt", header=T);
HPV[,5] <- "NA";

## HPV gene labels
HPV.gene   <- data.frame(name=HPV[,1], po=(HPV[,2]+HPV[,3])/2, gene=HPV[,4]);
HPV.gene.l <- HPV;
HPV.gene.l[,5] <- sample(c(1:nrow(HPV)), nrow(HPV));

## generating db with genomes of human and HPV 
hs.db  <- segAnglePo(UCSC.hg19.chr, seg=unique(UCSC.hg19.chr[,1]), angle.start=1, angle.end=270);
hpv.db <- segAnglePo(HPV, seg="HPV", angle.start=271, angle.end=360);
hpv.db[1,4] <- as.numeric(hpv.db[1,4]) + as.numeric(hs.db[24,5])+1;
hpv.db[1,5] <- as.numeric(hpv.db[1,4]) + as.numeric(hpv.db[1,5]);
db          <- rbind(hs.db, hpv.db);

## human chromosome label
ref     <- db;
ref[,1] <- gsub("chr", "", ref[,1]);
chr     <- unique(ref[,1]);
chr.l   <- c();
for (ch in chr){
  dat.i <- which(ref[,1]==ch);
  m     <- max(as.numeric(ref[dat.i,7]))/2;
  chr.l <- rbind(chr.l, c(ch, m, ch));
}

## input data
dat  <- read.table("HPV_simu.txt", header=T);
## human insert positions
dat2 <- dat[,c(4:6)];

## colors
col1    <- c(rainbow(24), rgb(0.5, 0.5, 0.5, 0.5));
col2    <- rainbow(nrow(HPV));
col3    <- rainbow(10, alpha=0.9);

## link colors
hpv.c   <- data.frame(HPV, col=rainbow(nrow(HPV), alpha=0.3));
out.s   <- c();
out  <- sapply(hpv.c$end, '>=', dat$start1) & sapply(hpv.c$start, '<=', dat$end1);
for (i in 1:ncol(out)){
  n.j <- which(out[,i]);
  if (length(n.j)==0){
    next;
  }
  for (k in 1:length(n.j)){
    j <- n.j[k];
    out.s  <- rbind(out.s, c(hpv.c[i,], dat[j,])); 
  }
}
dat.s <- out.s[,c(7:12)];
col4  <- as.character(out.s[,6]);

## plot using OmicCircos
pdf("HPV_OmicCircos20160809_v2.pdf", 8, 8)
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

## chromosomes
circos(R=310, type="chr", cir=db, W=20, scale=F, print.chr.lab=F, col=col1);
## human chr label
circos(R=320, cir=db, type="label2", W=30, mapping=chr.l, side="out", cex=1, col=c("black"));
## HPV gene region
circos(R=260, cir=db, W=50, mapping=HPV.gene.l[-c(1,10),], col.v=5, type="arc", col=col2[-c(1,10)], lwd=8);
## human insert position
circos(R=260, cir=db, W=50, mapping=dat2, type="b3", col=col3[7], lwd=1);
## links of HPV segments and human insert positions
circos(R=250, cir=db, W=40, mapping=dat.s, type="link.pg", lwd=2, col=col4);

legend("topleft", legend=HPV[-c(1,10),4], col=col2[-c(1,10)], lwd=6, bty="n", title="HPV Genes");
dev.off();
