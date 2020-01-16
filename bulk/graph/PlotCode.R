source("AddFIMOtoMACS2NarrowPeaks.R")
#chrom1
#if(!exists("tbl.narrowPeaks")){
tbl.narrowPeaks <- read.table("ctcf__peaks.narrowPeak", sep="\t", as.is=TRUE)

colnames(tbl.narrowPeaks) <- c("chrom", "start", "end", "name", "score", "strand", "foldChange",
                               "pScore", "qScore", "summitPos")
#}

if(!file.exists("ctcf-human.meme")){
  motif <- query(MotifDb, c("ctcf", "sapiens", "jaspar2018", "MA0139"))
  export(motif, con="ctcf-human.meme", format="meme")
}

tbl.chrom1 <- addFimoToNarrowPeaks(tbl.narrowPeaks, "chr1", 1e-2, -1)

#chrom2
tbl.chrom2<- addFimoToNarrowPeaks(tbl.narrowPeaks, "chr2", 1e-2, -1) 

tbl.combined<- rbind(tbl.chrom1, tbl.chrom2) #combines the two tables into one table

chrom="chr2"
for (i in 3:22) { #how to combine the rest of the chromosomes in the table
  tbl.chrom.i <- addFimoToNarrowPeaks(tbl.narrowPeaks, gsub("2",i,chrom), 1e-2, -1)
  tbl.combined<- rbind(tbl.combined, tbl.chrom.i)
} 
tbl.combined

yaxis.before=tbl.combined$motif.p.value
yaxis=-log10(yaxis.before)
plot(tbl.combined$score, yaxis, main="GSM749704 ChIPseq CTCF- best FIMO match", ylab="-log10FIMO", xlab="MACS2 ChIP narrow peak score", col= ifelse(yaxis < 4, "red", ifelse(yaxis >= 4, "black", "black")))
