library(MotifDb)
source("fimoBatchTools.R")
library(igvR)
library(phastCons100way.UCSC.hg38); phast.100 <- phastCons100way.UCSC.hg38
library(phastCons7way.UCSC.hg38); phast.7 <- phastCons7way.UCSC.hg38

library(phastCons100way.UCSC.hg19); phast.100.hg19 <- phastCons100way.UCSC.hg19 #just in case you need hg19

#------------------------------------------------------------------------------------------------------------------------
# a convenience function
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.narrowPeaks")){
  tbl.narrowPeaks <- read.table("ctcf__peaks.narrowPeak", sep="\t", as.is=TRUE)
  colnames(tbl.narrowPeaks) <- c("chrom", "start", "end", "name", "score", "strand", "foldChange",
                                 "pScore", "qScore", "summitPos")
}

if(!file.exists("ctcf-human.meme")){
   motif <- query(MotifDb, c("ctcf", "sapiens", "jaspar2018", "MA0139"))
   export(motif, con="ctcf-human.meme", format="meme")
   }
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_addFimoToNarrowPeaks()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
# 1e-4 is the traditional FIMO pValue threshold: only motif-to-sequence matches better than this
# are returned.  we find that loosening this to at least 1e-3 is needed to find motifs in high-scoring
# MACS2-called narrow peaks.
addFimoToNarrowPeaks <- function(tbl.narrowPeaks, chromosome, fimo.threshold, numberOfPeaksToConsider=-1, plot=FALSE)
{
   tbl.np <- subset(tbl.narrowPeaks, chrom==chromosome)

   if(numberOfPeaksToConsider > 0)   # this makes for quick testing
      tbl.np <- tbl.np[seq_len(numberOfPeaksToConsider),]

   printf("narrowPeak rows: %d", nrow(tbl.np))
   print(system.time(tbl.match <- fimoBatch(tbl.np, matchThreshold=fimo.threshold, genomeName="hg19",
                                            pwmFile="ctcf-human.meme")))
   printf(" matches from fimo: %d", nrow(tbl.match))

   tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.match), GRanges(tbl.np), type="within"))
   printf("overlaps: %d", nrow(tbl.ov))
   colnames(tbl.ov) <- c("fimo", "narrowPeak")

   tbl.noMatch <- data.frame(motif.chrom=NA, motif.start=NA_integer_, motif.end=NA_integer_,
                             motif.tf=NA, motif.strand=NA, motif.score=NA_real_, motif.p.value=1,
                             motif.matched_sequence=NA, motif.motif_id=NA)


     # a motif often matches at multiple points within a sequence.  we want just the best match
     # this function is used in the lapply which follows

   bestFimoRow <- function(narrowPeakRow){
      fimoRows=subset(tbl.ov, narrowPeak==narrowPeakRow)$fimo;
      tbl.matchThisPeak <- tbl.match[fimoRows,];
      tbl.bestMatch <- subset(tbl.matchThisPeak, score==max(tbl.matchThisPeak$score))[1,]
      grep(tbl.bestMatch$start, tbl.match$start)[1]
      }

   x <- lapply(unique(tbl.ov$narrowPeak), bestFimoRow)
   tbl.ov.trimmed <- data.frame(fimo=unlist(x), narrowPeak=unique(tbl.ov$narrowPeak))
   unmatched.peaks <- setdiff(seq_len(nrow(tbl.np)), tbl.ov.trimmed$narrowPeak)
   length(unmatched.peaks)

   tbl.matchFiltered <- tbl.match[tbl.ov.trimmed$fimo,]
   dim(tbl.matchFiltered)
   colnames(tbl.matchFiltered) <- paste0("motif.", colnames(tbl.matchFiltered))

   tbl.both <- cbind(tbl.np[tbl.ov.trimmed$narrowPeak,], tbl.matchFiltered, row.names=NULL)
   dim(tbl.both)

   if(length(unmatched.peaks) > 0){
      printf("unmatched.peaks: %d", length(unmatched.peaks))
      tbl.unmatched <- cbind(tbl.np[unmatched.peaks,], tbl.noMatch, row.names=NULL)
      tbl.both <- rbind(tbl.both, tbl.unmatched)
      }

   dim(tbl.both)
   new.order <- order(tbl.both$start, decreasing=FALSE)
   tbl.both <- tbl.both[new.order,]

   checkEquals(tbl.both$start, tbl.np$start)
   if(plot)
      with(tbl.both, plot(score, -log10(motif.p.value), main=chromosome))

   tbl.both

} # addFimoToNarrowPeaks
#------------------------------------------------------------------------------------------------------------------------
test_addFimoToNarrowPeaks <- function()
{
   printf("--- test_addFimoToNarrowPeaks")

   tbl.10.3 <- addFimoToNarrowPeaks(tbl.narrowPeaks, "chr19", 1e-2, 10)

   checkEquals(dim(tbl.10.3), c(10,19))


} # test_addFimoToNarrowPeaks
#------------------------------------------------------------------------------------------------------------------------
addConserverationScores <- function(chrom, start, end, tbl)
{
   tbl.try2 <- subset(tbl, !is.na(motif.start))
   gr.try2 <- with(tbl.try2, GRanges(seqnames=chrom, IRanges(start=motif.start, end=motif.end)))
   gscores(phast.7, gr.try2)
   tbl.phast7 <- as.data.frame(gscores(phast.7, gr.try2))
   tbl.merged <- merge(tbl, tbl.phast7, by.x="motif.start", by.y="start", all.x=TRUE)
   grep("default", colnames(tbl.merged))
   colnames(tbl.merged)[24] <- "phast7"
   tbl.merged <- tbl.merged[order(tbl.merged$start),]
   
 #add the phast.100 column
  gscores(phast.100, gr.try2)
  tbl.phast100 <- as.data.frame(gscores(phast.100, gr.try2))
  tbl.merged <- merge(tbl.merged, tbl.phast100, by.x="motif.start", by.y="start", all.x=TRUE)
  grep("default", colnames(tbl.merged))
  colnames(tbl.merged)[29] <- "phast100" #renaming the column from "default" to "phast100" 
  tbl.merged <- tbl.merged[order(tbl.merged$start),]
  head(tbl.merged) 
   
} # addConservationScores
#------------------------------------------------------------------------------------------------------------------------
#chrom1
tbl.10.4 <- addFimoToNarrowPeaks(tbl.narrowPeaks, "chr1", 1e-2, 10)
tbl.chrom1<-addConserverationScores(chrom="chr1", start=1, end=NA,tbl.10.4 )

#chrom2
tbl.10.5 <- addFimoToNarrowPeaks(tbl.narrowPeaks, "chr2", 1e-2, 10)
tbl.chrom2<-addConserverationScores(chrom="chr2", start=1, end=NA,tbl.10.5 )

tbl.combined<- rbind(tbl.chrom1, tbl.chrom2) #combines the two tables into one table

chrom="chr2"
for (i in 3:22) { #how to combine the rest of the chromosomes in the table
  tbl.10.i <- addFimoToNarrowPeaks(tbl.narrowPeaks, gsub("2",i,chrom), 1e-2, 10)
  tbl.chrom.i<-addConserverationScores(chrom=gsub("2",i,chrom), start=1, end=NA,tbl.10.i)
  tbl.combined<- rbind(tbl.combined, tbl.chrom.i)
}
tbl.combined

new.tbl<-na.exclude(tbl.combined)
cor(new.tbl$phast7, new.tbl$phast100)
