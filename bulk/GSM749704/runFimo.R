library(MotifDb)
source(file.path(Sys.getenv("HOME"), "github", "fimoService", "batchMode", "fimoBatchTools.R"))
library(igvR)
#------------------------------------------------------------------------------------------------------------------------
enableMotifDisplay <- function(igv)
{
  #body <- sprintf("return('%s')", "<h3>fobo</h3>")
  url <- "http://jaspar.genereg.net/static/logos/svg/MA0803.1.svg"
  imageTag <- sprintf("<img src=\"%s\" width=\"200\" />", url)
  imageTag <- sprintf("<img src='%s' width='200' />", url)
  body <- sprintf('{return("%s")}', imageTag)

  body.parts <- c(
     'var returnValue = undefined;',
     'popoverData.forEach(function(i){',
     '   if(i.name=="name" && i.value.startsWith("http:")){',
     '      var url = i.value;',
     '      console.log(url);',
     '      var tag = "<img src=\'" + url + "\' width=300\'/>";',
     '      console.log(tag);',
     '      returnValue=tag;',
     '      };',
     '   });',
     '   console.log("--- returnValue:");',
     '   console.log(returnValue);',
     '   return(returnValue);'
     )

  body <- paste(body.parts, collapse=" ")
  x <- list(arguments="track, popoverData", body=body)
  setTrackClickFunction(igv, x)

} # enableMotifDisplay
#------------------------------------------------------------------------------------------------------------------------
if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg19")
   pfm <- query(MotifDb, c("sapiens", "CTCF", "MA0139", "jaspar2018"))
   export(pfm, "ctcf-human.meme", 'meme')
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.narrowPeaks")){
   tbl.narrowPeaks <- read.table("GSM749704_peaks.narrowPeak", sep="\t", as.is=TRUE)
   colnames(tbl.narrowPeaks) <- c("chrom", "start", "end", "name", "score", "strand", "foldChange",
                                  "pScore", "qScore", "summitPos")
   }
#------------------------------------------------------------------------------------------------------------------------
displayRegion <- function(chromosome, start.loc, end.loc, fimo.threshold)
{
   showGenomicRegion(igv, sprintf("%s:%d-%d", chromosome, start.loc, end.loc))
   tbl.regions <- subset(tbl.narrowPeaks, chrom==chromosome & start >= start.loc & end <= end.loc)
   dim(tbl.regions)

   tbl.match <- fimoBatch(tbl.regions, matchThreshold=fimo.threshold, genomeName="hg19", pwmFile="ctcf-human.meme")
   dim(tbl.match)
   if(nrow(tbl.match) > 0){
      tbl.matchScored <- tbl.match[, c("chrom", "start", "end", "p.value")]
      tbl.matchScored$p.value <- -log10(tbl.matchScored$p.value)
      track <- DataFrameQuantitativeTrack("motif.score", tbl.matchScored, autoscale=TRUE, color="brown")
      displayTrack(igv, track)

      tbl.match$tf <- sprintf("MotifDb::%s", names(pfm))
      tbl.match <- tbl.match[, c("chrom", "start", "end", "tf")]
      colnames(tbl.match)[4] <- "name"
      dim(tbl.match)
      track <- DataFrameAnnotationTrack("motif", tbl.match, color="brown")
      displayTrack(igv, track)
      }

   start.loc <- min(tbl.regions$start - 1000)
   end.loc <- max(tbl.regions$end + 1000)
   chromosome <- tbl.regions$chrom[1]
   #showGenomicRegion(igv, sprintf("%s:%d-%d", chromosome, start.loc, end.loc))

   track <- DataFrameQuantitativeTrack("np", tbl.regions[, c(1,2,3,5)], autoscale=TRUE, color="darkGreen")
   displayTrack(igv, track)
   tbl.summits.roi <- subset(tbl.summits, chrom==chromosome & start>=start.loc & end<=end.loc)[, c(1,2,3,5)]
   dim(tbl.summits.roi)
   track <- DataFrameQuantitativeTrack("summit", tbl.summits.roi, autoscale=TRUE, color="blue")
   displayTrack(igv, track)


   bamFile <- "GSM749704_hg19_wgEncodeUwTfbsGm12878CtcfStdAlnRep1.bam"
   stopifnot(file.exists(bamFile))

   which <- GRanges(seqnames=chromosome, IRanges(start.loc, end.loc))
   param <- ScanBamParam(which=which, what = scanBamWhat())
   x <- readGAlignments(bamFile, use.names=TRUE, param=param)
   track <- GenomicAlignmentTrack("ChIP", x, visibilityWindow=1000000)
   displayTrack(igv, track)

   tbl.summits <- read.table("GSM749704_summits.bed", sep="\t", as.is=TRUE)
   colnames(tbl.summits) <- c("chrom", "start", "end", "name", "score")
   dim(tbl.summits)
   head(tbl.summits)

   #loc <- getGenomicRegion(igv)
   #with(loc, {chromosome <- chrom; start.loc<-start; end.loc<-end})

} # displayRegion
#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
   chromosome <- "chr14"
   start.loc <- 21560668
   end.loc <- 21561002
   #displayRegion("chr14", 21560668, 21561002)
   loc <- getGenomicRegion(igv)
   showGenomicRegion(igv, sprintf("%s:%d-%d", loc$chrom, loc$start, loc$end))
   with(loc, displayRegion(chrom, start, end, fimo.threshold=1e-2))

} # run
#------------------------------------------------------------------------------------------------------------------------
# interesting locations:
#  1) chr14:21,913,061-21,951,891
#       peak     np       fimo (-log10 pval
#          1    297         5.4
#          2    768         2.1
