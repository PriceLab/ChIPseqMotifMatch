library(ChIPseqMotifMatch)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
test.bamFilename <- system.file(package="ChIPseqMotifMatch", "extdata", "GSM749704_hg19_Ctcf_chr19.bam")
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_constructorMissingBamFile()
   test_getSetRegionsOfInterest()
   test_calculateNarrowPeaks()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   csmm <- ChIPseqMotifMatch(test.bamFilename, quiet=FALSE)
   checkTrue(all(c("ChIPseqMotifMatch") %in% is(csmm)))
   checkEquals(getBamFilename(csmm), test.bamFilename)

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_constructorMissingBamFile <- function()
{
   message(sprintf("--- test_constructorMissingBamFile"))

   checkException(ChIPseqMotifMatch("bogus.bam"), silent=TRUE)

} # test_constructorMissingBamFile
#------------------------------------------------------------------------------------------------------------------------
test_getSetRegionsOfInterest <- function()
{
   message(sprintf("--- test_getSetRegionsOfInterest"))

   csmm <- ChIPseqMotifMatch(test.bamFilename, quiet=FALSE)
   checkEquals(getRegionsOfInterest(csmm), data.frame())

     # set should fail if data.frame is empty, or if the columns are wrong, if the chrom column is a factor
   checkException(setRegionsOfInterest(csmm, data.frame()), silent=TRUE)
   checkException(setRegionsOfInterest(csmm, data.frame(chrom="chr1", start=20, mistake=40, stringsAsFactors=FALSE)), silent=TRUE)
   checkException(setRegionsOfInterest(csmm, data.frame(chrom="chr1", start=20, end=40, stringsAsFactors=TRUE)), silent=TRUE)

   tbl.regions <- data.frame(chrom=rep("chr19", 6),
                             start=c(42257038, 43033314, 19256647, 42210396, 35632204, 42937212),
                             end=  c(42258277, 43034146, 19257529, 42211437, 35634113, 42945062),
                             stringsAsFactors=FALSE)
   setRegionsOfInterest(csmm, tbl.regions)
   checkEquals(getRegionsOfInterest(csmm), tbl.regions)

} # test_constructorMissingBamFile
#------------------------------------------------------------------------------------------------------------------------
test_calculateNarrowPeaks <- function()
{
   message(sprintf("--- test_calculateNarrowPeaks"))
   csmm <- ChIPseqMotifMatch(test.bamFilename, quiet=FALSE)
   tbl.regions <- data.frame(chrom="chr19", start=42257038, end=42258277, stringsAsFactors=FALSE)
   setRegionsOfInterest(csmm, tbl.regions)
   tbl.narrowPeaks <- calculateNarrowPeaks(csmm)

   checkTrue(is.data.frame(tbl.narrowPeaks))
      # add more tests as you implement the method

} # test_calculateNarrowPeaks
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
