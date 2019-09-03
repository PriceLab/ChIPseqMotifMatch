library(ChIPseqMotifMatch)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("csmm")){
   bamFilename <- system.file(package="ChIPseqMotifMatch", "extdata", "GSM749704_hg19_Ctcf_chr19.bam")
   csmm <- ChIPseqMotifMatch(bamFilename, quiet=FALSE)
   }
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("ChIPseqMotifMatch") %in% is(csmm)))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
