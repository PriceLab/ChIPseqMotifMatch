#------------------------------------------------------------------------------------------------------------------------
#' @import methods
#'
#' @title ChIPseqMotifMatch-class
#'
#' @name ChIPseqMotifMatch-class
#' @rdname ChIPseqMotifMatch-class
#' @aliases ChIPseqMotifMatch
#' @exportClass ChIPseqMotifMatch
#'

.ChIPseqMotifMatch <- setClass("ChIPseqMotifMatch",
                               representation=representation(
                                  bamFilename="character",
                                  state="environment",
                                  quiet="logical"
                                  ))

#------------------------------------------------------------------------------------------------------------------------
setGeneric('getBamFilename', signature='obj', function(obj) standardGeneric('getBamFilename'))
setGeneric('setRegionsOfInterest', signature='obj', function(obj, tbl.regions) standardGeneric('setRegionsOfInterest'))
setGeneric('getRegionsOfInterest', signature='obj', function(obj) standardGeneric('getRegionsOfInterest'))
setGeneric('calculateNarrowPeaks', signature='obj', function(obj) standardGeneric('calculateNarrowPeaks'))
#------------------------------------------------------------------------------------------------------------------------
#' Define an object of class ChIPseqMotifMatch
#'
#' @description
#' Survey some or all of a bam file, scoring pwm motif match to sequence
#'   where ChIP hits have been called by MACS2
#'
#' @rdname ChIPseqMotifMatch-class
#'
#' @param bamFilename character
#' @param quiet logical
#'
#' @export
#'
#' @return An object of the ChIPseqMotifMatch class
#'
ChIPseqMotifMatch <- function(bamFilename, quiet=TRUE)

{
   stopifnot(file.exists(bamFilename))
   state <- new.env(parent=emptyenv())
   state$roi <- data.frame()

   .ChIPseqMotifMatch(bamFilename=bamFilename, state=state, quiet=quiet)

} # ChIPseqMotifMatch, the constructor
#------------------------------------------------------------------------------------------------------------------------
#' Get all the name of the bamfile
#'
#' @rdname getBamFilename
#' @aliases getBamFilename
#'
#' @param obj An object of class ChIPseqMotifMatch
#'
#' @export

setMethod('getBamFilename',  'ChIPseqMotifMatch',

     function(obj){
       obj@bamFilename
       })

#------------------------------------------------------------------------------------------------------------------------
#' specify the chromosomal regions you wish to examine and operat upon
#'
#' @rdname setRegionsOfInterest
#' @aliases setRegionsOfInterest
#'
#' @param obj An object of class ChIPseqMotifMatch
#' @param tbl.regions a data.frame with chrom, start, end columns
#'
#' @export

setMethod('setRegionsOfInterest',  'ChIPseqMotifMatch',

     function(obj, tbl.regions){
        stopifnot(is.data.frame(tbl.regions))
        stopifnot(nrow(tbl.regions) >= 1)
        stopifnot(colnames(tbl.regions) == c("chrom", "start", "end"))
        stopifnot(is.character(tbl.regions$chrom))
        obj@state$roi <- tbl.regions
        })

#------------------------------------------------------------------------------------------------------------------------
#' get the chromosomal regions of interest
#'
#' @rdname getRegionsOfInterest
#' @aliases getRegionsOfInterest
#'
#' @param obj An object of class ChIPseqMotifMatch
#'
#' @return a (possibly empty) data.frame
#
#' @export

setMethod('getRegionsOfInterest',  'ChIPseqMotifMatch',

     function(obj){
        obj@state$roi
        })

#------------------------------------------------------------------------------------------------------------------------
#' run MACS2 on the current regions of interest in the bam file using the narrowPeaks options
#'
#' @rdname calculateNarrowPeaks
#' @aliases calculateNarrowPeaks
#'
#' @param obj An object of class ChIPseqMotifMatch
#'
#' @return a (possibly empty) data.frame
#
#' @export

setMethod('calculateNarrowPeaks',  'ChIPseqMotifMatch',

     function(obj){
        stopifnot(nrow(obj@state$roi) >= 1)
           # mariam: use all your magic and good existing code to implement this method!
           # i offer a simple placeholder for now, returning an empty data.frame
        return(data.frame())
        })

#------------------------------------------------------------------------------------------------------------------------
