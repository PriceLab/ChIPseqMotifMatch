#----------------------------------------------------------------------------------------------------
#' @import methods
#' @import org.Hs.eg.db
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
                                  quiet="logical"
                                  ))

#----------------------------------------------------------------------------------------------------
#' Define an object of class ChIPseqMotifMatch
#'
#' @description
#' Survey some or all of a bam file, scoring pwm motif match to sequence
#'   where ChIP hits have been called by MACS2
#'
#' @rdname ChIPseqMotifMatch-class
#'
#' @export
#'
#' @return An object of the ChIPseqMotifMatch class
#'
ChIPseqMotifMatch <- function(bamFilename, quiet=TRUE)

{
   .ChIPseqMotifMatch(bamFilename=bamFilename, quiet=quiet)

} # ChIPseqMotifMatch, the constructor
#----------------------------------------------------------------------------------------------------
