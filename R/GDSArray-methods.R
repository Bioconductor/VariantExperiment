#' @description \code{gdsfile}: on-disk location of GDS file
#'     represented by this object.
#' @param x GDSArray, GDSMatrix, GDSArraySeed or SummarizedExperiment object.
#' @rdname GDSArray
setGeneric("gdsfile", function(x) standardGeneric("gdsfile"))
#' @rdname GDSArray
#' @exportMethod gdsfile
setMethod("gdsfile", "GDSArraySeed", function(x) x@file)

#' @description \code{dim}, \code{dimnames}: dimension and dimnames of
#'     object contained in the GDS file.
#' @rdname GDSArray
#' @exportMethod dim
setMethod("dim", "GDSArraySeed", function(x) x@dim)

#' @rdname GDSArray
#' @exportMethod dimnames
setMethod("dimnames", "GDSArraySeed", function(x) x@dimnames)

#' @rdname GDSArray
#' @importFrom DelayedArray seed
#' @exportMethod gdsfile
setMethod("gdsfile", "GDSArray", function(x) gdsfile(seed(x)))

#' @rdname GDSArray
#' @exportMethod gdsfile
setMethod("gdsfile", "DelayedArray", function(x) gdsfile(seed(x)))

#' @rdname GDSArray
#' @exportMethod gdsfile
#' @import SeqArray
setMethod("gdsfile", "SummarizedExperiment", function(x) {
    vapply(assays(x), gdsfile, character(1))
})
