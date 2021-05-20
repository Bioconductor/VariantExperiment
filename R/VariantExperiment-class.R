#' VariantExperiment class and slot getters and setters.
#' @title VariantExperiment-class
#' @rdname VariantExperiment-class
#' @description VariantExperiment could represent big genomic data in
#'     RangedSummarizedExperiment object, with on-disk GDS back-end
#'     data. The assays are represented by \code{DelayedArray}
#'     objects; \code{rowData} and \code{colData} could be represented
#'     by \code{DelayedDataFrame} or \code{DataFrame} objects.
#' @import SummarizedExperiment
#' @export

setClass(
    "VariantExperiment",
    contains="RangedSummarizedExperiment",
    )

###-------------- constructor --------------
#' @rdname VariantExperiment-class
#' @param assays A ‘list’ or ‘SimpleList’ of matrix-like elements, or
#'     a matrix-like object. All elements of the list must have the
#'     same dimensions, and dimension names (if present) must be
#'     consistent across elements and with the row names of
#'     ‘rowRanges’ and ‘colData’.
#' @param rowRanges A GRanges or GRangesList object describing the
#'     ranges of interest. Names, if present, become the row names of
#'     the SummarizedExperiment object. The length of the GRanges or
#'     GRangesList must equal the number of rows of the matrices in
#'     ‘assays’.
#' @param colData An optional DataFrame describing the samples. Row
#'     names, if present, become the column names of the
#'     VariantExperiment.
#' @param metadata An optional ‘list’ of arbitrary content describing
#'     the overall experiment.
#' @return a \code{VariantExperiment} object.
#' @details check "?RangedSummarizedExperiment" for more details.
#' @importFrom GenomicRanges GRangesList
#' @export VariantExperiment
#' 

VariantExperiment <- function(assays, rowRanges=GRangesList(),
                              colData=DelayedDataFrame(),
                              metadata=list())
{
    if (missing(assays))
        assays <- SimpleList()
    result <- SummarizedExperiment(
        assays=assays,
        rowRanges=rowRanges,
        colData = colData,
        metadata=as.list(metadata)
    )
    new("VariantExperiment", result)
}

## VariantExperiment <- function(assays, rowRanges=GRangesList(), colData=DelayedDataFrame(), metadata=list(), ...)
## {
##     ## se <- SummarizedExperiment:::.new_RangedSummarizedExperiment(
##     se <- VariantExperiment(
##         rowRanges=rowRanges,
##         colData=colData,
##         assays=assays,
##         ...)
##     ## rowData(se) <- as(rowData(se), "DelayedDataFrame")
##     ## colData(se) <- as(colData(se), "DelayedDataFrame")
##     ## mcols(rowRanges(se)) <- rowData(se)
##     as(se, "VariantExperiment")
## }

## setAs("SummarizedExperiment", "VariantExperiment", function(from)
## {
##     ## gf <- gdsfile(from)
##     ## .VariantExperiment(from, gdsfile = gf)
##     .VariantExperiment(from)
## })

###----------------
### Coercion
###----------------

#' @importFrom methods is
.validate_VariantExperiment <- function(x)
{
    ## ## DelayedDataFrame
    ## if(!all(is(rowData(x), "DelayedDataFrame"), is(colData(x), "DelayedDataFrame")))
    ##     return(wmsg("'rowData(x)' and 'colData(x)' must be DelayedDataFrame object"))

    ## GDSArray/DelayedArray for assay data.
    if(!all(vapply(assays(x), is, logical(1), "DelayedArray")))
        return(wmsg("'assays(x)' must be DelayedArray object"))
    
    ## ## gdsfile correlated with assay data
    ## if(is.character(gdsfile(x)))
    ##     return(wmsg("There should be an gds file correlated with 'x'"))
    TRUE
}


## Here only check @assay slot for 'DelayedArray' (can be any extensions)
#' @import S4Vectors 
setValidity2("VariantExperiment", .validate_VariantExperiment)

###--------------------
### getter and setter
###--------------------

## Still keep the `gdsfile()` function here, so that
## save/loadVariantExperiment() function works.

## `gdsfile()` function assumes that the VE comes from a single gds
## file or vcf file that was internally represented by a single gds
## file.

#' @export gdsfile
#' @rdname VariantExperiment-class
#' @param object a \code{VariantExperiment} object.

setMethod("gdsfile", "VariantExperiment", function(object)
    ## vapply(assays(object), gdsfile, character(1)))
    gdsfile(assays(object)[[1]])   ## here we assume all assay data are
                          ## correlated with the same gds file.
)

## ?? disable the "gdsfile" setter for now. Use
## "saveVariantExperiment" to save to a new file path.

#' @export "gdsfile<-"
#' @param value the new gds file path for VariantExperiment object.
#' @rdname VariantExperiment-class
setReplaceMethod("gdsfile", "VariantExperiment", function(object, value) {
    new_filepath <- tools::file_path_as_absolute(value)
    assays(object) <- lapply(assays(object), function(assay)
        BiocGenerics:::replaceSlots(seed(assay), file=value, check=FALSE))
    if (is(colData(object), "DelayedDataFrame")) {
        colData(object) <- DelayedDataFrame(lapply(colData(object), function(cols)
            BiocGenerics:::replaceSlots(seed(cols), file=value, check=FALSE)))
    }
    if (is(rowData(object), "DelayedDataFrame")) {
        rowData(object) <- DelayedDataFrame(lapply(rowData(object), function(cols)
            BiocGenerics:::replaceSlots(seed(cols), file=value, check=FALSE)))
    }
    object
})
