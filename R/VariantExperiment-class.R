#' VariantExperiment
#' 
#' @name VariantExperiment-class
#' @description VariantExperiment could represent big genomic data in RangedSummarizedExperiment object, with on-disk GDS back-end data. The assays are represented by \code{DelayedArray} objects; \code{rowData} and \code{colData} could be represented by \code{DelayedDataFrame} objects.
#' @importClassesFrom SummarizedExperiment SummarizedExperiment RangedSummarizedExperiment 
#' @exportClass VariantExperiment
#' @rdname VariantExperiment-class
#' 
.VariantExperiment <- setClass(
    "VariantExperiment",
    contains="RangedSummarizedExperiment",
)

###--------------
### constructor
###--------------
#' @export VariantExperiment
#' @importFrom GenomicRanges GRangesList
#' @importFrom SummarizedExperiment Assays
#' @rdname VariantExperiment-class
VariantExperiment <- function(assays, rowRanges=GRangesList(), colData=DelayedDataFrame(), metadata=list())
{
    result <- SummarizedExperiment(
        rowRanges=rowRanges,
        colData = colData,
        assays=assays,
        metadata=as.list(metadata)
    )
    .VariantExperiment(result)
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
##     ## gf <- gdsfile(from)[[1]]
##     ## .VariantExperiment(from, gdsfile = gf)
##     .VariantExperiment(from)
## })

###----------------
### Coercion
###----------------

## 

.validate_VariantExperiment <- function(x)
{
    ## DelayedDataFrame
    if(!all(is(rowData(x), "DelayedDataFrame"), is(colData(x), "DelayedDataFrame")))
        return(wmsg("'rowData(x)' and 'colData(x)' must be DelayedDataFrame object"))

    ## GDSArray for assay data.
    if(!all(vapply(assays(x), is, logical(1), "GDSArray")))
        return(wmsg("'assays(x)' must be GDSArray object"))
    
    ## ## gdsfile correlated with assay data
    ## if(is.character(gdsfile(x)[[1]]))
    ##     return(wmsg("There should be an gds file correlated with 'x'"))
    TRUE
}

#' @import S4Vectors 
## setValidity2("VariantExperiment", .validate_VariantExperiment)

###--------------------
### getter and setter
###--------------------
#' @export gdsfile
#' @rdname VariantExperiment-class
setMethod("gdsfile", "VariantExperiment", function(object)
    vapply(assays(object), gdsfile, character(1)))

#' @export "gdsfile<-"
#' @param value the new gds file path for VariantExperiment object.
#' @rdname Variantexperiment-class
setReplaceMethod("gdsfile", "VariantExperiment", function(object, value) {
    new_filepath <- tools::file_path_as_absolute(value)
    assays(object) <- lapply(assays(object), function(assay)
        BiocGenerics:::replaceSlots(seed(assay), file=value, check=FALSE))
    object
})


