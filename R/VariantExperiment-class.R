#' @importClassesFrom SummarizedExperiment SummarizedExperiment RangedSummarizedExperiment 
#' @exportClass VariantExperiment
.VariantExperiment <- setClass(
    "VariantExperiment",
    contains="RangedSummarizedExperiment",
    ## slots = c(gdsfile = "character")
    ## slots=c(gdsfile="character",
    ##         x="integer")
)

###--------------
### constructor
###--------------
#' @export VariantExperiment
VariantExperiment <- function(assays, rowRanges=GRangesList(), colData=DelayedDataFrame(), metadata=list())
{
    elementMetadata <- S4Vectors:::make_zero_col_DataFrame(length(rowRanges))
    if (!is(assays, "Assays"))
        assays <- Assays(assays)
    new("VariantExperiment", rowRanges=rowRanges,
        colData=colData,
        assays=assays,
        elementMetadata=elementMetadata,
        metadata=as.list(metadata))
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
setValidity2("VariantExperiment", .validate_VariantExperiment)

###--------------------
### getter and setter
###--------------------
setMethod("gdsfile", "VariantExperiment", function(x)
    vapply(assays(x), gdsfile, character(1)))

setReplaceMethod("gdsfile", "VariantExperiment", function(x, value) {
    new_filepath <- tools::file_path_as_absolute(value)
    assays(x) <- lapply(assays(x), function(assay)
        BiocGenerics:::replaceSlots(seed(assay), file=value, check=FALSE))
    x
})


