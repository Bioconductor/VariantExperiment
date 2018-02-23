#' @exportClass VariantExperiment
.VariantExperiment <- setClass(
    "VariantExperiment",
    contains="SummarizedExperiment",
    ## slots = c(gdsfile = "character")
    ## slots=c(gdsfile="character",
    ##         x="integer")
)

## constructor
#' @export VariantExperiment
VariantExperiment <- function(assays, rowRanges=GRangesList(), colData=DelayedDataFrame(), metadata=list(), ...)
{
    ## se <- SummarizedExperiment:::.new_RangedSummarizedExperiment(
    se <- SummarizedExperiment(
        rowRanges=rowRanges,
        colData=colData,
        assays=assays,
        ...)
    ## rowData(se) <- as(rowData(se), "DelayedDataFrame")
    ## colData(se) <- as(colData(se), "DelayedDataFrame")
    ## mcols(rowRanges(se)) <- rowData(se)
    se
    ## as(se, "VariantExperiment")
}

## setAs("SummarizedExperiment", "VariantExperiment", function(from)
## {
##     gf <- gdsfile(from)[[1]]
##     .VariantExperiment(from, gdsfile = gf)
## })
###----------------------------------------------------------------
### Coercion
###

## 

.validate_VariantExperiment <- function(x)
{
    ## DelayedDataFrame
    if(!all(is(rowData(x), "DelayedDataFrame"), is(colData(x), "DelayedDataFrame")))
        return(wmsg("'rowData(x)' and 'colData(x)' must be DelayedDataFrame object"))

    ## GDSArray for assay data.
    ## if(!all(vapply(assays(x), function(x) is(x, "GDSArray")), logical(1)))
    if(!all(vapply(assays(x), is, logical(1), "GDSArray")))
        return(wmsg("'assays(x)' must be GDSArray object"))
    
    ## gdsfile correlated with assay data
    if(is.character(gdsfile(x)[[1]]))
        return(wmsg("There should be an gds file correlated with 'x'"))
    TRUE
}

#' @importFrom S4Vectors setValidity2
setValidity2("VariantExperiment", .validate_VariantExperiment)


