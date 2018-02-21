#' @exportClass VariantExperiment 
setClass("VariantExperiment",
         contains="RangedSummarizedExperiment",
         slots = c(gdsfile = "character"))

## constructor
VariantExperiment <- function(assays, rowRanges=GRangesList(), colData=DelayedDataFrame(), metadata=list(), ...)
{
    ## se <- SummarizedExperiment:::.new_RangedSummarizedExperiment(
    ##                                  rowRanges=rowRanges,
    ##                                  colData=colData,
    ##                                  assays=assays,
    ##                                  ...)
    gf <- gdsfile(se)
    new("VariantExperiment",
        rowRanges = rowRanges,
        colData = colData,
        assays = assays,
        metadata = metadata,
        gdsfile = gf)
}

.validate_VariantExperiment <- function(x)
{
    ## DelayedDataFrame
    if(!all(is(rowData(x), "DelayedDataFrame"), is(colData(x), "DelayedDataFrame")))
        return(wmsg("'rowData(x)' and 'colData(x)' must be DelayedDataFrame object"))

    ## GDSArray for assay data.
    if(!all(vapply(assays(x), function(x) is(x, "GDSArray")), logical(1)))
        return(wmsg("'assays(x)' must be GDSArray object"))

    ## gdsfile correlated with assay data
    if(is.character(gdsfile(x)[[1]]))
        return(wmsg("There should be an gds file correlated with 'x'"))
    TRUE
}

#' @importFrom S4Vectors setValidity2
setValidity2("VariantExperiment", .validate_VariantExperiment)


