setMethod("[", c("DelayedDataFrame", "ANY", "ANY", "ANY"),
    function(x, i, j, ..., drop = TRUE)
{
    lazyIndex(x) <- lazyIndex[i,j]
    ##
    x@listData <- x@listData[seq_along(.index(lazyIndex(x)))]
    if (!is.null(length(lazyIndex(x)))) 
        slot(x, "nrows", check=FALSE) <- length(lazyIndex(x)) 
    if (!is.null(rownames(x))) {
        slot(x, "rownames", check = FALSE) <-
            make.unique(extractROWS(rownames(x), i))  ## FIXME
    }
