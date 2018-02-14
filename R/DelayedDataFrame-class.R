## library(DelayedArray)

.lazyList <- setClass(
    "lazyList",
    slots = c(
        indexes = "list",
        has_index = "integer"
    )
)

update_index <- function(lazyList, j, value) {
    for (i in seq_along(lazyList@indexes)){
        index <- lazyList@indexes[[i]]$index[[1]]
        if(!is.null(index)){
            if (all(value == index)) {
                lazyList@has_index[j] = i
                return(lazyList)
            }
        }
    }
    ## clone (refClass clone, different from ordinary S4 class)
    new_index_id <- length(lazyList@indexes) + 1L
    lazyList@has_index[j] = new_index_id
    index = IndexList(index=list(value))
    lazyList@indexes[[new_index_id]] = index
    lazyList
}

update_row <- function(lazyList, value) {
    lazyList@indexes <- lapply(lazyList@indexes, function(index) {
        index <- DelayedArray:::.clone(index)
        if (is.null(index$index[[1]]))
            index$index[[1]] <- value
        else {
            if (length(value) > length(index$index[[1]]) )
                stop("the subscripts are out of bound")
            index$index[[1]] <- index$index[[1]][value]
        }
        index
    })
    lazyList
}

.DelayedDataFrame = setClass(
    "DelayedDataFrame",
    contains = "DataFrame",
    slots = c(lazyIndex = "lazyList")
)

DelayedDataFrame <- function(x) {
    ## check column seed dimension
    ## (DelayedDF only works for DelayedArray with different backend)
    ndim <- length(dim(x[[1]]))
    if (ndim == 1)
        indexList = IndexList(index=list(NULL))
    else if (ndim == 2)
        indexList = IndexList(index = list(NULL, NULL))
    for (i in seq_along(x))
        x[[i]]@index = indexList
    lazyIndex <- .lazyList(
        indexes = list(indexList),
        has_index = rep(1L, length(x))
    )
    .DelayedDataFrame(x, lazyIndex = lazyIndex)
}

.validate_delayedDF <- function(x)
{
    ## indexes class must be "IndexList"
    indexes <- x@lazyIndex@indexes
    indexClass <- vapply(indexes,
                         function(index) is(index, "IndexList"),
                         logical(1))
    if (!all(indexClass))
        return(wmsg("'x@lazyIndex' must be a list of 'IndexList'"))

    ## indexes length must be same
    indexLength <- vapply(indexes,
                          function(index) length(index$index[[1]]),
                          integer(1))
    if (length(unique(indexLength)) > 1 & all(unique(indexLength) > 0))
        return(wmsg("'x@lazyIndex' must be of same length or 'NULL'"))

    ## col classes (must be DelayedArray with different back-end)
    colClass <- vapply(x,
                       function(col) is(col, "DelayedArray"),
                       logical(1))
    if (!all(colClass))
        return(wmsg("`DelayedDataFrame` columns must be",
                    "inherited from 'DelayedArray'"))

    TRUE
}

#' @importFrom S4Vectors setValidity2
setValidity2("DelayedDataFrame", .validate_delayedDF)

.get_index <- function(x, j) {
    j <- x@lazyIndex@has_index[[j]]
    x@lazyIndex@indexes[[j]]
}

#' @exportMethod extractROWS
.extractROWS_DelayedDataFrame <- function(x, i)
{
    i <- normalizeSingleBracketSubscript(i, x, exact = FALSE, 
                                         allow.NAs = TRUE
                                         , as.NSBS = TRUE
                                         )
    x@lazyIndex <- update_row(x@lazyIndex, i@subscript)
    slot(x, "nrows", check = FALSE) <- length(i)
    if (!is.null(rownames(x))) {
        slot(x, "rownames", check = FALSE) <- make.unique(extractROWS(rownames(x), 
            i))
    }
    ## copy updated lazy index to each column
    for (j in seq_along(x))
        DelayedArray:::.index(x[[j]]) <- .get_index(x, j)
    x
}
setMethod("extractROWS", "DelayedDataFrame", .extractROWS_DelayedDataFrame)


#################### VariantExperiment example  ####################################
## file <- SeqArray::seqExampleFileName("gds")
## library(VariantExperiment)
## se <- makeSummarizedExperimentFromGDS(file)
## rd <- rowData(se)
## rd[[1]]@index  ## ordinary list.
## ddf <- DelayedDataFrame(rd)
## ddf@lazyIndex
## a <- ddf[c(1,3,5), c(TRUE, FALSE), drop=FALSE]

## ############################################ bk ####################################
## object.size(rd[TRUE, ])
## ## 1407656 bytes

## se1 <- makeSummarizedExperimentFromGDS(file, rowDataOnDisk=FALSE)
## object.size(rowData(se1)[TRUE,])
## ## 335152 bytes

## object.size(as.list(rowData(se1)[["ALT"]]))
## ## 7300816 bytes


