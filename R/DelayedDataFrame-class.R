## library(DelayedArray)

.lazyList <- setClass(
    "lazyList",
    slots = c(
        indexes = "list",
        has_index = "integer"
    )
)

## "update_index" not actually used in constructing "DelayedDataFrame"...
.update_index <- function(lazyList, j, value) {
    for (i in seq_along(lazyList@indexes)){
        index <- lazyList@indexes[[i]]$index
        if (identical(value, index)){
            lazyList@has_index[j] <- i
            return(lazyList)
        }
    }
    ## clone (refClass clone, different from ordinary S4 class)
    new_index_id <- length(lazyList@indexes) + 1L
    lazyList@has_index[j] = new_index_id
    index = IndexList(index=value)
    lazyList@indexes[[new_index_id]] = index
    lazyList
}

## FIXME: how to simplify these code? 
.update_lazyIndex <- function(lazyIndex, from){
    if(!is(from, "DelayedDataFrame")){
        for (i in seq_len(length(from))){
            index <- from[[i]]@index
            lazyIndex <- .update_index(lazyIndex, i, index)
        }
    }
    index_inuse <- seq_len(length(lazyIndex@indexes)) %in% unique(lazyIndex@has_index)
    lazyIndex@indexes <- lazyIndex@indexes[index_inuse]
    
    for (i in seq_len(length(from))){
        index <- from[[i]]@index
        lazyIndex <- .update_index(lazyIndex, i, index)
    }
    lazyIndex
}
    
.update_row <- function(lazyList, value) {
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

.get_index <- function(x, j) {
    j <- x@lazyIndex@has_index[[j]]
    x@lazyIndex@indexes[[j]]
}

#' @exportClass DelayedDataFrame
.DelayedDataFrame = setClass(
    "DelayedDataFrame",
    contains = "DataFrame",
    slots = c(lazyIndex = "lazyList")
)

#' @exportMethod extractROWS
.extractROWS_DelayedDataFrame <- function(x, i)
{
    i <- normalizeSingleBracketSubscript(i, x, exact = FALSE, 
                                         allow.NAs = TRUE
                                         , as.NSBS = FALSE
                                         )
    x@lazyIndex <- .update_row(x@lazyIndex, i)
    slot(x, "nrows", check = FALSE) <- length(i)
    if (!is.null(rownames(x))) {
        slot(x, "rownames", check = FALSE) <- make.unique(extractROWS(rownames(x), 
            i))
    }
    ## copy updated lazy index to each column
    for (j in seq_along(x))
        DelayedArray:::.index(x[[j]]) <- .get_index(x, j)$index
    x
}
setMethod("extractROWS", "DelayedDataFrame", .extractROWS_DelayedDataFrame)

## DelayedDataFrame constructor

DelayedDataFrame <- function(..., row.names=NULL, check.names=TRUE){
    df <- DataFrame(..., row.names=row.names, check.names=check.names)
    as(df, "DelayedDataFrame")
}

## validity check
.validate_DelayedDataFrame <- function(x)
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
setValidity2("DelayedDataFrame", .validate_DelayedDataFrame)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

## DelayedDataFrame has inherited from DataFrame, so it inherits coercion methods of DataFrame to matrix/data.frame/list (as.matrix, as.list, as.data.frame/as(x, "data.frame/list")). Will only need to define set("ANY", "DelayedDataFrame"). 

## .as.DelayedDataFrame.DataFrame <- function(x){
setAs("DataFrame", "DelayedDataFrame", function(from){
    delayed_ops <- vapply(from, is, logical(1), "DelayedArray")
    
    ## (DelayedDF only works for DelayedArray with different backend)
    if(any(!delayed_ops)){
        if (sum(!delayed_ops) == 1){
            words <- c("column", "is", "object")
        } else {
            words <- c("columns", "are", "objects")
        }
        stop("\n", "All columns should be inherited from DelayedArray. ",
             "\n", "The ", words[1], " of '",
             paste(names(from)[!delayed_ops], collapse=", "), "' ",
             words[2], " not DelayedArray ",
             words[3], ". \n")
    }
    
    ## 
    lazyIndex <- .lazyList(indexes = list(IndexList(index=vector("list", 1))),
                           has_index = rep(1L, length(from)))
    lazyIndex <- .update_lazyIndex(lazyIndex, from)

    .DelayedDataFrame(from, lazyIndex = lazyIndex)
})

## setMethod("DelayedDataFrame", "DataFrame", .as.DelayedDataFrame.DataFrame)

###
setAs("ANY", "DelayedDataFrame", function(from){
    df <- as(from, "DataFrame")
    as(df, "DelayedDataFrame")
})


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


