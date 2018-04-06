###----------
## LazyIndex
###----------

.LazyIndex <- setClass(
    "LazyIndex",
    contains = "SimpleList",
    slots = c(
        index = "integer"
    )
)

LazyIndex <-
    function(listData = list(), index = integer())
{
    .LazyIndex(listData, index = index)
}

.listData <- function(x)
    x@listData

.index <- function(x)
    x@index

.validate_LazyIndex <- function(x)
{
    ## indexes length must be same
    indexes <- .listData(x)
    indexLength <- lengths(indexes)
    uniqLen <- unique(indexLength)
    if (length(uniqLen) == 1)
        return(TRUE)
    if (length(uniqLen[uniqLen != 0]) > 1)
        return(wmsg("'.index(x)' must be of same length or 'NULL'"))
    TRUE
}

#' @importFrom S4Vectors setValidity2
setValidity2("LazyIndex", .validate_LazyIndex)

setMethod("concatenateObjects", "LazyIndex",
          function(x, objects=list(), use.names = TRUE,
                   ignore.mcols = FALSE, check = TRUE)
{
    listData <- c(.listData(x), lapply(objects, slot, "listData"))
    indexes <- c(list(.index(x)), lapply(objects, slot, "index"))

    ## index offsets
    offsets <- head(cumsum(c(0L, lengths(indexes))), -1L)
    index <- unlist(indexes) + rep(offsets, lengths(indexes))

    .lazyIndex_compose(listData, index)
})

setMethod("[", c("LazyIndex", "ANY"),
    function(x, i, j, ..., drop = TRUE)
{
    listData <- .listData(x)
    index <- .index(x)

    ## browser(); browser()
    if (!isTRUEorFALSE(drop)) 
        stop("'drop' must be TRUE or FALSE")
    if (length(list(...)) > 0L) 
        warning("parameters in '...' not supported")
    list_style_subsetting <- (nargs() - (!missing(drop))) < 3L
    if (list_style_subsetting || !missing(j)) {
        if (list_style_subsetting) {
            if (!missing(drop)) 
                warning("'drop' argument ignored by list-style subsetting")
            if (missing(i)) 
                return(x)
            j <- i
        }
        if (!is(j, "IntegerRanges")) {
            xstub <- setNames(seq_along(index), seq_along(index))
            j <- normalizeSingleBracketSubscript(j, xstub)
        }
        index <- extractROWS(index, j)
        x <- .lazyIndex_compose(listData, index)
        if (list_style_subsetting)
            return(x)
    }
    if (!missing(i)) {
        x <- .update_row(x, i)
        ## new_listData <- extractROWS(listData, i)
        ## LazyIndex(new_listData, index)
    }
    x
})

###---------------------------------------
## Utility functions for .LazyIndex
###---------------------------------------

.lazyIndex_compose <-
    function(listData0, index0)
{
    ## 1. check if any duplicate in @listData. Modify @index
    listData1 <- unique(listData0)
    index1 <- match(listData0, listData1)[index0]

    ## 2. check if all @listData in use in @index. remove if not.
    listData2 <- listData1[seq_along(listData1) %in% index1]
    index2 <- match(listData1, listData2)[index1]

    ## 3. reorder indexes (index for 1st columns in @listData[[1]])
    index3 <- unique(index2)
    listData <- listData2[index3]
    index <- match(index2, index3)

    LazyIndex(listData, index=index)
}

.update_index <- function(lazyList, j, value)
{
    listData <- c(.listData(lazyList), list(value))
    index <- .index(lazyList)
    index[j] <- length(listData)
    .lazyIndex_compose(listData, index)
}

.update_row <- function(lazyList, i)
{
    ## browser(); browser()
    listData <- .listData(lazyList)
    isNull <- vapply(listData, is.null, logical(1))
    if (any(lengths(listData[!isNull]) < length(i)))
        stop("subscripts are out of bound")

    listData[isNull] <- list(i)
    listData[!isNull] <- lapply(listData[!isNull], `[`, i = i)

    .lazyIndex_compose(listData, .index(lazyList))
}

setMethod("show", "LazyIndex", function(object)
{
    lo <- length(object)
    cat(classNameForDisplay(object), " of length ", lo, "\n",
        sep = "")
    ## cat("Indexes: ", "\n", sep="")
    print(.listData(object))
    cat("index of each column: ", "\n", sep="")
    print(.index(object))
})
