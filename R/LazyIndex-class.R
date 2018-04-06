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
    function(lazyData = list(), index = integer())
{
    ## order index c(4, 2, 1) --> 1, 2, 3
    uindex <- unique(index)
    index <- match(index, uindex)
    lazyData <- lazyData[uindex]
    .LazyIndex(lazyData, index = index)
}

setMethod("show", "LazyIndex", function(object)
{
    lo <- length(object)
    cat(classNameForDisplay(object), " of length ", lo, "\n",
        sep = "")
    ## cat("Indexes: ", "\n", sep="")
    print(object@listData)
    cat("index of each column: ", "\n", sep="")
    print(object@index)
})

.validate_LazyIndex <- function(x)
{
    ## indexes length must be same
    indexes <- x@listData
    indexLength <- lengths(indexes)
    uniqLen <- unique(indexLength)
    if (length(uniqLen) == 1)
        return(TRUE)
    if (length(uniqLen[uniqLen != 0]) > 1)
        return(wmsg("'x@lazyIndex' must be of same length or 'NULL'"))
    TRUE
}
#' @importFrom S4Vectors setValidity2
setValidity2("LazyIndex", .validate_LazyIndex)

setMethod("[", c("LazyIndex", "ANY", "missing", "ANY"),
    function(x, i, j, ..., drop = TRUE)
{
    LazyIndex(x@listData, index = x@index[i])
})

###---------------------------------------
## Utility functions for .LazyIndex
###---------------------------------------

.lazyIndex_inuse <- function(lazyList)
{
    ## browser(); browser()
    ## 1. check if any duplicate in @listData. Modify @index
    ## correspondingly.
    listData <- unique(lazyList@listData)
    index <- match(lazyList@listData, listData)[lazyList@index]

    ## 2. check if all @listData in use in @index. remove if not.
    listData_index <- seq_along(listData)
    inUse <- listData_index %in% index
    listData <- listData[inUse]
    index <- cumsum(inUse)[index]

    ## 3. reorder indexes (index for 1st columns in @listData[[1]])
    LazyIndex(listData, index=index)
}

.update_index <- function(lazyList, j, value)
{
    lazyList@listData <- c(lazyList@listData, list(value))
    lazyList@index[j] <- length(lazyList@listData)
    .lazyIndex_inuse(lazyList)
}

.update_row <- function(lazyList, value)
{
    ## browser(); browser()
    lazyList@listData <- lapply(lazyList@listData, function(index) {
        if (is.null(index)) {
            index <- value
        } else {
            if (length(value) > length(index) )
                stop("the subscripts are out of bound")
            index <- index[value]
        }
        index
    })
    .lazyIndex_inuse(lazyList)
}
