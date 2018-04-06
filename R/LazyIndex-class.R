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

setMethod("concatenateObjects", "LazyIndex",
          function(x, objects=list(), use.names = TRUE,
                   ignore.mcols = FALSE, check = TRUE)
{
    ## browser(); browser()
    listData <- c(x@listData, lapply(objects, slot, "listData"))
    indexes <- c(list(x@index), lapply(objects, slot, "index"))

    ## index offsets
    offsets <- cumsum(c(0L, head(lengths(indexes), -1)))
    index <- unlist(indexes) + rep(offsets, lengths(indexes))

    .lazyIndex_compose(listData, index)
})

setMethod("[", c("LazyIndex", "ANY", "missing", "ANY"),
    function(x, i, j, ..., drop = TRUE)
{
    LazyIndex(x@listData, index = x@index[i])
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
    listData <- c(lazyList@listData, list(value))
    index <- lazyList@index
    index[j] <- length(listData)
    .lazyIndex_compose(listData, index)
}

.update_row <- function(lazyList, value)
{
    ## browser(); browser()
    listData <- lapply(lazyList@listData, function(index) {
        if (is.null(index)) {
            index <- value
        } else {
            if (length(value) > length(index) )
                stop("the subscripts are out of bound")
            index <- index[value]
        }
        index
    })
    index <- lazyList@index
    .lazyIndex_compose(listData, index)
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
