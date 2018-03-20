###----------  
## lazyList
###----------
.lazyList <- setClass(
    "lazyList",
    slots = c(
        indexes = "list",
        has_index = "integer"
    )
)

###---------------------------------------
## Utility functions for .lazyList
###---------------------------------------

.update_index <- function(lazyList, j, value)
{
    for (i in seq_along(lazyList@indexes)) {
        index <- lazyList@indexes[[i]]
        if (identical(value, index)) {
            lazyList@has_index[j] <- i
            return(lazyList)
        }
    }
    new_index_id <- length(lazyList@indexes) + 1L
    lazyList@has_index[j] <- new_index_id
    index <- value
    lazyList@indexes[[new_index_id]] <- index
    lazyList
}

.update_row <- function(lazyList, value)
{
    lazyList@indexes <- lapply(lazyList@indexes, function(index) {
        if (is(index, "list")) {
            if (is.null(index[[1]])) {
                index[[1]] <- value
            } else {
            if (length(value) > length(index[[1]]) )
                stop("the subscripts are out of bound")
            index[[1]] <- index[[1]][value]
            }
        }
        index
    })
    lazyList
}

.lazyIndex_inuse <- function(lazyIndex)
{
    orig <- seq_len(length(lazyIndex@indexes)) 
    index_inuse <- orig %in% unique(lazyIndex@has_index)
    old <- orig[index_inuse]
    lazyIndex@indexes <- lazyIndex@indexes[index_inuse]
    new <- seq_len(length(lazyIndex@indexes))

    lazyIndex@has_index <- new[match(lazyIndex@has_index, old)]
    lazyIndex
}

###------------------
## DelayedDataFrame
###------------------ 

#' DelayedDataFrame-class
#' @name DelayedDataFrame
#' @exportClass DelayedDataFrame
#' @importFrom methods as initialize is new "slot<-"
#' @aliases DelayedDataFrame-class
#' @description The \code{DelayedDataFrame} class extends the
#'     \code{DataFrame} class and supports the storage of any type of
#'     object (with ‘length’ and ‘[’ methods) as columns.
#' @rdname DelayedDataFrame-class
#' @details The \code{DelayedDataFrame} inherits from \code{DataFrame}
#'     and behaves very similarily in terms of construction,
#'     subsetting, splitting, combining, etc. The most notable
#'     exception is that The additional slot of \code{lazyIndex},
#'     enables \code{DelayedArray} (with different back-ends) columns
#'     to share indexes when possible.
## refer ?DataFrame 
## #' \code{Constructor}: \code{DelayedDataFrame(..., row.names = NULL,
## #' check.names = TRUE)} constructs a ‘DelayedDataFrame’ in similar
## #' fashion to ‘DataFrame’. Each argument in \code{...} is coerced to a
## #' ‘DelayedDataFrame’ and combined column-wise.

.DelayedDataFrame <- setClass(
    "DelayedDataFrame",
    contains = "DataFrame",
    slots = c(lazyIndex = "lazyList")
)

###---------------------------------------
## Utility functions for DelayedDataFrame
###---------------------------------------
.get_index <- function(x, j)
{
    j <- x@lazyIndex@has_index[[j]]
    x@lazyIndex@indexes[[j]]
}

#' @description \code{update_lazyIndex}: make sure the indexes are
#'     consistent with columns and being used. update the
#'     \code{indexes} slot and \code{@has_index} slot in lazyIndex.
#' @rdname DelayedDataFrame-class
update_lazyIndex <- function(from)
{
    lazyIndex <- from@lazyIndex
    delayedClass <- vapply(from, is, logical(1), "DelayedArray")

    for (i in seq_len(length(from))) {
        if (delayedClass[i]) {
            index <- from[[i]]@index
        } else {
            index <- NULL
        }
        lazyIndex <- .update_index(lazyIndex, i, index)
    }
    ## remove any index if not in use
    from@lazyIndex <- .lazyIndex_inuse(lazyIndex)
    ## reorder indexes
    from@lazyIndex <- from@lazyIndex[TRUE]
    ## assign new indexes into DelayedArray columns.
    for (i in seq_len(length(from))) {
        if(delayedClass[i])
            from[[i]]@index <- .get_index(from, i)
    }
    from
}

###-----------------
## subsetting
###----------------
.extractROWS_DelayedDataFrame <- function(x, i)
{
    ## x <- update_lazyIndex(x)  
    i <- normalizeSingleBracketSubscript(
        i, x, exact = FALSE, allow.NAs = TRUE, as.NSBS = FALSE)
    x@lazyIndex <- .update_row(x@lazyIndex, i)
    slot(x, "listData") <- lapply(
        setNames(seq_len(ncol(x)), names(x)), function(j)
        {
            if(!is(x[[j]], "DelayedArray"))  ## non-delayed_ops object.
                extractROWS(x[[j]], i)
            else {                    ## DelayedArray objects.
                a <- x[[j]]
                a@index <- .get_index(x, j)
                a
            }
        })
    slot(x, "nrows", check = FALSE) <- length(i)
    if (!is.null(rownames(x))) {
        slot(x, "rownames", check = FALSE) <-
            make.unique(extractROWS(rownames(x), i))
    }
    x
}
#' @exportMethod extractROWS
#' @aliases extractROWS,DelayedDataFrame-method
#' @rdname DelayedDataFrame-class
setMethod("extractROWS", "DelayedDataFrame", .extractROWS_DelayedDataFrame)

setMethod("[", c("lazyList", "ANY", "missing", "ANY"),
          function(x, i, j, ..., drop = TRUE)
          {
              has_index <- x@has_index[i]
              uhas_index <- unique(has_index)
              has_index <- match(has_index, uhas_index)
              indexes <- x@indexes[uhas_index]
              .lazyList(indexes = indexes, has_index = has_index)
          }
          )

#' @importFrom methods callNextMethod
#' @exportMethod [
#' @aliases [,DelayedDataFrame-method
#' @rdname DelayedDataFrame-class
#' @param x input
#' @param i row subscript
#' @param j col subscript
#' @param drop if drop with reduced dimension, default is TRUE.
#' @param row.names rownames
#' @param check.names if check names.
#' @param ... other arguments to pass.

setMethod("[", c("DelayedDataFrame", "ANY", "ANY", "ANY"),
    function(x, i, j, ..., drop = TRUE)
    {
        list_style_subsetting <- (nargs() - (!missing(drop))) < 3L
        if (list_style_subsetting || !missing(j)) {
            if (list_style_subsetting) {
                if (!missing(drop)) 
                    warning("'drop' argument ignored by list-style subsetting")
                if (missing(i)) 
                    return(x)
                j <- i
            }
            
            if(!missing(j) && !is(j, "IntegerRanges")) {
                xstub <- setNames(seq_along(x), names(x))
                j <- normalizeSingleBracketSubscript(j, xstub)
            }
            x@lazyIndex <- x@lazyIndex[j]
        }
        callNextMethod()
    })

setMethod(
    "concatenateObjects", "DelayedDataFrame",
    function(x, objects = list(), use.names = TRUE, ignore.mcols = FALSE, check = TRUE)
{
    y <- as(objects, "DelayedDataFrame")
    x@lazyIndex <- concatenateObjects(x@lazyIndex, y@lazyIndex)
    callNextMethod()
})

setMethod("concatenateObjects", "lazyList",
    function(x, objects, use.names = TRUE, ignore.mcols = FALSE, check = TRUE) 
{
    ## objects <- unlist(objects)
    if (!is(objects, "lazyList")) 
        stop("'objects' must be a \"lazyList\" class")
    if (!isTRUEorFALSE(use.names)) 
        stop("'use.names' must be TRUE or FALSE")

    for (i in seq_len(length(objects@has_index))) {
        new_has_index <- objects@has_index[i]
        new_index <- objects@indexes[[new_has_index]]
        x@has_index <- c(x@has_index, new_has_index)
        x <- .update_index(x, length(x@has_index), value=new_index)
    }
    x
})

###-------------
## constructor
###-------------

#' @export
#' @rdname DelayedDataFrame-class
DelayedDataFrame <- function(..., row.names=NULL, check.names=TRUE)
{
    df <- DataFrame(..., row.names=row.names, check.names=check.names)
    as(df, "DelayedDataFrame")
}

###-------------
### Coercion
###-------------

## DelayedDataFrame has inherited from DataFrame, so it inherits
## coercion methods of DataFrame to matrix/data.frame/list (as.matrix,
## as.list, as.data.frame/as(x, "data.frame/list")). Will only need to
## define set("ANY", "DelayedDataFrame").

#' @name coerce
#' @exportMethod coerce
#' @aliases coerce,DataFrame,DelayedDataFrame-method
#' @rdname DelayedDataFrame-class
#' @param from a \code{DataFrame} object
setAs("DataFrame", "DelayedDataFrame", function(from){
    ## initial lazyIndex
    lazyIndex <- .lazyList(
        indexes = vector("list", 1),
        has_index = rep(1L, length(from)))
    ans <- .DelayedDataFrame(from, lazyIndex = lazyIndex)
    update_lazyIndex(ans)
})

###
setAs("ANY", "DelayedDataFrame", function(from){
    df <- as(from, "DataFrame")
    as(df, "DelayedDataFrame")
})

###--------------
## slot setters
###--------------

## replace method for lazyIndex(DDF)
setGeneric(
    "lazyIndex<-",
    function(x, value) standardGeneric("lazyIndex<-"),
    signature="x")

#' @exportMethod "lazyIndex<-"
#' @rdname DelayedDataFrame-class
#' @description the setter for the \code{lazyIndex} slot of \code{DelayedDataFrame} object.
#' @return the new value of \code{lazyIndex} slot for \code{DelayedDataFrame} object.
setReplaceMethod( "lazyIndex", "DelayedDataFrame", function(x, value) {
    BiocGenerics:::replaceSlots(x, lazyIndex=value, check=FALSE)
})

###-----------------
## validity check
###----------------
.validate_DelayedDataFrame <- function(x)
{
    ## indexes length must be same
    indexes <- x@lazyIndex@indexes
    indexLength <- vapply(indexes, function (index) {
        if (is.list(index)) length(index[[1]])
        else length(index)
    }, integer(1))
    
    if (length(unique(indexLength)) > 1 && all(unique(indexLength) > 0))
        return(wmsg("'x@lazyIndex' must be of same length or 'NULL'"))

    ## has_index must have same length of ncol(x)
    if(length(x@lazyIndex@has_index) != ncol(x))
        return(wmsg("'x@has_index' must be of same length of 'ncols(x)'"))
    TRUE
}

#' @importFrom S4Vectors setValidity2
setValidity2("DelayedDataFrame", .validate_DelayedDataFrame)

