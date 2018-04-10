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
    slots = c(lazyIndex = "LazyIndex")
)

###---------------------------------------
## Utility functions for DelayedDataFrame
###---------------------------------------
.get_index <- function(x, j)
{
    j <- .index(lazyIndex(x))[[j]]
    .listData(lazyIndex(x))[[j]]
}

###-------------
## constructor
###-------------

#' @export
#' @rdname DelayedDataFrame-class
DelayedDataFrame <- function(..., row.names=NULL, check.names=TRUE)
{
    ## no-op for DelayedDataFrame input
    if (length(list(...)) == 1L && is(list(...)[[1]], "DelayedDataFrame"))
        return(list(...)[[1]])
    df <- DataFrame(..., row.names=row.names, check.names=check.names)
    as(df, "DelayedDataFrame")
}

###-------------
## accessor
###-------------

setGeneric("lazyIndex", function(x) standardGeneric("lazyIndex"), signature="x")

setMethod("lazyIndex", "DelayedDataFrame", function(x) x@lazyIndex)

###-------------
## methods
###-------------

setMethod("getListElement", "DelayedDataFrame", function(x, i, exact=TRUE)
{
    i2 <- normalizeDoubleBracketSubscript(
        i, x, exact = exact,
        allow.NA = TRUE,
        allow.nomatch = TRUE)
    if (is.na(i2)) 
        return(NULL)
    index <- .get_index(x, i2)
    elt <- x@listData[[i2]]
    if (!is.null(index))
        elt <- extractROWS(elt, index)
    elt
})

## "as.list" function is called in lapply("DelayedDataFrame", ) and names("DelayedDataFrame")...
setMethod("as.list", "DelayedDataFrame", function(x, use.names=TRUE)  
{
    ans <- lapply(seq_along(x), function(j) x[[j]])
    if (use.names)
        names(ans) <- names(x)
    ans
})

setMethod("names", "DelayedDataFrame", function(x)
{
    names(x@listData)
})

## FIXME: GDSArray coerced into character vector...
## setMethod("rbind", "DelayedDataFrame", function(..., deparse.level=1)
## {
##     DelayedDataFrame(callNextMethod())
## })

setMethod("cbind", "DelayedDataFrame", function(..., deparse.level=1)
{
    df <- callNextMethod()
    DelayedDataFrame(df)
})


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
setAs("DataFrame", "DelayedDataFrame", function(from)
{
    if (identical(dim(from), c(0L, 0L))) {
    lazyIndex <- .LazyIndex()
    } else {     
        lazyIndex <- .LazyIndex(vector("list", 1), index=rep(1L, length(from)))
    }
    .DelayedDataFrame(from, lazyIndex = lazyIndex)
})

setAs("DelayedDataFrame", "DataFrame", function(from)
{
    new_listData <- as.list(from)
    ans <- S4Vectors:::new_DataFrame(listData=new_listData, nrows=nrow(from))
    if (!is.null(rownames(from))) {
       ans@rownames <- rownames(from)
    }
    ans
})

## setAs("DelayedDataFrame", "DelayedDataFrame", function(from) from)  ## no-op

###
setAs("ANY", "DelayedDataFrame", function(from){
    df <- as(from, "DataFrame")
    as(df, "DelayedDataFrame")
})

###-----------------
## validity check
###----------------
.validate_DelayedDataFrame <- function(x)
{
    .validate_LazyIndex(lazyIndex(x))
    ## @index must have same length of ncol(x)
    if(length(.index(lazyIndex(x))) != ncol(x))
        return(wmsg("'.index(x)' must be of same length of 'ncols(x)'"))
    TRUE
}

#' @importFrom S4Vectors setValidity2
setValidity2("DelayedDataFrame", .validate_DelayedDataFrame)

###-----------------
## subsetting
###----------------
.extractROWS_DelayedDataFrame <- function(x, i)
{
    i <- normalizeSingleBracketSubscript(
        i, x, exact = FALSE, allow.NAs = TRUE, as.NSBS = FALSE)
    ## lazyIndex(x) <- .update_row(lazyIndex(x), i)
    lazyIndex(x) <- lazyIndex(x)[i,]
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

setReplaceMethod(
    "[", c("DelayedDataFrame", "ANY"),
    function(x, i, j, ..., value)
{
    xstub <- setNames(seq_along(x), names(x))
    if (missing(j)) {
        i <- normalizeSingleBracketSubscript(i, xstub)
        lazyIndex(x) <- .update_index(lazyIndex(x), i, NULL)
    } else {
        j <- normalizeSingleBracketSubscript(j, xstub)
        x@listData[j] <- lapply(j, function(j, x) x[[j]], x)
        lazyIndex(x) <- .update_index(lazyIndex(x), j, NULL)
    }
    callNextMethod()
})

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
            lazyIndex(x) <- lazyIndex(x)[j]
        }
        callNextMethod()
    })

## constructing a new DelayedDataFrame
setMethod(
    "concatenateObjects", "DelayedDataFrame",
    function(x, objects = list(), use.names = TRUE, ignore.mcols = FALSE, check = TRUE)
{
    for (i in seq_len(length(objects))) {
        y <- as(objects[[i]], "DelayedDataFrame")
        lazyIndex(x) <- c(lazyIndex(x), lazyIndex(y))
    }
    callNextMethod()
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

