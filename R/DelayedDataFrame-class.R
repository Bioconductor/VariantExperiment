###----------  
## lazyList
###----------
.LazyList <- setClass(
    "LazyList",
    contains = "SimpleList",
    slots = c(
        index = "integer"
    )
)

setMethod("show", "LazyList", function(object)
{
    lo <- length(object)
    cat(classNameForDisplay(object), " of length ", lo, "\n",
        sep = "")
    ## cat("Indexes: ", "\n", sep="")
    print(object@listData)
    cat("index of each column: ", "\n", sep="")
    print(object@index)   
})

.validate_LazyList <- function(x)
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
setValidity2("LazyList", .validate_LazyList)

###---------------------------------------
## Utility functions for .LazyList
###---------------------------------------

.lazyIndex_inuse <- function(LazyList)  
{
    ## browser(); browser()
    ## 1. check if any duplicate in @listData. Modify @index
    ## correspondingly.
    dupIndex <- duplicated(LazyList@listData)
    if (any(dupIndex)) {
        dupOldIndex <- which(dupIndex)
        for (i in dupOldIndex) {
            dupNewIndex <- match(LazyList@listData[i], LazyList@listData)  ## replace
            ## multiple
            ## index.
            LazyList@index[LazyList@index %in% i] <- dupNewIndex
        }
    }
    ## 2. check if all @listData in use in @index. remove if not.
    old_listData_ind <- seq_len(length(LazyList@listData)) 
    index_inuse <- old_listData_ind %in% unique(LazyList@index)

    new_listData <- LazyList@listData[index_inuse]
    new_index <- match(LazyList@index, old_listData_ind[index_inuse])

    ans <- .LazyList(new_listData, index=new_index)
    ## 3. reorder indexes (index for 1st columns in @listData[[1]])
    ans[TRUE]
}

.update_index <- function(LazyList, j, value)
{
    ## browser(); browser()
    for (i in seq_along(LazyList@listData)) {
        index <- LazyList@listData[[i]]
        if (identical(value, index)) {
            LazyList@index[j] <- i
            return(LazyList)
        }
    }
    new_index_id <- length(LazyList@listData) + 1L
    LazyList@index[j] <- new_index_id
    LazyList@listData[[new_index_id]] <- value
    .lazyIndex_inuse(LazyList)
}

.update_row <- function(LazyList, value)
{
    ## browser(); browser()
    LazyList@listData <- lapply(LazyList@listData, function(index) {
        if (is.null(index)) {
            index <- value
        } else {
            if (length(value) > length(index) )
                stop("the subscripts are out of bound")
            index <- index[value]
        }
        index
    })
    .lazyIndex_inuse(LazyList)
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
    slots = c(lazyIndex = "LazyList")
)

###---------------------------------------
## Utility functions for DelayedDataFrame
###---------------------------------------
.get_index <- function(x, j)
{
    j <- x@lazyIndex@index[[j]]
    x@lazyIndex@listData[[j]]
}

## #' @description \code{update_lazyIndex}: make sure the indexes are
## #'     consistent with columns and being used. update the
## #'     \code{indexes} slot and \code{@has_index} slot in lazyIndex.
## #' @rdname DelayedDataFrame-class
## update_lazyIndex <- function(from)
## {
##     if (identical(dim(from), c(0L, 0L))) {
##         from@lazyIndex <- .LazyList()
##         return(from)
##     }     
##     lazyIndex <- from@lazyIndex
##     delayedClass <- vapply(from, is, logical(1), "DelayedArray")

##     for (i in seq_len(length(from))) {  ## FIXME, do not distinguish delayedClass or not.
##         if (delayedClass[i]) {
##             index <- from[[i]]@index  ## FIXME. index <- from[[i]]@index[[1]] ??
##         } else {
##             index <- NULL
##         }
##         lazyIndex <- .update_index(lazyIndex, i, index)
##     }
##     from
## }

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
    names <- names(x@listData)
    ans <- lapply(
        setNames(seq_len(ncol(x)), names), function(j) x[[j]]
    )
    if (!use.names)
        names(ans) <- NULL
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
    DelayedDataFrame(callNextMethod())
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
    lazyIndex <- .LazyList()
    } else {     
        lazyIndex <- .LazyList(vector("list", 1), index=rep(1L, length(from)))
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
    .validate_LazyList(x@lazyIndex)
    ## @index must have same length of ncol(x)
    if(length(x@lazyIndex@index) != ncol(x))
        return(wmsg("'x@index' must be of same length of 'ncols(x)'"))
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
    x@lazyIndex <- .update_row(x@lazyIndex, i)
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

setMethod("[", c("LazyList", "ANY", "missing", "ANY"),
          function(x, i, j, ..., drop = TRUE)
          {
              has_index <- x@index[i]
              uhas_index <- unique(has_index)
              has_index <- match(has_index, uhas_index)
              indexes <- x@listData[uhas_index]
              .LazyList(indexes, index = has_index)
          }
          )

setReplaceMethod(
    "[", c("DelayedDataFrame", "ANY"),
    function(x, i, j, ..., value)
{
    xstub <- setNames(seq_along(x), names(x))
    if (missing(j)) {
        i <- normalizeSingleBracketSubscript(i, xstub)
        x@lazyIndex <- .update_index(x@lazyIndex, i, NULL)
    } else {
        j <- normalizeSingleBracketSubscript(j, xstub)
        x@listData[j] <- lapply(j, function(j, x) x[[j]], x)
        x@lazyIndex <- .update_index(x@lazyIndex, j, NULL)
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
            x@lazyIndex <- x@lazyIndex[j]
        }
        callNextMethod()
    })


setMethod("concatenateObjects", "LazyList",
          function(x, objects=list(), use.names = TRUE,
                   ignore.mcols = FALSE, check = TRUE) 
{
    ## browser(); browser()
    for (i in seq_len(length(objects))){
        y <- objects[[i]]
        x@listData <- c(x@listData, y@listData)
    }
    x@listData <- c(x@listData, lapply(objects, function(obj) obj@listData))
    for (j in seq_len(length(objects))) {
        y <- objects[[j]]
        for (i in seq_len(length(y@index))) {
            new_index <- y@index[i]
            new_indexes <- y@listData[[new_index]]
            ind <- match(list(new_indexes), x@listData)
            x@index <- c(x@index, ind)
        }
    }
    .lazyIndex_inuse(x)
    
})

## constructing a new DelayedDataFrame
setMethod(
    "concatenateObjects", "DelayedDataFrame",
    function(x, objects = list(), use.names = TRUE, ignore.mcols = FALSE, check = TRUE)
{
    for (i in seq_len(length(objects))) {
        y <- as(objects[[i]], "DelayedDataFrame")
        x@lazyIndex <- c(x@lazyIndex, y@lazyIndex)
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

