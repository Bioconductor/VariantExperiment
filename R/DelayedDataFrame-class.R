###----------  
## lazyList
###----------
## .lazyList <- setClass(
##     "lazyList",
##     slots = c(
##         indexes = "list",
##         has_index = "integer"
##     )
## )

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
    indexLength <- vapply(indexes, function (index) {
        length(index[[1]])  ## FIXME: do not use list for indexes. only keep i subscripts.
    }, integer(1))
    
    if (length(unique(indexLength)) > 1 && all(unique(indexLength) > 0))
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
    orig <- seq_len(length(LazyList@listData)) 
    index_inuse <- orig %in% unique(LazyList@index)
    old <- orig[index_inuse]
    LazyList@listData <- LazyList@listData[index_inuse]
    new <- seq_len(length(LazyList@listData))

    LazyList@index <- new[match(LazyList@index, old)]
    LazyList
}

.update_index <- function(LazyList, j, value)
{
    for (i in seq_along(LazyList@listData)) {
        index <- LazyList@listData[[i]]
        if (identical(value, index)) {
            LazyList@index[j] <- i
            return(LazyList)
        }
    }
    new_index_id <- length(LazyList@listData) + 1L
    LazyList@index[j] <- new_index_id
    index <- value
    LazyList@listData[[new_index_id]] <- index
    .lazyIndex_inuse(LazyList)
}

.update_row <- function(LazyList, value)
{
    LazyList@listData <- lapply(LazyList@listData, function(index) {
        if (is(index, "list")) {
            if (is.null(index[[1]])) {
                index[[1]] <- value
            } else {
            if (length(value) > length(index[[1]]) )
                stop("the subscripts are out of bound")
            index[[1]] <- index[[1]][value]
            }
        } else {
            if (!is.null(index))
                index <- index[value]
        }
        index
    })
    LazyList
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

#' @description \code{update_lazyIndex}: make sure the indexes are
#'     consistent with columns and being used. update the
#'     \code{indexes} slot and \code{@has_index} slot in lazyIndex.
#' @rdname DelayedDataFrame-class
update_lazyIndex <- function(from)
{
    if (identical(dim(from), c(0L, 0L))) {
        from@lazyIndex <- .LazyList()
        return(from)
    }     
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
    ## ## assign new indexes into DelayedArray columns.
    ## for (i in seq_len(length(from))) {
    ##     if(delayedClass[i])
    ##         from[[i]]@index <- .get_index(from, i)
    ## }
    from
}

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
## methods
###-------------

setMethod("getListElement", "DelayedDataFrame", function(x, i)
{
    index <- .get_index(x, i)[[1]]  ## FIXME after modifying the indexes to be not list.
    extractROWS(x@listData[[i]], index)
})

## setMethod("as.list", "DelayedDataFrame", function(x, use.names=TRUE)
## {
##     ans <- as.list(DataFrame(x))
##     if (!use.names)
##         names(ans) <- NULL
## })

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
    ## initialize lazyIndex
    ## if (identical(dim(from), c(0L, 0L))) {
    ## browser(); browser()
    lazyIndex <- .LazyList()
    ## } else {     
    delayedClass <- vapply(from, is, logical(1), "DelayedArray")

    for (i in seq_len(length(from))) {
        if (delayedClass[i]) {
            index <- vector("list", length(dim(from[[i]])))
            ## index <- from[[i]]@index
        } else {
            index <- vector("list", 1)  ## ordinary vector columns also has list index.
        }
        lazyIndex <- .update_index(lazyIndex, i, index)
    }
    ## ## remove any index if not in use
    ## from@lazyIndex <- .lazyIndex_inuse(lazyIndex)
    ## ## reorder indexes
    ## from@lazyIndex <- from@lazyIndex[TRUE]
    ## lazyIndex <- .LazyList(
    ##     listData = vector("list", 1),
    ##     index = rep(1L, length(from)))
    ans <- .DelayedDataFrame(from, lazyIndex = lazyIndex)
    ans
    ## update_lazyIndex(ans)
})

setAs("DelayedDataFrame", "DataFrame", function(from){
    ## realize the lazyIndex to each column.
    ## browser();browser()
    indexes <- from@lazyIndex@listData
    rowSubscripts <- lapply(indexes, function(index) index[[1]])
    nrowSubscript <- nrow(from)
    isNullRowSubscripts <- vapply(rowSubscripts, is.null, logical(1))
    if(all( isNullRowSubscripts ))
        return(S4Vectors:::new_DataFrame(listData=from@listData))  ## is it expensive??
    listData <- lapply(
        setNames(seq_len(ncol(from)), names(from)), function(j)
        {
            rowSubscript <- .get_index(from, j)[[1]]
            if(is.null(rowSubscript))
                rowSubscript <- TRUE   
            extractROWS(from[[j]], rowSubscript)
        })
    ans <- S4Vectors:::new_DataFrame(listData=listData, nrows=nrowSubscript)
    if (!is.null(rownames(from))) {
        rnmSubscript <- rowSubscripts[[ which(!isNullRowSubscripts)[1] ]]
        rownames <- make.unique(extractROWS(rownames(from), rnmSubscript))
        ans@rownames <- rownames
    }
    ans
})

###
setAs("ANY", "DelayedDataFrame", function(from){
    df <- as(from, "DataFrame")
    as(df, "DelayedDataFrame")
})

###-----------------
## show method
###----------------

setMethod("show", "DelayedDataFrame", function (object) 
{
    ## browser(); browser()
    nhead <- get_showHeadLines()
    ntail <- get_showTailLines()
    ## rowSubscripts <- object@lazyIndex@listData[[1]][[1]]  ## FIXME?
    ## nr <- ifelse(is.null(rowSubscripts), nrow(object), length(rowSubscripts))
    nr <- nrow(object)
    nc <- ncol(object)
    cat(class(object), " with ", nr, ifelse(nr == 1, " row and ", 
        " rows and "), nc, ifelse(nc == 1, " column\n", " columns\n"), 
        sep = "")
    if (nr > 0 && nc > 0) {
        nms <- rownames(object)
        if (nr <= (nhead + ntail + 1L)) {
            object <- DataFrame(object)
            out <- as.matrix(format(as.data.frame(lapply(object, 
                showAsCell), optional = TRUE)))
            if (!is.null(nms)) 
                rownames(out) <- nms
        }
        else {
            out <- rbind(
                as.matrix(format(as.data.frame(lapply(
                    setNames(seq_len(length(object)), names(object)), 
                    function(j) {
                        index <- .get_index(object, j)[[1]]
                        ## FIXME: if null, should get the length of
                        ## the other index, not nrow(object)
                        if (is.null(index))
                            index <- seq_len(nrow(object))
                        headSubscript <- head(index, nhead)
                        showAsCell(extractROWS(object[[j]], headSubscript))
                    })),
                    optional = TRUE)),
                rbind(rep.int("...", nc)),
                as.matrix(format(as.data.frame(lapply(
                    setNames(seq_len(length(object)), names(object)), 
                    function(j) {
                        index <- .get_index(object, j)[[1]]
                        if (is.null(index))
                            index <- seq_len(nrow(object))
                        tailSubscript <- tail(index, ntail)
                        showAsCell(extractROWS(object[[j]], tailSubscript))
                    })),
                    optional = TRUE)))
            rownames(out) <- S4Vectors:::.rownames(nms, nr, nhead, ntail)  ## FIXME, rownames
        }
        classinfo <- matrix(unlist(lapply(object, function(x) {
            paste0("<", classNameForDisplay(x)[1], ">")
        }), use.names = FALSE), nrow = 1, dimnames = list("", 
            colnames(out)))
        out <- rbind(classinfo, out)
        print(out, quote = FALSE, right = TRUE)
    }
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
    ## x <- update_lazyIndex(x)  
    i <- normalizeSingleBracketSubscript(
        i, x, exact = FALSE, allow.NAs = TRUE, as.NSBS = FALSE)
    x@lazyIndex <- .update_row(x@lazyIndex, i)
    ## slot(x, "listData") <- lapply(
    ##     setNames(seq_len(ncol(x)), names(x)), function(j)
    ##     {
    ##         if(!is(x[[j]], "DelayedArray"))  ## non-delayed_ops object.
    ##             extractROWS(x[[j]], i)
    ##         else {                    ## DelayedArray objects.
    ##             a <- x[[j]]
    ##             a@index <- .get_index(x, j)
    ##             a
    ##         }
    ##     })
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
    ## x@lazyIndex <- do.call("c", c(list(x), objects))
    ## browser()
    x@lazyIndex <- c(x@lazyIndex, y@lazyIndex)
    callNextMethod()
})

setMethod("concatenateObjects", "LazyList",
          function(x, objects=list(), use.names = TRUE,
                   ignore.mcols = FALSE, check = TRUE) 
{
    if (!isTRUEorFALSE(use.names)) 
        stop("'use.names' must be TRUE or FALSE")
    for (j in seq_len(length(objects))) {
        lazyIndex <- objects[[j]]
        for (i in seq_len(length(lazyIndex@index))) {
            new_index <- lazyIndex@index[i]
            new_indexes <- lazyIndex@listData[[new_index]]
            x@index <- c(x@index, new_index)
            x <- .update_index(x, length(x@index), value=new_indexes)
        }
    }
    x
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

