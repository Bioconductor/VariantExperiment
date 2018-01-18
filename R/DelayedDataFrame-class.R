## #' GDSArraySeed
## #' Generate the seed for gds data for the GDSArray
## #' @import S4Vectors
## #' @import methods

### =========================================================================
### GDSArray objects
### -------------------------------------------------------------------------
## #' @importClassesFrom S4Vectors DataFrame
## ## #' @export 
## setClass("DelayedDataFrame",
##          contains = "DataFrame", ## from S4VectorsA data.frame-like interface for
##                                  ## S4 objects that implement length() and `[`
##          slots = c(
##              file="character",      ## Absolute path to the gds file so the object
##                                ## doesn't break when the user changes the working
##                                ## directory (e.g. with setwd()).
##              name="character", ## character vector, gds node names in the gds file.
##              dim = "integer",
##              dimnames = "list",  
##              permute = "logical",   ## logical vector.
##              first_val = "ANY"
##          )
##          )
## )

## library(DelayedArray)

## .DataFrameSeed <- setClass(
##     "DataFrameSeed",
##     contains = "Array",
##     slots = c(
##         name = "character",
##         dim = "integer",
##         dimnames = "list"
##     )
## )

## setMethod(
##     "extract_array", "DataFrameSeed",
##     function(x, index)
##     {
##         data(list=x@name)
##         object <- get(x@name)
##         rowidx <- index[[1]]
##         colidx <- index[[2]]
##         answer <- lapply(object[colidx], `[`, rowidx)
##         ## as(answer, "DataFrame")  ## DataFrame / data.frame both work.
##                                     ## only DelayedDataFrame(seed) does not
##                                     ## show properly. 
##         as.data.frame(answer)
##         ## array(as.data.frame(answer))
##         ## answer
##     })

## obj <- new("DataFrameSeed", name="mtcars", dim=dim(mtcars), dimnames=dimnames(mtcars))
## ## obj <- .DataFrameSeed(name="mtcars", dim=dim(mtcars), dimnames=dimnames(mtcars))
## a <- extract_array(obj, list(1:5, 2:5))
## DelayedDataFrame(obj)

## obj <- new("DataFrameSeed", name="iris", dim=dim(iris), dimnames=dimnames(iris))
## extract_array(obj, list(1:5, 2:4))


## ## .DelayedDataFrame <- setClass("DelayedDataFrame", contains = "DelayedArray")
## ## setClass("GDSArray", contains="DelayedArray")
## ## setClass("GDSMatrix", contains=c("DelayedMatrix", "GDSArray"))
## setClass("DelayedDataFrame", contains = "DelayedArray")
## setGeneric("DelayedDataFrame", function(x) standardGeneric("DelayedDataFrame"))
## setMethod("DelayedDataFrame", "DataFrameSeed",
##           function(x)
##               DelayedArray:::new_DelayedArray(x, Class="DelayedDataFrame")
##           )
## ddf <- DelayedDataFrame(obj)

## ddf = DelayedDataFrame(DelayedArray(obj))
## as(obj, "DelayedDataFrame")

## str(DelayedArray(obj))
## str(.DelayedDataFrame(DelayedArray(obj)))

## obj1 <- .DataFrameSeed(name="iris", dim=dim(iris), dimnames=dimnames(iris))
## extract_array(obj1, list(1:5, 2:5))
## ddf1 <- .DelayedDataFrame(DelayedArray(obj1))
## ddf1 <- .DelayedDataFrame(obj1)
## ddf1 <- DelayedArray(obj1)
