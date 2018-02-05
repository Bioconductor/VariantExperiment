## ## set class 
## setClass(
##     "GDSDataFrameSeed",
##     contains = "Array",
##     slots = c(
##         file = "character",
##         names = "character",
##         dim = "integer",
##         dimnames = "list",
##         first_val = "ANY"
##     )
## 0)

## ## methods
## setMethod("dim", "GDSDataFrameSeed", function(x) x@dim)
## setMethod("dimnames", "GDSDataFrameSeed", function(x) x@dimnames)
## setMethod("gdsfile", "GDSDataFrameSeed", function(x) x@file)

## ## show method
## setMethod("show", "GDSDataFrameSeed",
##           function(object){
##               cat("GDSDataFrameSeed\n",
##                   "gds file: ", object@file, "\n",
##                   "annotation data: ", paste(object@names, collapse=", "), "\n",
##                   ## "dim: ", nrow(object), " x ", ncol(object), "\n",
##                   "dim: ", paste(dim(object), collapse=" x "), "\n",
##                   "first value: ", object@first_val, "\n",
##                   sep="")
##           }
##           )

## ## seed constructor
## GDSDataFrameSeed <- function(file, names){
##     if (!isSingleString(file))
##         stop(wmsg("'file' must be a single string specifying the path to ",
##                   "the gds file where the dataset is located."))
##     if (!is.character(names))
##         stop("'type' must be a single string or NA")
##     file <- file_path_as_absolute(file)
    
##     ff <- .get_gdsdata_fileFormat(file)
##     if(ff == "SNP_ARRAY"){
##         f <- snpgdsOpen(file)
##         on.exit(snpgdsClose(f))
##     }else if(ff == "SEQ_ARRAY"){
##         f <- seqOpen(file)
##         on.exit(seqClose(f))
##     }else{
##         f <- openfn.gds(file)
##         on.exit(closefn.gds(f))
##     }
##     dim <- c(.get_gdsdata_dim(f, names[1]), length(names))
##     dimnames <- .get_gdsdata_dimnames(f, names[1], "SEQ_ARRAY")
##     dimnames$annotation.id <- names
##     first_val <- .read_gdsdata_first_val(f, names[1])

##     new("GDSDataFrameSeed",
##         file = file,
##         names = names,
##         dim = dim,
##         dimnames = dimnames,
##         first_val = first_val
##         )
## }

## .extract_DF_from_GDSDataFrameSeed <- function(x, index)
## {
##     ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
##     if (any(ans_dim == 0L)){
##         ans <- x@first_val[0]
##         dim(ans) <- ans_dim
##     } else {
##         f <- openfn.gds(x@file)
##         on.exit(closefn.gds(f))
##         ## take the cidx for gdsnodes to apply the subsetting.
##         ## cidx <- index[[2]]
##         ## names <- x@names[cidx]
##         ## only take the row index, and apply to all names.
##         ridx <- index[[1]]  
##         ans <- lapply(x@names, function(node) readex.gdsn(index.gdsn(f, node), ridx))
##         ans <- setNames(DataFrame(ans), x@names)
##     }
##     ans
## }
## setMethod("extract_array", "GDSDataFrameSeed", .extract_DF_from_GDSDataFrameSeed)


## #######################################################
## ### example
## #######################################################

## file <- SeqArray::seqExampleFileName("gds")
## f <- seqOpen(file)
## allnodes <- .get_gdsdata_allNodes(f)
## alldims <- lapply(allnodes, function(x) .get_gdsdata_dim(f, x))
## var.id <- vapply(alldims, function(x) length(x) == 1 & all(x == 1348), logical(1))
## var.node <- allnodes[var.id]
## names <- c("variant.id", "annotation/id", "annotation/info/AA")

## seed <- GDSDataFrameSeed(file, names) 
## extract_array(seed, list(1:5, 1:3))
## DelayedArray(seed)

## seed1 <- GDSDataFrameSeed(file, c("variant.id", "position", "chromosome"))
## extract_array(seed1, list(1:5, 1:2))
## DelayedArray(seed1)
## ## <1348 x 3> DelayedMatrix object of type "integer":
## ## Error in prettyNum(.Internal(format(x, trim, digits, nsmall, width, 3L,  : 
## ##   first argument must be atomic

## seed2 <- GDSDataFrameSeed(file, c("variant.id", "position"))
## extract_array(seed2, list(1:5, 1:3))
## DelayedArray(seed1)

#######################################
## Martin's example
#######################################

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
