#' @description \code{dim}, \code{dimnames}: dimension and dimnames of
#'     object contained in the GDS file.
#' @rdname GDSArray
#' @exportMethod dim
setMethod("dim", "GDSArraySeed", function(x) x@dim)

#' @rdname GDSArray
#' @exportMethod dimnames
setMethod("dimnames", "GDSArraySeed", function(x) x@dimnames)

#' @description \code{gdsfile}: on-disk location of GDS file
#'     represented by this object.
#' @param x GDSArray, GDSMatrix, GDSArraySeed or SummarizedExperiment object.
#' @rdname GDSArray
setGeneric("gdsfile", function(x) standardGeneric("gdsfile"))
#' @rdname GDSArray
#' @exportMethod gdsfile
setMethod("gdsfile", "GDSArraySeed", function(x) x@file)

#' @rdname GDSArray
#' @importFrom DelayedArray seed
## #' @exportMethod gdsfile
setMethod("gdsfile", "GDSArray", function(x) gdsfile(seed(x)))

#' @rdname GDSArray
setMethod("gdsfile", "DelayedArray", function(x) gdsfile(seed(x)))

#' @rdname GDSArray
#' @import SeqArray
setMethod("gdsfile", "SummarizedExperiment", function(x) {
    vapply(assays(x), gdsfile, character(1))
})

#' @rdname GDSArray
#' @exportMethod "gdsfile<-"
setGeneric("gdsfile<-",
           function(x, value) standardGeneric("gdsfile<-"),
                      signature="x")
setReplaceMethod("gdsfile", "GDSArraySeed",
                 function(x, value){
                     new_filepath <- tools::file_path_as_absolute(value)
                     ## ## Check dim compatibility.
                     ## gfile_new <- openfn.gds(value)  ## error.
                     ## on.exit(closefn.gds(gfile_new))
                     ## new_dim <- .get_gdsdata_dim(gfile_new, x@name)
                     ## if(!x@permute) object_dim <- dim(x)
                     ## else object_dim <- rev(dim(x))
                     ## if (!identical(new_dim, object_dim)) {
                     ##     new_dim_in1string <- paste0(new_dim, collapse=" x ")
                     ##     dim_in1string <- paste0(object_dim, collapse=" x ")
                     ##     stop(wmsg("dimensions (", new_dim_in1string, ") ",
                     ##               "of GDS dataset '", x@name, "' ",
                     ##               "from file '", value, "' are not ",
                     ##               "as expected (", dim_in1string, ")"))
                     ## }
                     
                     ## ## Check first val compatibility.
                     ## new_first_val <- .read_gdsdata_first_val(gfile_new, x@name)
                     ## gfile_object <- openfn.gds(x@file)
                     ## on.exit(closefn.gds(gfile_object))
                     ## first_val <- .read_gdsdata_first_val(gfile_object, x@name)
                     ## if (!identical(new_first_val, first_val))
                     ##     stop(wmsg("first value in GDS dataset '", x@name, "' ",
                     ##               "from file '", value, "' is not ",
                     ##               "as expected"))
                     ## Set new path.
                     BiocGenerics:::replaceSlots(x, file=value, check=FALSE)
                     ## x@file <- new_filepath
                     ## x
                 }
                 )

setReplaceMethod("gdsfile", "GDSArray",
                 function(x, value){
                     new_filepath <- tools::file_path_as_absolute(value)
                     BiocGenerics:::replaceSlots(seed(x), file=value, check=FALSE)
                 })
