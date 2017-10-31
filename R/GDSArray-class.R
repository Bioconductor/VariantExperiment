### =========================================================================
### GDSArray objects
### -------------------------------------------------------------------------


setClass("GDSArraySeed",
    contains = "Array", ## ?? or "gds.class"
    slots = c(
        file="character",   # Absolute path to the gds file so the object
                            # doesn't break when the user changes the working
                            # directory (e.g. with setwd()).
        dim = "integer",
        dimnames = "list"
        ## first_val="ANY"  # First value in the dataset. Needed for GDSArraySeed???? 
    )
)

setMethod("dim", "GDSArraySeed", function(x) x@dim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subset_seed_as_array()  ???????
###

.subset_GDSArraySeed_as_array <- function(seed, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(seed))
    if (any(ans_dim == 0L)) {
        ans <- seed@first_val[0]
        dim(ans) <- ans_dim
    } else {
        ans <- h5read2(seed@file, seed@name, index)
    }
    ans
}

setMethod("subset_seed_as_array", "GDSArraySeed",
    .subset_GDSArraySeed_as_array
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GDSArraySeed constructor
###

## library(SeqArray)
## file <- seqExampleFileName("gds")
## file <- seqOpen(fn)
#' @import SeqArray

.get_gdsdata_dim <- function(file)
{
    f <- seqOpen(file)
    sample.id <- seqGetData(f, "sample.id")
    variant.id <- seqGetData(f, "variant.id")
    seqClose(f)
    dim <- c(length(sample.id), length(variant.id))
    if (!is.integer(dim)) {
        if (any(dim > .Machine$integer.max)) {
            dim_in1string <- paste0(dim, collapse=" x ")
            stop(wmsg("The dimensions of GDS dataset '", file, "' are: ",
                      dim_in1string, "\n\nThe GDSArray package only ",
                      "supports datasets with all dimensions <= 2^31-1",
                      " (this is ", .Machine$integer.max, ") at the moment."))
        }
        dim <- as.integer(dim)
    }
    dim
}

### Will fail if the dataset is empty (i.e. if at least one of its dimensions
### is 0).
## .read_gdsdata_first_val <- function(file, name, ndim)
.read_gdsdata_first_val <- function(file)
{
    f <- seqOpen(file)
    sample.id <- seqGetData(f, "sample.id")
    variant.id <- seqGetData(f, "variant.id")
    suppressMessages(
        seqSetFilter(f, sample.id=sample.id[1], variant.id=variant.id[1])
        )   ## doesn't work. 
    first_val <- paste(seqGetData(f, "genotype")[, 1, ], collapse="/")
    suppressMessages(seqResetFilter(f))
    seqClose(f)
    first_val
    ## index <- rep.int(list(1L), ndim)
    ## ans <- h5read2(file, name, index)
    ## stopifnot(length(ans) == 1L)  # sanity check
    ## ans[[1L]]  # drop any attribute
}

### Return a GDSArraySeed object with NO dimnames!
### FIXME: Investigate the possiblity to store the dimnames in the HDF5 file
### and make dimnames() on the object returned by GDSArraySeed() bring them
### back.
GDSArraySeed <- function(file, name, type=NA)
{
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the HDF5 file where the dataset is located"))
    file <- file_path_as_absolute(file)
    if (!isSingleString(name))
        stop(wmsg("'name' must be a single string specifying the name ",
                  "of the dataset in the HDF5 file"))
    if (name == "")
        stop(wmsg("'name' cannot be the empty string"))
    if (!isSingleStringOrNA(type))
        stop("'type' must be a single string or NA")
    dim <- .get_gdsdata_dim(file, name)
    if (any(dim == 0L)) {
        if (is.na(type))
            stop(wmsg("This HDF5 dataset is empty! Don't know how to ",
                      "determine the type of an empty HDF5 dataset at the ",
                      "moment. Please use the 'type' argument to help me ",
                      "(see '?GDSArray' for more information)."))
        first_val <- match.fun(type)(1)  # fake value
        if (!is.atomic(first_val))
            stop(wmsg("invalid type: ", type))
    } else {
        first_val <- .read_gdsdata_first_val(file, name, length(dim))
        detected_type <- typeof(first_val)
        if (!(is.na(type) || type == detected_type))
            warning(wmsg("The type specified via the 'type' argument (",
                         type, ") doesn't match the type of this HDF5 ",
                         "dataset (", detected_type, "). Ignoring the ",
                         "former."))
    }
    new2("GDSArraySeed", file=file,
                          name=name,
                          dim=dim,
                          first_val=first_val)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GDSArray and HDF5Matrix objects
###
### We define these classes only for cosmetic reasons i.e. to hide the
### DelayedArray and DelayedMatrix classes from the user. The user will see
### and manipulate GDSArray and HDF5Matrix objects instead of DelayedArray
### and DelayedMatrix objects.
###

setClass("GDSArray", contains="DelayedArray")

setClass("HDF5Matrix", contains=c("DelayedMatrix", "GDSArray"))

### Automatic coercion method from GDSArray to HDF5Matrix silently returns
### a broken object (unfortunately these dummy automatic coercion methods don't
### bother to validate the object they return). So we overwrite it.
setAs("GDSArray", "HDF5Matrix", function(from) new("HDF5Matrix", from))

### For internal use only.
setMethod("matrixClass", "GDSArray", function(x) "HDF5Matrix")

.validate_GDSArray <- function(x)
{
    if (!is(x@seed, "GDSArraySeed"))
        return(wmsg("'x@seed' must be a GDSArraySeed object"))
    if (!DelayedArray:::is_pristine(x))
        return(wmsg("'x' carries delayed operations"))
    TRUE
}

setValidity2("GDSArray", .validate_GDSArray)

setAs("ANY", "HDF5Matrix",
    function(from) as(as(from, "GDSArray"), "HDF5Matrix")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

setMethod("DelayedArray", "GDSArraySeed",
    function(seed) DelayedArray:::new_DelayedArray(seed, Class="GDSArray")
)

### Works directly on a GDSArraySeed object, in which case it must be called
### with a single argument.
GDSArray <- function(file, name, type=NA)
{
    if (is(file, "GDSArraySeed")) {
        if (!(missing(name) && identical(type, NA)))
            stop(wmsg("GDSArray() must be called with a single argument ",
                      "when passed a GDSArraySeed object"))
        seed <- file
    } else {
        seed <- GDSArraySeed(file, name, type=type)
    }
    DelayedArray(seed)
}

