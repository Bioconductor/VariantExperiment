#' @import SNPRelate
#' @import gdsfmt
#' @importFrom tools file_path_as_absolute
#' @importFrom S4Vectors isSingleString isSingleStringOrNA
#' 

## library(gdsfmt)
## library(SNPRelate)
## library(S4Vectors)
## library(tools)

### =========================================================================
### GDSArray objects
### -------------------------------------------------------------------------

setClass("GDSArraySeed",
    contains = "array", ## ?? or "gds.class"
    slots = c(
        file="character",   # Absolute path to the gds file so the object
                            # doesn't break when the user changes the working
                            # directory (e.g. with setwd()).
        dim = "integer",
        dimnames = "list",
        first_val="ANY"  # First value in the dataset. Needed for GDSArraySeed???? 
    )
)

###
## accessors for GDSArraySeed object
###
setMethod("dim", "GDSArraySeed", function(x) x@dim)
setMethod("dimnames", "GDSArraySeed", function(x) x@dimnames)

###
## show method for GDSArraySeed object
###
setMethod(
    "show", "GDSArraySeed",
    function(object){
        cat("GDSArraySeed\n")
        cat("gds file:", object@file, "\n")
        cat("The total number of samples: ", unname(dim(object)["sample.id"]), "\n")
        cat("The total number of snps: ", unname(dim(object)["snp.id"]), "\n")
        cat("The first value of assay:", object@first_val, "\n")
    }
)


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

## ## library(SeqArray)
## ## file <- seqExampleFileName("gds")
## ## library(SNPRelate)
## ## file <- seqOpen(fn)

## file <- system.file("extdata", "hapmap_geno.gds", package = "SNPRelate")
## snpgdsSummary(file)
## f <- snpgdsOpen(file)

.get_gdsdata_dim <- function(file){
    dim <- sapply(snpgdsSummary(file, show=FALSE), length)
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

.get_gdsdata_dimnames <- function(file){
    dimnames <- snpgdsSummary(file, show=FALSE)
    if (!is.list(dimnames)) {
        stop(wmsg("The dimnames of GDS dataset '", file, "' should be a list!"))
    }
    dimnames
}
## if dimnames is shorter than the corresponding dimensions, use NULL to extend.??

.read_gdsdata_first_val <- function(file){
    f <- snpgdsOpen(file)
    first_val <- read.gdsn(index.gdsn(f, "genotype"), start=c(1,1), count=c(1,1))
    snpgdsClose(f)
    first_val
}

## vcf.fn <- system.file("extdata", "sequence.vcf", package="SNPRelate")
## snpgdsVCF2GDS(vcf.fn, "test2.gds", method="biallelic.only", snpfirstdim=TRUE)
## snpgdsSummary("test2.gds")  ## 2snp X 3samples
## f2 <- snpgdsOpen("test2.gds")
## snpgdsClose(f2)

### Return a GDSArraySeed object with NO dimnames!
### FIXME: Investigate the possiblity to store the dimnames in the GDS file
### and make dimnames() on the object returned by GDSArraySeed() bring them
### back.
GDSArraySeed <- function(file, type=NA){
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the gds file where the dataset is located."))
    file <- file_path_as_absolute(file)
    if (!isSingleStringOrNA(type))
        stop("'type' must be a single string or NA")
    dim <- .get_gdsdata_dim(file)
    dimnames <- .get_gdsdata_dimnames(file)
    first_val <- .read_gdsdata_first_val(file)
    ## if (any(dim == 0L)) {
    ## if (is.na(type))
    ##     stop(wmsg("This gds dataset is empty! Don't know how to ",
    ##               "determine the type of an empty gds dataset at the ",
    ##               "moment. Please use the 'type' argument to help me ",
    ##               "(see '?GDSArray' for more information)."))
    
    ## first_val <- match.fun(type)(1)  # fake value
    ## if (!is.atomic(first_val))
    ##         stop(wmsg("invalid type: ", type))
    ## } else {
    ## first_val <- .read_gdsdata_first_val(file, length(dim))
    ## detected_type <- typeof(first_val)
    ## if (!(is.na(type) || type == detected_type))
    ##     warning(wmsg("The type specified via the 'type' argument (",
    ##                  type, ") doesn't match the type of this GDS ",
    ##                  "dataset (", detected_type, "). Ignoring the ",
    ##                  "former."))
    ## }
    new2("GDSArraySeed", file=file,
         dim=dim,
         dimnames = dimnames,
         first_val=first_val)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GDSArray and GDSMatrix objects
###
### We define these classes only for cosmetic reasons i.e. to hide the
### DelayedArray and DelayedMatrix classes from the user. The user will see
### and manipulate GDSArray and GDSMatrix objects instead of DelayedArray
### and DelayedMatrix objects.
###

setClass("GDSArray", contains="DelayedArray")

setClass("GDSMatrix", contains=c("DelayedMatrix", "GDSArray"))

### Automatic coercion method from GDSArray to GDSMatrix silently returns
### a broken object (unfortunately these dummy automatic coercion methods don't
### bother to validate the object they return). So we overwrite it.
setAs("GDSArray", "GDSMatrix", function(from) new("GDSMatrix", from))

### For internal use only.
setMethod("matrixClass", "GDSArray", function(x) "GDSMatrix")

.validate_GDSArray <- function(x)
{
    if (!is(x@seed, "GDSArraySeed"))
        return(wmsg("'x@seed' must be a GDSArraySeed object"))
    if (!DelayedArray:::is_pristine(x))
        return(wmsg("'x' carries delayed operations"))
    TRUE
}

setValidity2("GDSArray", .validate_GDSArray)

setAs("ANY", "GDSMatrix",
    function(from) as(as(from, "GDSArray"), "GDSMatrix")
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

