#' GDSArraySeed
#' Generate the seed for gds data for the GDSArray
#' @import SNPRelate
#' @import gdsfmt
#' @importFrom tools file_path_as_absolute
#' @import S4Vectors
#' 

## library(gdsfmt)
## library(SNPRelate)
## library(S4Vectors)
## library(tools)

### =========================================================================
### GDSArray objects
### -------------------------------------------------------------------------

#' @importClassesFrom DelayedArray Array
#' @importFrom DelayedArray subset_seed_as_array
#' @export 
setClass("GDSArraySeed",
         contains = "Array", ## from DelayedArray: A virtual class with no slots
                             ## to be extended by concrete subclasses with
                             ## an array-like semantic.
    slots = c(
        file="character",   # Absolute path to the gds file so the object
                            # doesn't break when the user changes the working
                            # directory (e.g. with setwd()).
        name="character",   # Name of the dataset in the gds file.
        dim = "integer",
        dimnames = "list",
        permute = "logical",
        first_val = "ANY"
    )
)

###
## accessors for GDSArraySeed object
###

setGeneric("gdsfile", function(x) standardGeneric("gdsfile"))
setMethod("gdsfile", "GDSArraySeed", function(x) x@file)
setMethod("dim", "GDSArraySeed", function(x) x@dim)
setMethod("dimnames", "GDSArraySeed", function(x) x@dimnames)

###
## show method for GDSArraySeed object
###
setMethod(
    "show", "GDSArraySeed",
    function(object){
        cat("GDSArraySeed\n",
            "gds file: ", object@file, "\n",
            "array data: ", object@name, "\n",
            "dim: ", nrow(object), " x ", ncol(object), "\n",
            sep="")
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subset_seed_as_array()
###

.subset_GDSArraySeed_as_array <- function(seed, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(seed))
    if (any(ans_dim == 0L)) {
        ans <- seed@first_val[0]
        dim(ans) <- ans_dim
    } else {
        f <- openfn.gds(seed@file)
        on.exit(closefn.gds(f))
        if(seed@permute){
            index <- index[rev(seq_len(length(index)))] ## always 2 dimensional? SeqArray...
            ans <- t(readex.gdsn(index.gdsn(f, seed@name), index))
        }else{
            ans <- readex.gdsn(index.gdsn(f, seed@name), index)
        }
    }
    ans
}

setMethod("subset_seed_as_array", "GDSArraySeed",
    .subset_GDSArraySeed_as_array
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GDSArraySeed constructor
###

## file <- system.file("extdata", "hapmap_geno.gds", package = "SNPRelate")
## snpgdsSummary(file)
## f <- openfn.gds(file)

.get_gdsdata_arraynodes <- function(file){
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    names.gdsn <- ls.gdsn(f)
    a <- lapply(names.gdsn, function(x) ls.gdsn(index.gdsn(f, x)))
    a[lengths(a)==0] <- ""
    n <- rep(names.gdsn, lengths(a))
    all.gdsn <- paste(n, unlist(a), sep="/")
    all.gdsn <- sub("/$", "", all.gdsn)

    isarray <- sapply(all.gdsn, function(x)objdesp.gdsn(index.gdsn(f, x))$is.array)
    dims <- lapply(all.gdsn, function(x)objdesp.gdsn(index.gdsn(f, x))$dim)
    names(dims) <- all.gdsn
    names(dims)[unlist(isarray) & lengths(dims) > 1]
}

.read_gdsdata_sampleInCol <- function(file, node){
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    rd <- names(get.attr.gdsn(index.gdsn(f, node)))
    if ("snp.order" %in% rd) sampleInCol <- TRUE   ## snpfirstdim (in row)
    if ("sample.order" %in% rd) sampleInCol <- FALSE
    sampleInCol
}

.get_gdsdata_dim <- function(file, node){
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    dim <- objdesp.gdsn(index.gdsn(f, node))$dim
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

.get_gdsdata_dimnames <- function(file, node){
    summ <- SNPRelate::snpgdsSummary(file, show=FALSE)
    sample.id <- summ$sample.id
    snp.id <- as.character(summ$snp.id)
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    rd <- names(get.attr.gdsn(index.gdsn(f, node))) ## ?
    if ("snp.order" %in% rd){
        dimnames <- list(snp.id = snp.id, sample.id = sample.id)
    }else{
        dimnames <- list(sample.id = sample.id, snp.id = snp.id)
    }
    dimnames
    ## dimnames <- SNPRelate::snpgdsSummary(file, show=FALSE)
    if (!is.list(dimnames)) {
        stop(wmsg("The dimnames of GDS dataset '", file, "' should be a list!"))
    }
    dimnames
}
## if dimnames is shorter than the corresponding dimensions, use NULL to extend.??

.read_gdsdata_first_val <- function(file, node){
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    first_val <- read.gdsn(index.gdsn(f, node), start=c(1,1), count=c(1,1))
    first_val
}


## vcf.fn <- system.file("extdata", "sequence.vcf", package="SNPRelate")
## snpgdsVCF2GDS(vcf.fn, "test2.gds", method="biallelic.only", snpfirstdim=TRUE)
## SNPRelate::snpgdsSummary("test2.gds")  ## 2snp X 3samples
## f2 <- openfn.gds("test2.gds")
## closefn.gds(f2)

#' @param gfile the gds file
#' @param node the gds nodes to be read into GDSArray
#' @export
#' 
###
## GDSArraySeed constructor
###
GDSArraySeed <- function(file, name){
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the gds file where the dataset is located."))
    file <- file_path_as_absolute(file)
    if (!isSingleStringOrNA(name))
        stop("'type' must be a single string or NA")
    ## get the nodes that are array data. 
    arrayNodes <- .get_gdsdata_arraynodes(file)
    if(!name %in% arrayNodes){
        stop(wmsg("the `name` must be a rectangular array data."))
    }
    ## force to read in data with variant * sample order.
    ## file <- .write_gdsdata_sampleInCol(file, node = name)
    dim <- .get_gdsdata_dim(file, node = name)
    dimnames <- .get_gdsdata_dimnames(file, node = name)
    first_val <- .read_gdsdata_first_val(file, node = name)

    permute = !.read_gdsdata_sampleInCol(file, node = name)
    if(permute){
        dim <- rev(dim)
        dimnames <- dimnames[rev(seq_len(length(dimnames)))]
    }
    new2("GDSArraySeed", file=file,
         name=name,
         dim=dim,
         dimnames = dimnames,
         permute = permute,
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

#' @importClassesFrom DelayedArray DelayedArray DelayedMatrix
setClass("GDSArray", contains="DelayedArray")
setClass("GDSMatrix", contains=c("DelayedMatrix", "GDSArray"))

### Automatic coercion method from GDSArray to GDSMatrix
setAs("GDSArray", "GDSMatrix", function(from) new("GDSMatrix", from))

### For internal use only.
## setMethod("matrixClass", "GDSArray", function(x) "GDSMatrix")

.validate_GDSArray <- function(x)
{
    if (!is(x@seed, "GDSArraySeed"))
        return(wmsg("'x@seed' must be a GDSArraySeed object"))
    if (!DelayedArray:::is_pristine(x))
        return(wmsg("'x' carries delayed operations"))
    TRUE
}

#' @importFrom S4Vectors setValidity2
setValidity2("GDSArray", .validate_GDSArray)

setAs("ANY", "GDSMatrix",
      function(from) as(as(from, "GDSArray"), "GDSMatrix")
    )


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @importFrom DelayedArray DelayedArray
setMethod("DelayedArray", "GDSArraySeed",
              function(seed) DelayedArray:::new_DelayedArray(seed, Class="GDSArray")
          )
setMethod("gdsfile", "DelayedArray", function(x) gdsfile(seed(x)))

#' @export
### Works directly on a GDSArraySeed object, in which case it must be called
### with a single argument.
GDSArray <- function(file, name=NA){
    if (is.na(name)){
        name <- "genotype"
    } 
    if (is(file, "GDSArraySeed")) {
        seed <- file
    } else {
        seed <- GDSArraySeed(file, name)
    }
    as(DelayedArray(seed), "GDSMatrix")
}

setMethod("gdsfile", "GDSArray", function(x) gdsfile(seed(x)))
