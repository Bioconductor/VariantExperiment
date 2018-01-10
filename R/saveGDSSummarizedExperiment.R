### =========================================================================
### Save/load a gdsfmt-based SummarizedExperiment object
### -------------------------------------------------------------------------

.create_dir <- function(dir, replace)
{
    if (dir.exists(dir)) {
        if (!replace)
            stop(wmsg("Directory \"", dir, "\" already exists. ",
                      "Use 'replace=TRUE' to replace it. ",
                      "Its content will be lost!"))
        if (unlink(dir, recursive=TRUE) != 0L)
            stop("failed to delete directory \"", dir, "\"")
    } else if (file.exists(dir)) {
        stop(wmsg("\"", dir, "\" already exists and is a file, ",
                  "not a directory"))
    }
    if (!suppressWarnings(dir.create(dir)))
        stop("cannot create directory \"", dir, "\"")
}

.create_seqgds <- function(se, gds_path){
    gfile <- createfn.gds(filename=gds_path)
    on.exit(closefn.gds(gfile))
    put.attr.gdsn(gfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")
    dscp <- addfolder.gdsn(gfile, "description", replace=TRUE)
    put.attr.gdsn(dscp, "source.format", "GDSArray-based SummarizedExperiment")
    add.gdsn(gfile, "sample.id", colnames(se), compress="LZMA_RA", closezip=TRUE)
    add.gdsn(gfile, "variant.id", as.integer(rownames(se)),
             compress="LZMA_RA", closezip=TRUE)
    add.gdsn(gfile, "chromosome", as.character(seqnames(rowRanges(se))),
             compress="LZMA_RA", closezip=TRUE)
    add.gdsn(gfile, "position", ranges(rowRanges(se))@start,
             compress="LZMA_RA", closezip=TRUE) ## ??
    if(all(c("REF", "ALT") %in% names(rowData(se)))){
        ref <- mcols(rowRanges(se))$REF
        alt <- mcols(rowRanges(se))$ALT
        if(is(ref, "List"))  ## "DNAStringSet" in-memory way
            ref <- unlist(lapply(ref, function(x) paste(x, collapse=",")))
        if(is(alt, "List"))  ## "DNAStringSetList" in-memory way
            alt <- unlist(lapply(alt, function(x) paste(x, collapse=",")))
        add.gdsn(gfile, "allele", val = paste(ref, alt, sep=","),
                 compress="LZMA_RA", closezip=TRUE)
    }
}

.create_snpgds <- function(se, gds_path){
    ## use "LZMA_RA" as default "compress" argument
    gfile <- createfn.gds(filename=gds_path)
    on.exit(closefn.gds(gfile))
    put.attr.gdsn(gfile$root, "FileFormat", "SNP_ARRAY")
    dscp <- addfolder.gdsn(gfile, "description", replace=TRUE)
    put.attr.gdsn(dscp, "source.format", "GDSArray-based SummarizedExperiment")
    ## "description" not required for "SNP_ARRAY" gds file.
    add.gdsn(gfile, "sample.id", colnames(se), compress="LZMA_RA", closezip=TRUE)
    add.gdsn(gfile, "snp.id", as.integer(rownames(se)), compress="LZMA_RA",
             closezip=TRUE, replace=T)
    ## if("ID" %in% names(rowData(se)))
    ##     add.gdsn(gfile, "snp.rs.id", rowData(se)$ID,
    ##              compress="LZMA_RA", closezip=TRUE)
    ## add.gdsn(gfile, "snp.rs.id", as.integer(rownames(se)), compress="LZMA_RA", replace=T)
    add.gdsn(gfile, "snp.chromosome", as.character(seqnames(SummarizedExperiment::rowRanges(se))),
             compress="LZMA_RA", closezip=TRUE)
    add.gdsn(gfile, "snp.position", ranges(rowRanges(se))@start,
             compress="LZMA_RA", closezip=TRUE) ## ??
    if(all(c("ALLELE1", "ALLELE2") %in% names(rowData(se))))
        add.gdsn(gfile, "snp.allele",
                 val = paste(mcols(rowRanges(se))$ALLELE1,
                             mcols(rowRanges(se))$ALLELE2,
                             sep="/"),
                 compress="LZMA_RA", closezip=TRUE
                 )  ## "," separated.
}
 
.write_sedata_as_gds <- function(data, name, ff, gds_path, chunk_size = 1000, nrow, ncol) {
    gfile <- openfn.gds(gds_path, readonly=FALSE)
    on.exit(closefn.gds(gfile))

    if(is(data, "GDSArray")){
        name <- seed(data)@name
        perm <- seed(data)@permute
    } else{
        if(length(data) == nrow) {
            name <- paste0("annotation/", sub("_", "/", name))
        } else if(length(data) == ncol){
            pre <- ifelse(ff=="SEQ_ARRAY", "sample.annotation", "sample.annot")
            name <- paste(pre, name, sep="/")
        }
        perm <- FALSE
        ## convert XList into array with dim.
        if(is(data, "List"))
            data <- unlist(lapply(data, function(x) paste(x, collapse=",")))
        data <- as.array(data)
    }
    datapath <- strsplit(name, "/")[[1]]  ## gds node path.
    directories <- head(datapath, -1)
    dirnode <- index.gdsn(gfile, "")
    for (directory in directories) {
        if (!directory %in% ls.gdsn(index.gdsn(dirnode, ""))) {
            dirnode <- addfolder.gdsn(dirnode, name=directory, type="directory")
        } else {
            dirnode <- index.gdsn(dirnode, directory)
        }
    }
    ## FIXME: if a matrix / array, do this in 'chunks', e.g., of 10,000 rows
    idx <- seq(1, nrow(data), by = chunk_size)
    for (start in idx) {
        ridx = (start - 1) +
            seq_len( min(start + chunk_size, nrow(data) + 1) - start )
        if(length(dim(data)) == 1){
            value = unname(as.array(data[ridx, drop=FALSE]))
        }else if(length(dim(data)) == 2){
            value = unname(as.array(data[ridx,, drop=FALSE]))
        }else if(length(dim(data)) == 3){
            value = unname(as.array(data[ridx,,, drop=FALSE]))
        }
        if (perm) value = aperm(value)
        if(start == 1){
            datanode <- add.gdsn(dirnode, tail(datapath, 1),
                                 val = value, compress = "LZMA_RA", check=T)
        } else {
            datanode <- index.gdsn(gfile, name)
            append.gdsn(datanode, val=value, check=T) ##?? FIXME?
        }
    }
}
## .annodata_ondisk <- function(gds_path, name){
##     gfile <- openfn.gds(gds_path)
##     on.exit(closefn.gds(gfile))
##     ## anno.id <- read.gdsn(index.gdsn(gfile, paste0(rowORsample, ".id")))
##     ## node <- ifelse(rowORsample == "row", paste0("rowData/", name), paste0("colData/", name))
##     seed <- new("GDSArraySeed",
##                 file=gds_path,
##                 name=node,
##                 dim=.get_gdsdata_dim(gfile, node),
##                 dimnames=list(anno.id),  ## for snp/variant nodes only.
##                 permute=FALSE,
##                 first_val="ANY")
##     DelayedArray(seed)  ## return a DelayedArray/GDSArray object without names.
## }

## .shorten_gds_paths <- function(assays){
##     nassay <- length(assays)
##     for (i in seq_len(nassay)) {
##         a <- assays[[i]]
##         a@seed@file <- basename(a@seed@file)
##         assays[[i]] <- a
##     }
##     assays
## }

#' saveGDSSummarizedExperiment
#' Save all the assays in GDS format, including in-memory assays. Delayed assays with delayed operations on them are realized while they are written to disk.
#' @param x A SummarizedExperiment object, with the array data being ordinary array structure.
#' @param dir The directory to save the gds format of the array data, and the newly generated SummarizedExperiment object with array data in GDSArray format.
#' @param replace Whether to replace the directory if it already exists. The default is FALSE.
#' @param rowDataOnDisk whether to save the \code{rowData} as DelayedArray object. The default is TRUE.
#' @param colDataOnDisk whether to save the \code{colData} as DelayedArray object. The default is TRUE.
#' @param verbose whether to print the process messages. The default is FALSE.
#' @importFrom SummarizedExperiment assays assay "assays<-"
#' @export
saveGDSSummarizedExperiment <- function(se, dir="my_gds_se", replace=FALSE, rowDataOnDisk=TRUE, colDataOnDisk=TRUE, verbose=FALSE){
    if (!is(se, "SummarizedExperiment"))
        stop("'se' must be a SummarizedExperiment object")
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory where to save the ", class(se),
                  " object (the directory will be created)"))
    if (!isTRUEorFALSE(replace))
        stop("'replace' must be TRUE or FALSE")
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    
    .create_dir(dir, replace)
    gds_path <- file.path(dir, "assays.gds")
    ff <- .get_gdsdata_fileFormat(gdsfile(se)[1])

    ## write and save se data as gds format.
    if(ff == "SEQ_ARRAY") .create_seqgds(se, gds_path)
    if(ff == "SNP_ARRAY") .create_snpgds(se, gds_path)
    ## error: unable to find an inherited method for function ‘rowRanges’ for signature ‘"RangedSummarizedExperiment"
   
    for (elt in seq_along(assays(se))){
        if(verbose) message("Assay ", elt, "/", length(assays(se)), ": ",
                            names(assays(se))[elt] )
        .write_sedata_as_gds(assays(se)[[elt]], elt, ff, gds_path, chunk_size, nrow(se), ncol(se))
    }
    rnms <- names(rowData(se))
    if(ff == "SEQ_ARRAY") id.al <- match(c("REF", "ALT"), rnms) ## test needed.
    if(ff == "SNP_ARRAY") id.al <- grep("ALLELE", rnms) 
    for (elt in rnms[-id.al]){
        if(verbose) message( "rowData: ", elt )
        .write_sedata_as_gds(rowData(se)[[elt]], elt, ff, gds_path, chunk_size, nrow(se), ncol(se))
    }
    for (elt in names(colData(se))){
        if(verbose) message( "colData: ", elt )
        .write_sedata_as_gds(colData(se)[[elt]], elt, ff, gds_path, chunk_size, nrow(se), ncol(se))
    }

    ## todo: refer to makeSEfromGDS:: .
    colData <- .colData_gdsdata(gds_path, ff, names(colData(se)), colDataOnDisk)
    rowRange <- .rowRanges_gdsdata(gds_path, ff, names(rowData(se)), rowDataOnDisk)
    ## todo: remove info columns for rowData(se)
    ## error: Invalid SNP GDS file: invalid dimension of 'genotype'.
    
    if(ff == "SEQ_ARRAY"){
        infocols <- .info_seqgds(file, infoColumns, rowDataOnDisk)
        mcols(rowRange) <- DataFrame(mcols(rowRange), infocols)
    }
    ## colData(se) <- colData
    
    
    nassay <- length(assays(se))
    namesAssay <- names(assays(se))
    ## namesAssay_new <- gsub("/", "_", namesAssay)

    ## save assay data as GDSArray, skip when assay is already GDSArray or on-disk.
    ## by default, the input SE should have in-memory array data for all assays.
    for (i in seq_len(nassay)){
        if(is.array(assays(se)[[i]]))
            assays(se)[[i]] <- GDSArray(gds_path, name=namesAssay[i])
        ## assays(se)[[i]] <- GDSArray(gds_path, name=namesAssay_new[i])
    }
    names(assays(se)) <- namesAssay

    ## if(rowDataOnDisk){
    ##     .write_annoData_into_gds(se, gds_path, "row", verbose)
    ##     res <- setNames(
    ##         lapply(names(rowData(se)), function(x)
    ##             .annodata_ondisk(gds_path, "row", name=x)),
    ##         names(rowData(se)))
    ##     resDF <- DataFrame(lapply(res, I))
    ##     rowData(se) <- resDF
    ## }
    
    ## if(colDataOnDisk){
    ##     .write_annoData_into_gds(se, gds_path, "col", verbose)
    ##     res <- setNames(
    ##         lapply(names(colData(se)), function(x)
    ##             .annodata_ondisk(gds_path, "sample", name=x)),
    ##         names(colData(se)))
    ##     resDF <- DataFrame(lapply(res, I))
    ##     colData(se) <- resDF
    ## }

    rds_path <- file.path(dir, "se.rds")
    ans <- se
    ## se@assays <- .shorten_gds_paths(se@assays)
    saveRDS(se, file=rds_path)
    invisible(ans)
}

.THE_EXPECTED_STUFF <- c(
    "a GDS-based SummarizedExperiment object previously ",
    "saved with saveGDSSummarizedExperiment()"
)

.stop_if_bad_dir <- function(dir)
    stop(wmsg("directory \"", dir, "\" does not seem to contain ",
              .THE_EXPECTED_STUFF))

#' loadGDSSummarizedExperiment
#' to load the GDS back-end SummarizedExperiment object into R console. 
#' @param dir The directory to save the gds format of the array data, and the newly generated SummarizedExperiment object with array data in GDSArray format.
#' @export
loadGDSSummarizedExperiment <- function(dir="my_gds_se")
{
    ## library(rhdf5)  # for h5ls()
    ## library(HDF5Array)  # for the HDF5Array class
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory containing ", .THE_EXPECTED_STUFF))
    gds_path <- file.path(dir, "assays.gds")
    rds_path <- file.path(dir, "se.rds")
    if (!file.exists(gds_path) || !file.exists(rds_path))
        .stop_if_bad_dir(dir)
    ## h5_content <- try(rhdf5::h5ls(gds_path), silent=TRUE)
    ## if (inherits(h5_content, "try-error"))
    ##     .stop_if_bad_dir(dir)
    ## h5_datasets <- h5_content[ , "name"]
    ans <- readRDS(rds_path)
    if (!is(ans, "SummarizedExperiment"))
        .stop_if_bad_dir(dir)
    for (i in seq_along(assays(ans))) {
        a <- assay(ans, i, withDimnames=FALSE)
        if (!is(a, "GDSArray"))
            ## if(!identical(basename(gdsfile(a)), "assays.gds"))
            ## || !(a@seed@name %in% h5_datasets))
            .stop_if_bad_dir(dir)
        ## a@seed@file <- file_path_as_absolute(file.path(dir, a@seed@file))
        ## assay(ans, i, withDimnames=FALSE) <- a
    }
    ans
}

