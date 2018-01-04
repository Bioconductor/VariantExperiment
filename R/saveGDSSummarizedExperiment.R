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

## do not return value, just rewrite the assay data into gds file.
.write_assay_as_gds <- function(se, nassay, namesAssay, gds_path, verbose){
    gfile <- createfn.gds(filename=gds_path)
    on.exit(closefn.gds(gfile))
    put.attr.gdsn(gfile$root, "FileFormat", "SE_ARRAY")
    add.gdsn(gfile, "sample.id", val=colnames(se))
    add.gdsn(gfile, "row.id", val=rownames(se))
    
    for (i in seq_len(nassay)) {
        ## only consider the in-memory array in assay data.
        a <- assays(se)[[i]]
        if (verbose)
            message("Start writing assay ", i, "/", nassay, " to '",
                    gds_path, "':")
        ## if(inherits(a, "GDSArray")){
        ##     a <- as.array(a)
        ## }
        add.gdsn(gfile, namesAssay[i], val=a)
        ## a <- HDF5Array::writeHDF5Array(a, h5_path, h5_name, chunk_dim, level, verbose=verbose)
        if (verbose)
            message("Finished writing assay ", i, "/", nassay, " to '",
                    gds_path, "'.")
    }
}

.write_rowData_into_gds <- function(se, gds_path, verbose){
    gfile <- openfn.gds(gds_path, readonly=FALSE)
    on.exit(closefn.gds(gfile))
    rdnode <- addfolder.gdsn(gfile, name="rowData", type="directory")
    ## lapply(rowData(se), function(x) add.gdsn(rdnode, val=x))
    rd <- rowData(se)
    for(i in seq_len(ncol(rd))){
        name <- names(rd)[i]
        a <- rd[[i]]
        if (verbose)
            message("Start writing rowData column ", i, "/", ncol(rd), " to '",
                    gds_path, "':")
        ## if(isS4(a))
        ## Error in add.gdsn(rdnode, name, val = a) :
        ## No support of the storage mode 'S4'.
        if(is(a, "XStringSetList"))  ## DNAStringSetList
            a <- sapply(a, function(x) paste(x, collapse="/"))
        else if(is(a, "XStringSet") | is(a, "AtomicList"))
            ## DNAStringSet|CharacterList is(a)
            a <- as.character(a)
        suppressWarnings(add.gdsn(rdnode, name, val=a, replace=TRUE))
        if (verbose)
            message("Finished writing rowData column ", i, "/", ncol(rd), " to '",
                    gds_path, "':")
    }
}

.write_annoData_into_gds <- function(se, gds_path, rowORcol, verbose){
    gfile <- openfn.gds(gds_path, readonly=FALSE)
    on.exit(closefn.gds(gfile))
    annonode <- addfolder.gdsn(gfile, name=paste0(rowORcol, "Data"),
                               type="directory", replace=TRUE)
    if(rowORcol == "row") data <- rowData(se) else data <- colData(se)
    for(i in seq_len(ncol(data))){
        name <- names(data)[i]
        a <- data[[i]]
        if (verbose)
            message("Start writing ", rowORcol, "Data column ", i, "/", ncol(data),
                    " to '", gds_path, "':")
        ## if(isS4(a))
        ## Error in add.gdsn(datanode, name, val = a) :
        ## No support of the storage mode 'S4'.
        if(is(a, "XStringSetList"))  ## DNAStringSetList
            a <- sapply(a, function(x) paste(x, collapse="/"))
        else if(is(a, "XStringSet") | is(a, "AtomicList"))
            ## DNAStringSet|CharacterList is(a)
            a <- as.character(a)
        suppressWarnings(add.gdsn(annonode, name, val=a, replace=TRUE))
        if (verbose)
            message("Finished writing rowData column ", i, "/", ncol(data), " to '",
                    gds_path, "':")
    }
}

.annodata_ondisk <- function(gdsfile, rowORsample, name){
    gfile <- openfn.gds(gdsfile)
    on.exit(closefn.gds(gfile))
    anno.id <- read.gdsn(index.gdsn(gfile, paste0(rowORsample, ".id")))
    node <- ifelse(rowORsample == "row", paste0("rowData/", name), paste0("colData/", name))
    seed <- new("GDSArraySeed",
                file=gdsfile,
                name=node,
                dim=.get_gdsdata_dim(gfile, node),
                dimnames=list(anno.id),  ## for snp/variant nodes only.
                permute=FALSE,
                first_val="ANY")
    DelayedArray(seed)  ## return a DelayedArray/GDSArray object without names.
}

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
    
    ## We try library(HDF5Array) before deleting or creating directory 'dir'.
    ## That way if HDF5Array is not installed then we will stop without having
    ## made changes to the file system.
    ## library(HDF5Array)  # for writeHDF5Array()
    .create_dir(dir, replace)
  
    gds_path <- file.path(dir, "assays.gds")
    
    nassay <- length(assays(se))
    namesAssay <- names(assays(se))
    namesAssay_new <- gsub("/", "_", namesAssay)

    ## write and save assay data as gds format.
    .write_assay_as_gds(se, nassay, namesAssay_new, gds_path, verbose)

    ## save assay data as GDSArray, skip when assay is already GDSArray or on-disk.
    ## by default, the input SE should have in-memory array data for all assays.
    for (i in seq_len(nassay)){
        if(is.array(assays(se)[[i]]))
            assays(se)[[i]] <- GDSArray(gds_path, name=namesAssay_new[i])
    }
    names(assays(se)) <- namesAssay

    if(rowDataOnDisk){
        .write_annoData_into_gds(se, gds_path, "row", verbose)
        res <- setNames(
            lapply(names(rowData(se)), function(x)
                .annodata_ondisk(gds_path, "row", name=x)),
            names(rowData(se)))
        resDF <- DataFrame(lapply(res, I))
        rowData(se) <- resDF
    }
    
    if(colDataOnDisk){
        .write_annoData_into_gds(se, gds_path, "col", verbose)
        res <- setNames(
            lapply(names(colData(se)), function(x)
                .annodata_ondisk(gds_path, "sample", name=x)),
            names(colData(se)))
        resDF <- DataFrame(lapply(res, I))
        colData(se) <- resDF
    }

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

