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
.write_gds_assays <- function(assays, gds_path, verbose, allow.duplicate){
    nassay <- length(assays)
    namesAssay <- names(assays)
    namesAssay <- sub("/", "_", namesAssay)
    for (i in seq_len(nassay)) {
        a <- assays[[i]]
        ## gds_name <- sprintf("assay%03d", i)
        if (verbose)
            message("Start writing assay ", i, "/", nassay, " to '",
                    gds_path, "':")
        gfile <- createfn.gds(filename=gds_path, allow.duplicate=allow.duplicate)
        on.exit(closefn.gds(gfile))
        put.attr.gdsn(t$root, "FileFormat", "SE_ARRAY")
        add.gdsn(gfile, "sample.id", val=)
        add.gdsn(gfile "row.id", val=)
        add.gdsn(gfile, namesAssay[i], val=a)
        ## a <- HDF5Array::writeHDF5Array(a, h5_path, h5_name, chunk_dim, level, verbose=verbose)
        if (verbose)
            message("Finished writing assay ", i, "/", nassay, " to '",
                    gds_path, "'.")
    }
}

.shorten_gds_paths <- function(assays){
    nassay <- length(assays)
    for (i in seq_len(nassay)) {
        a <- assays[[i]]
        a@seed@file <- basename(a@seed@file)
        assays[[i]] <- a
    }
    assays
}

### Save all the assays in HDF5 format, including in-memory assays.
### Delayed assays with delayed operations on them are realized while they
### are written to disk..
saveGdsfmtSummarizedExperiment <- function(x, dir="my_gds_se", replace=FALSE,
                                           allow.duplicate=FALSE,
                                           verbose=FALSE){
    if (!is(x, "SummarizedExperiment"))
        stop("'x' must be a SummarizedExperiment object")
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory where to save the ", class(x),
                  " object (the directory will be created)"))
    if (!isTRUEorFALSE(replace))
        stop("'replace' must be TRUE or FALSE")
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    
    ## We try library(HDF5Array) before deleting or creating directory 'dir'.
    ## That way if HDF5Array is not installed then we will stop without having
    ## made changes to the file system.
    library(HDF5Array)  # for writeHDF5Array()
    .create_dir(dir, replace)
  
    gds_path <- file.path(dir, "assays.gds")
    nassay <- length(assays)
    namesAssay <- names(assays)
    namesAssay_new <- sub("/", "_", namesAssay)
    .write_gds_assays(x@assays, gds_path, allow.duplicate, verbose)
    for (i in seq_len(nassay)){
        assays(x)[[i]] <- GDSArray(gds_path, name=namesAssay_new[i])
    }
    names(assays(x)) <- namesAssay
    ## todo: write the gdsfile "makeSummarizeExperimentFromGdsfmt" and save as. 
    rds_path <- file.path(dir, "se.rds")
    ans <- x
    x@assays <- .shorten_gds_paths(x@assays)
    saveRDS(x, file=rds_path)
    
    invisible(ans)
}

.THE_EXPECTED_STUFF <- c(
    "a GDS-based SummarizedExperiment object previously ",
    "saved with saveGdsfmtSummarizedExperiment()"
)

.stop_if_bad_dir <- function(dir)
    stop(wmsg("directory \"", dir, "\" does not seem to contain ",
              .THE_EXPECTED_STUFF))

### Does a lot of checking and tries to fail graciously if the content
### of 'dir' doesn't look as expected.
loadHDF5SummarizedExperiment <- function(dir="my_h5_se")
{
    library(rhdf5)  # for h5ls()
    library(HDF5Array)  # for the HDF5Array class
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory containing ", .THE_EXPECTED_STUFF))
    h5_path <- file.path(dir, "assays.h5")
    rds_path <- file.path(dir, "se.rds")
    if (!file.exists(h5_path) || !file.exists(rds_path))
        .stop_if_bad_dir(dir)
    h5_content <- try(rhdf5::h5ls(h5_path), silent=TRUE)
    if (inherits(h5_content, "try-error"))
        .stop_if_bad_dir(dir)
    h5_datasets <- h5_content[ , "name"]
    ans <- readRDS(rds_path)
    if (!is(ans, "SummarizedExperiment"))
        .stop_if_bad_dir(dir)
    for (i in seq_along(assays(ans))) {
        a <- assay(ans, i, withDimnames=FALSE)
        if (!is(a, "HDF5Array") || !identical(a@seed@file, "assays.h5") ||
            !(a@seed@name %in% h5_datasets))
            .stop_if_bad_dir(dir)
        a@seed@file <- file_path_as_absolute(file.path(dir, a@seed@file))
        assay(ans, i, withDimnames=FALSE) <- a
    }
    ans
}

