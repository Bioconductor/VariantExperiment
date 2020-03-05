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

.initiate_seqgds <- function(ve, gds_path, compress)
{
    gfile <- createfn.gds(filename=gds_path)
    on.exit(closefn.gds(gfile))
    put.attr.gdsn(gfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")
    dscp <- addfolder.gdsn(gfile, "description", replace=TRUE)
    put.attr.gdsn(dscp, "source.format",
                  "GDSArray-based SummarizedExperiment")
    add.gdsn(gfile, "sample.id", colnames(ve), compress=compress,
             closezip=TRUE)
    add.gdsn(gfile, "variant.id", as.integer(rownames(ve)),
             compress=compress, closezip=TRUE)
    add.gdsn(gfile, "chromosome",
             as.character(seqnames(SummarizedExperiment::rowRanges(ve))),
             compress=compress, closezip=TRUE)
    add.gdsn(gfile, "position",
             ranges(SummarizedExperiment::rowRanges(ve))@start,
             compress=compress, closezip=TRUE) ## ??
    if (all(c("REF", "ALT") %in% names(rowData(ve)))){
          ref <- rowData(ve)$REF
        alt <- rowData(ve)$ALT
        if (is(ref, "List"))  ## "DNAStringSet" in-memory way
            ref <- unlist(lapply(ref, function(x) paste(x, collapse=",")))
        if (is(alt, "List"))  ## "DNAStringSetList" in-memory way
            alt <- unlist(lapply(alt, function(x) paste(x, collapse=",")))
        add.gdsn(gfile, "allele", val = paste(ref, alt, sep=","),
                 compress=compress, closezip=TRUE)
    }
}

.initiate_snpgds <- function(ve, gds_path, compress)
{
    ## use "LZMA_RA" as default "compress" argument
    gfile <- createfn.gds(filename=gds_path)
    on.exit(closefn.gds(gfile))
    put.attr.gdsn(gfile$root, "FileFormat", "SNP_ARRAY")
    dscp <- addfolder.gdsn(gfile, "description", replace=TRUE)
    put.attr.gdsn(dscp, "source.format", "GDSArray-based SummarizedExperiment")
    ## "description" not required for "SNP_ARRAY" gds file.
    add.gdsn(gfile, "sample.id", colnames(ve), compress=compress, closezip=TRUE)
    add.gdsn(gfile, "snp.id", as.integer(rownames(ve)), compress=compress,
             closezip=TRUE, replace=TRUE)
    if ("ID" %in% names(rowData(ve)))
        add.gdsn(gfile, "snp.rs.id", as.character(rowData(ve)$ID),
                 compress=compress, closezip=TRUE)
    add.gdsn(gfile, "snp.chromosome",
             as.character(seqnames(SummarizedExperiment::rowRanges(ve))),
             compress=compress, closezip=TRUE)
    add.gdsn(gfile, "snp.position",
             ranges(SummarizedExperiment::rowRanges(ve))@start,
             compress=compress, closezip=TRUE) ## ??
    if (all(c("ALLELE1", "ALLELE2") %in% names(rowData(ve))))
        add.gdsn(gfile, "snp.allele",
                 val = paste(rowData(ve)$ALLELE1,
                             rowData(ve)$ALLELE2,
                             sep="/"),
                 compress=compress, closezip=TRUE
                 )  ## "," separated.
}

.write_vedata_as_gdsnode <- function(data, name, ff, gds_path,
                                     chunk_size, nrow, ncol, compress)
{
    gfile <- openfn.gds(gds_path, readonly=FALSE)
    on.exit(closefn.gds(gfile))

    if (is(data, "GDSArray")){
        name <- seed(data)@name
        perm <- seed(data)@permute
    } else {
        if (length(data) == nrow) {
            if (name %in% c("ID", "QUAL", "FILTER")) name <- tolower(name)
            name <- paste0("annotation/", sub("_", "/", name))
            ## sub for "info_" nodes.
        } else if (length(data) == ncol){
            pre <- ifelse(ff=="SEQ_ARRAY", "sample.annotation",
                          "sample.annot")
            name <- paste(pre, name, sep="/")
        }
        perm <- TRUE  ## for im-memory assay data, write to gds as reverse dimension.
        ## convert XList into array with dim.
        if (is(data, "List"))
            data <- unlist(lapply(data, function(x) paste(x, collapse=",")))
        data <- as.array(data)
    }
    datapath <- strsplit(name, "/")[[1]]  ## gds node path.
    directories <- head(datapath, -1)
    dirnode <- index.gdsn(gfile, "")
    for (directory in directories) {
        if (!directory %in% ls.gdsn(index.gdsn(dirnode, ""))) {
            dirnode <- addfolder.gdsn(dirnode, name=directory,
                                      type="directory")
        } else {
            dirnode <- index.gdsn(dirnode, directory)
        }
    }
    ## if a matrix / array, do this in 'chunks', e.g., of 10,000 rows
    idx <- seq(1, nrow(data), by = chunk_size)
    for (start in idx) {
        ridx = (start - 1) +
            seq_len( min(start + chunk_size, nrow(data) + 1) - start )
        if (length(dim(data)) == 1){
            value <- unname(as.array(data[ridx, drop=FALSE]))
        }else if (length(dim(data)) == 2){
            value <- unname(as.array(data[ridx,, drop=FALSE]))
        }else if (length(dim(data)) == 3){
            value <- unname(as.array(data[ridx,,, drop=FALSE]))
        }
        value <- aperm(value)  ## permute because adding rows in chunks.
        if (start == 1){
            datanode <- add.gdsn(dirnode, tail(datapath, 1),
                                 val = value, compress = compress, check=TRUE)
        } else {
            datanode <- index.gdsn(gfile, name)
            append.gdsn(datanode, val=value, check=TRUE) ##?? FIXME?
        }
    }
    ## permute back if permute == FALSE
    datanode <- index.gdsn(gfile, name)
    if (!perm){
        readmode.gdsn(datanode)
        permdim.gdsn(datanode, rev(seq_along(dim(data))))
    }
    ## add attribute for "snp_array" genotype node. 
    if (name == "genotype"){
        ord <- ifelse(perm, "sample.order", "snp.order")
        ## genotype should be a matrix. 
        put.attr.gdsn(datanode, ord)
    }
    ## add attribute for "seq_array" genotype folder node.
    if (length(directories) == 1L & all(directories == "genotype")){
        put.attr.gdsn(dirnode, "VariableName", "GT")
        put.attr.gdsn(dirnode, "Description", "Genotype")
        if (! "@data" %in% ls.gdsn(dirnode))
            add.gdsn(dirnode, "@data", val=rep(1L, nrow),
                     storage="uint8", compress=compress,
                     visible=FALSE)
        ## n <- .AddVar(storage.option, varGeno, "@data", storage="uint8", visible=FALSE)
    }
}

###
## Write gds file from VE with all contents
###

.write_ve_as_gds <- function(ve, fileFormat, gds_path, chunk_size,
                             compress, verbose)
{
    ### assays
    for (elt in seq_along(assays(ve))){
        if (verbose)
            message("Assay ", elt, "/", length(assays(ve)), ": ",
                    names(assays(ve))[elt] )
        .write_vedata_as_gdsnode(assays(ve)[[elt]],
                                 names(assays(ve))[elt],
                                 fileFormat,
                                 gds_path,
                                 chunk_size,
                                 nrow(ve),
                                 ncol(ve),
                                 compress)
    }
    ### rowRanges
    rnms <- names(rowData(ve))
    if (fileFormat == "SEQ_ARRAY")
        id.al <- match(c("REF", "ALT"), rnms)
    ## generated already with ".initiate_seqgds"
    ## if (fileFormat == "SNP_ARRAY") id.al <- grep("ALLELE", rnms)
    if (fileFormat == "SNP_ARRAY")
        id.al <- match(c("ID", "ALLELE1","ALLELE2"), rnms)
    ## generated already with ".initiate_snpgds"
    for (elt in rnms[-id.al]){
        if (verbose) message( "rowData: ", elt )
        .write_vedata_as_gdsnode(
            rowData(ve)[[elt]],
            elt,
            fileFormat,
            gds_path,
            chunk_size,
            nrow(ve),
            ncol(ve),
            compress)
    }
    ### colData
    for (elt in names(SummarizedExperiment::colData(ve))){
        if (verbose) message( "colData: ", elt )
        .write_vedata_as_gdsnode(
            SummarizedExperiment::colData(ve)[[elt]],
            elt,
            fileFormat,
            gds_path,
            chunk_size,
            nrow(ve),
            ncol(ve),
            compress)
    }
}

.write_ve_as_newve <- function(ve, gds_path, fileFormat,
                               colDataOnDisk, rowDataOnDisk)
{
    ### save assay data as GDSArray. if assay is already GDSArray or
    ### on-disk, only change the "file" slot to be newly generated gds
    ### file path.
    for (i in seq_along(assays(ve))){
    ##     if (is(assays(ve)[[i]], "DelayedArray")){
    ##         gdsfile(seed(assays(ve)[[i]])) <- gds_path
    ##     }else {
        assays(ve)[[i]] <- GDSArray(gds_path, name=names(assays(ve))[i])
    ##     }
    }
    names(assays(ve)) <- names(assays(ve))

    ## save GDSArray-based colData if previous in-memory. If colData
    ## already GDSArray, only change the "file" slot to be the newly
    ## generated gds file path.
    if (colDataOnDisk){
        if (all(vapply(SummarizedExperiment::colData(ve),
                       function(x) is(x, "DelayedArray"), logical(1)))){
            for (i in seq_len(ncol(SummarizedExperiment::colData(ve)))){
                gdsfile(seed(SummarizedExperiment::colData(ve)[[i]])) <- gds_path
            }
        }else {
            coldata <- .colData_gdsdata(
                gds_path,
                fileFormat,
                names(SummarizedExperiment::colData(ve)),
                colDataOnDisk)
            SummarizedExperiment::colData(ve) <- coldata
        }
    }
    ## save GDSArray-based rowData if previous in-memory. If rowData already GDSArray, only change the "file" slot to be the newly generated gds file path.
    if (rowDataOnDisk){
        if (all(vapply(rowData(ve), function(x) is(x, "DelayedArray"), logical(1)))){
            for (i in seq_len(ncol(rowData(ve)))){
                gdsfile(seed(rowData(ve)[[i]])) <- gds_path
            }
        }else {
            rowDataColumns <- sub("[0-9]", "", names(rowData(ve)))
            ## "snp_array", sub "ALLELE1/2" with "ALLELE"
            rowDataColumns <- rowDataColumns[!grepl("info_", rowDataColumns)]
            ## "seq_array", remove all "info_" columns.
            rowRange <- .rowRanges_gdsdata(gds_path, fileFormat,
                                           rowDataColumns,
                                           rowDataOnDisk)
            ## add "info_" columns for "seq_array".
            if (fileFormat == "SEQ_ARRAY"){
                infoColumns <- names(rowData(ve))[grepl("info_",
                                                        names(rowData(ve)))]
                infoColumns <- sub("info_", "", infoColumns)
                infocols <- .info_seqgds(gds_path, infoColumns,
                                         rowDataOnDisk)
                mcols(rowRange) <- cbind(mcols(rowRange), infocols)
            }
            rowRanges(ve) <- rowRange
        }
    }
    ve
}

#' saveVariantExperiment Save all the assays in GDS format, including
#' in-memory assays. Delayed assays with delayed operations on them
#' are realized while they are written to disk.
#' @param ve A SummarizedExperiment object, with the array data being
#'     ordinary array structure.
#' @param dir The directory to save the gds format of the array data,
#'     and the newly generated SummarizedExperiment object with array
#'     data in GDSArray format. The default is temporary directory
#'     within the R session.
#' @param replace Whether to replace the directory if it already
#'     exists. The default is FALSE.
#' @param fileFormat File format for the output gds file. See details.
#' @param compress the compression method for writing the gds
#'     file. The default is "LZMA_RA".
#' @param chunk_size The chunk size (number of rows) when reading
#'     GDSArray-based assays from input \code{ve} into memory and then
#'     write into a new gds file.
#' @param rowDataOnDisk whether to save the \code{rowData} as
#'     DelayedArray object. The default is TRUE.
#' @param colDataOnDisk whether to save the \code{colData} as
#'     DelayedArray object. The default is TRUE.
#' @param verbose whether to print the process messages. The default
#'     is FALSE.
#' @return An \code{VariantExperiment} object with the new
#'     \code{gdsfile()} \code{ve.gds} as specified in \code{dir}
#'     argument.
#' @export
#' @details If the input \code{SummarizedExperiment} object has
#'     GDSArray-based assay data, there is no need to specify the
#'     argument \code{fileFomat}. Otherwise, it takes values of
#'     \code{SEQ_ARRAY} for sequencing data or \code{SNP_ARRAY} SNP
#'     array data.
#' @examples
#' gds <- SeqArray::seqExampleFileName("gds")
#' ve <- makeVariantExperimentFromGDS(gds)
#' gdsfile(ve)
#' ve1 <- subsetByOverlaps(ve, GRanges("22:1-48958933"))
#' ve1
#' gdsfile(ve1)
#' aa <- tempfile()
#' obj <- saveVariantExperiment(ve1, dir=aa, replace=TRUE)
#' obj
#' gdsfile(obj)

saveVariantExperiment <-
    function(ve, dir=tempdir(), replace=FALSE, fileFormat=NULL,
             compress="LZMA_RA", chunk_size=10000, rowDataOnDisk=TRUE,
             colDataOnDisk=TRUE, verbose=FALSE)
{
    if (!is(ve, "SummarizedExperiment"))
        stop("'ve' must be a SummarizedExperiment object")
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory where to save the ", class(ve),
                  " object (the directory will be created)"))
    if (!isTRUEorFALSE(replace))
        stop("'replace' must be TRUE or FALSE")
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    
    .create_dir(dir, replace)
    gds_path <- file.path(dir, "ve.gds")
    
    if (is(assay(ve, 1), "DelayedArray"))
        fileFormat <- GDSArray:::.get_gds_fileFormat(gdsfile(ve))
    
    ## initiate gds file.
    if (fileFormat == "SEQ_ARRAY") .initiate_seqgds(ve, gds_path, compress)
    if (fileFormat == "SNP_ARRAY") .initiate_snpgds(ve, gds_path, compress)
    gds_path <- tools::file_path_as_absolute(gds_path)
    
    ## write se data into gds file.
    .write_ve_as_gds(ve, fileFormat, gds_path, chunk_size, compress,
                     verbose)
    
    ## save the new VE paired with new GDS file. 
    ve <- .write_ve_as_newve(ve, gds_path, fileFormat, colDataOnDisk,
                             rowDataOnDisk)
    
    ## save new ve file in ".rds"
    rds_path <- file.path(dir, "ve.rds")
    ## ve@assays <- .shorten_gds_paths(ve@assays)
    saveRDS(ve, file=rds_path)
    ## invisible(ve)
    return(loadVariantExperiment(dir = dir))
}

.THE_EXPECTED_STUFF <- c(
    "a GDS-based SummarizedExperiment object previously ",
    "saved with saveVariantExperiment()"
)

.stop_if_bad_dir <- function(dir)
    stop(wmsg("directory \"", dir, "\" does not seem to contain ",
              .THE_EXPECTED_STUFF))

#' loadVariantExperiment to load the GDS back-end SummarizedExperiment
#' object into R console.
#' @param dir The directory to save the gds format of the array data,
#'     and the newly generated SummarizedExperiment object with array
#'     data in GDSArray format.
#' @return An \code{VariantExperiment} object.
#' @export
#' @examples
#' gds <- SeqArray::seqExampleFileName("gds")
#' ve <- makeVariantExperimentFromGDS(gds)
#' ve1 <- subsetByOverlaps(ve, GRanges("22:1-48958933"))
#' aa <- tempfile()
#' saveVariantExperiment(ve1, dir=aa, replace=TRUE)
#' loadVariantExperiment(dir = aa)

loadVariantExperiment <- function(dir=tempdir())
{
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory containing ", .THE_EXPECTED_STUFF))
    gds_path <- tools::file_path_as_absolute(file.path(dir, "ve.gds"))
    rds_path <- tools::file_path_as_absolute(file.path(dir, "ve.rds"))
    if (!file.exists(gds_path) || !file.exists(rds_path))
        .stop_if_bad_dir(dir)
    ans <- readRDS(rds_path)
    if (!is(ans, "SummarizedExperiment"))
        .stop_if_bad_dir(dir)
    for (i in seq_along(assays(ans))) {
        a <- assay(ans, i, withDimnames=FALSE)
        if (!is(a, "GDSArray"))
            .stop_if_bad_dir(dir)
    }
    ans
}

