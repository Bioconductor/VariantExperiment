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
    ## read and add the $reference (4/5/21)
    add.gdsn(dscp, "reference", SeqArray::seqSummary(gdsfile(ve), verbose = FALSE)$reference)
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
    if ("snp.rs.id" %in% names(rowData(ve)))
        add.gdsn(gfile, "snp.rs.id", as.character(rowData(ve)$ID),
                 compress=compress, closezip=TRUE)
    add.gdsn(gfile, "snp.chromosome",
             as.character(seqnames(SummarizedExperiment::rowRanges(ve))),
             compress=compress, closezip=TRUE)
    add.gdsn(gfile, "snp.position",
             ranges(SummarizedExperiment::rowRanges(ve))@start,
             compress=compress, closezip=TRUE) ## ??
    if (all(c("snp.allele1", "snp.allele2") %in% names(rowData(ve))))
        add.gdsn(gfile, "snp.allele",
                 val = paste(rowData(ve)$snp.allele1,
                             rowData(ve)$snp.allele2,
                             sep="/"),
                 compress=compress, closezip=TRUE
                 )  ## "," separated.
}

.write_vedata_as_gdsnode <- function(data, name, ff, gds_path,
                                     chunk_size, nrow, ncol, compress)
{
    gfile <- openfn.gds(gds_path, readonly=FALSE)
    on.exit(closefn.gds(gfile))

    if (is(data, "DelayedArray")){
        name <- seed(data)@varname
    } else {  ## for row/colData that are in DataFrame format (in-memory)
        ## reserve the customization from SEQ_ARRAY and SNP_ARRAY. 
        name <- sub("\\.", "/", name)  ## e.g., annotation/id, annotation/qual, annotation/filter
        if (length(data) == nrow && grepl("info.", name)) {
            name <- paste0("annotation/", name)  ## e.g., annotation/info/AC, etc. 
        } else if (length(data) == ncol){ ## FIXME: need the prefix?? 
            pre <- ifelse(ff=="SEQ_ARRAY", "sample.annotation",
                          "sample.annot")
            name <- paste(pre, name, sep="/")  ## e.g., "sample.annotation/family"
        }
        ## convert XList into array with dim.
        if (is(data, "List"))
            data <- unlist(lapply(data, function(x) paste(x, collapse=",")))
        data <- as.array(data)
    }

    dirnode <- index.gdsn(gfile, "")
    datapath <- strsplit(name, "/")[[1]]  ## gds node path.
    if (length(datapath) > 1) {
        dirs <- head(datapath, -1)
        for (i in dirs) {
            if (!i %in% ls.gdsn(dirnode)) {
                dirnode <- addfolder.gdsn(dirnode, name=i,
                                          type="directory")
            } else {
                dirnode <- index.gdsn(dirnode, i)  ## e.g., annotation/format/DP/data, annotation/id
            }
        }
    }
    
    ## if a matrix / array, do this in 'chunks', e.g., of 10,000 rows
    to <- ifelse(length(dim(data)) ==1, length(data), ncol(data))
    idx <- seq(1, to, by = chunk_size)
    for (start in idx) {
        cidx = (start - 1) + seq_len( min(start + chunk_size, to + 1) - start )
        if (length(dim(data)) == 1){
            value <- unname(as.array(data[cidx, drop=FALSE]))
        }else if (length(dim(data)) == 2){
            value <- unname(as.array(data[,cidx, drop=FALSE]))
        }else if (length(dim(data)) == 3){
            value <- unname(as.array(data[,cidx,, drop=FALSE]))
        }
        if (start == 1){
            datanode <- add.gdsn(dirnode, tail(datapath, 1),
                                 val = value, compress = "", check = FALSE)
            ## check = FALSE to suppress warnings. see ?gdsfmt::add[/append].gdsn
        } else {
            datanode <- index.gdsn(dirnode, tail(datapath, 1))
            ndim <- length(dim(data))
            if (ndim > 2) { ## append.gdsn() only append on last dimension. 
                idx <- seq_len(ndim)
                perm1 <- c(idx[-2], 2)
                perm2 <- c(1, ndim, idx[-c(1,ndim)])
                permdim.gdsn(datanode, perm1)
                append.gdsn(datanode, val = aperm(value, perm1), check = FALSE)
                permdim.gdsn(datanode, perm2)
            } else {
                append.gdsn(datanode, val=value, check = FALSE)  
            }
        }
    }
    compression.gdsn(datanode, compress = compress)

    ## ## permute back if permute == FALSE
    ## datanode <- index.gdsn(dirnode, name)
    ## if (!perm){
    ##     readmode.gdsn(datanode)
    ##     permdim.gdsn(datanode, rev(seq_along(dim(data))))
    ## }
    
    ## ## add attributes
    ## if (name == "genotype/data") {  ## "SEQ_ARRAY" genotype data
    ##     put.attr.gdsn(datanode, "VariableName", "GT")
    ##     put.attr.gdsn(datanode, "Description", "Genotype")
    ##     ## if (! "@data" %in% ls.gdsn(datanode))
    ##     ##     add.gdsn(datanode, "@data", val=rep(1L, nrow),
    ##     ##              storage="uint8", compress=compress,
    ##     ##              visible=FALSE)
    ##     ## n <- .AddVar(storage.option, varGeno, "@data", storage="uint8", visible=FALSE)
    ## } else if (name == "genotype"){ ## "SNP_ARRAY" genotype DATA
    ##     ## ord <- ifelse(perm, "sample.order", "snp.order")
    ##     put.attr.gdsn(datanode, "snp.order") ## always in "snp.order" when saving from VE.
    ## }
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
    if (fileFormat == "SNP_ARRAY")
        id.al <- match(c("snp.rs.id", "snp.allele1", "snp.allele2"), rnms)
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

.write_ve_as_newve <- function(ve, gds_path, fileFormat, ftnode, smpnode,
                               colDataOnDisk, rowDataOnDisk)
{
    ### save assay data as GDSArray. 
    ### file path.
    for (i in seq_along(assays(ve))){
        assays(ve, withDimnames = FALSE)[[i]] <-
            GDSArray(gds_path, names(assays(ve))[i])
    }
    ## names(assays(ve)) <- names(assays(ve))

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
            coldata <- .colData_seqgds(gds_path, smpnode,
                                       colDataColumns = names(SummarizedExperiment::colData(ve)),
                                       colDataOnDisk = colDataOnDisk)
            SummarizedExperiment::colData(ve) <- coldata
        }
    }
    ## save GDSArray-based rowData if previous in-memory. If rowData
    ## already GDSArray, only change the "file" slot to be the newly
    ## generated gds file path.
    if (rowDataOnDisk){
        if (all(vapply(rowData(ve), function(x) is(x, "DelayedArray"), logical(1)))){
            for (i in seq_len(ncol(rowData(ve)))){
                gdsfile(seed(rowData(ve)[[i]])) <- gds_path
            }
        }else {
            ## gds_path has the same contents as in ve. so the row/colDataColumns should be the same. 
            rownodes <- showAvailable(gds_path, "rowDataColumns")[[1]]
            rowRange <- .rowRanges_seqgds(gds_path, ftnode, rownodes, rowDataOnDisk)
            ## add "info_" columns for "seq_array".
            if (fileFormat == "SEQ_ARRAY"){
                infoColumns <- showAvailable(gds_path, "infoColumns")[[1]]
                infocols <- .infoColumns_seqgds(gds_path, ftnode, infoColumns, rowDataOnDisk)
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
#' @param chunk_size The chunk size (number of columns) when reading
#'     GDSArray-based assays from input \code{ve} into memory and then
#'     write into a new gds file. Default is 1000. Can be modified to
#'     smaller value if chunk data is too big (e.g., when number of
#'     rows are large). 
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
             compress="LZMA_RA", chunk_size=1000, rowDataOnDisk=TRUE,
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
        fileFormat <- .get_gds_fileFormat(gdsfile(ve))
    
    ## initiate gds file.
    ## FIXME: need a function to initiate other types of gds file... 
    
    if (fileFormat == "SEQ_ARRAY") {
        .initiate_seqgds(ve, gds_path, compress)
        ftnode <- "variant.id"
        smpnode <- "sample.id"
        }
    if (fileFormat == "SNP_ARRAY") {
        .initiate_snpgds(ve, gds_path, compress)
        ftnode <- "snp.id"
        smpnode <- "sample.id"
    }
    gds_path <- tools::file_path_as_absolute(gds_path)
    
    ## write se data into gds file.
    .write_ve_as_gds(ve, fileFormat, gds_path, chunk_size, compress,
                     verbose)
    
    ## save the new VE paired with new GDS file. 
    ve <- .write_ve_as_newve(ve, gds_path, fileFormat, ftnode, smpnode,
                             colDataOnDisk, rowDataOnDisk)
    
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
#' ## ve <- makeVariantExperimentFromGDS(gds)
#' ## ve1 <- subsetByOverlaps(ve, GRanges("22:1-48958933"))
#' aa <- tempfile()
#' ## saveVariantExperiment(ve1, dir=aa, replace=TRUE)
#' ## loadVariantExperiment(dir = aa)

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

