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

.initiate_seqgds <- function(se, gds_path, compress){
    gfile <- createfn.gds(filename=gds_path)
    on.exit(closefn.gds(gfile))
    put.attr.gdsn(gfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")
    dscp <- addfolder.gdsn(gfile, "description", replace=TRUE)
    put.attr.gdsn(dscp, "source.format", "GDSArray-based SummarizedExperiment")
    add.gdsn(gfile, "sample.id", colnames(se), compress=compress, closezip=TRUE)
    add.gdsn(gfile, "variant.id", as.integer(rownames(se)),
             compress=compress, closezip=TRUE)
    add.gdsn(gfile, "chromosome", as.character(seqnames(SummarizedExperiment::rowRanges(se))),
             compress=compress, closezip=TRUE)
    add.gdsn(gfile, "position", ranges(SummarizedExperiment::rowRanges(se))@start,
             compress=compress, closezip=TRUE) ## ??
    if (all(c("REF", "ALT") %in% names(rowData(se)))){
        ## ref <- mcols(rowRanges(se))$REF
        ## alt <- mcols(rowRanges(se))$ALT
        ref <- rowData(se)$REF
        alt <- rowData(se)$ALT
        if (is(ref, "List"))  ## "DNAStringSet" in-memory way
            ref <- unlist(lapply(ref, function(x) paste(x, collapse=",")))
        if (is(alt, "List"))  ## "DNAStringSetList" in-memory way
            alt <- unlist(lapply(alt, function(x) paste(x, collapse=",")))
        add.gdsn(gfile, "allele", val = paste(ref, alt, sep=","),
                 compress=compress, closezip=TRUE)
    }
}

.initiate_snpgds <- function(se, gds_path, compress){
    ## use "LZMA_RA" as default "compress" argument
    gfile <- createfn.gds(filename=gds_path)
    on.exit(closefn.gds(gfile))
    put.attr.gdsn(gfile$root, "FileFormat", "SNP_ARRAY")
    dscp <- addfolder.gdsn(gfile, "description", replace=TRUE)
    put.attr.gdsn(dscp, "source.format", "GDSArray-based SummarizedExperiment")
    ## "description" not required for "SNP_ARRAY" gds file.
    add.gdsn(gfile, "sample.id", colnames(se), compress=compress, closezip=TRUE)
    add.gdsn(gfile, "snp.id", as.integer(rownames(se)), compress=compress,
             closezip=TRUE, replace=T)
    if ("ID" %in% names(rowData(se)))
        add.gdsn(gfile, "snp.rs.id", as.character(rowData(se)$ID),
                 compress=compress, closezip=TRUE)
    add.gdsn(gfile, "snp.chromosome", as.character(seqnames(SummarizedExperiment::rowRanges(se))),
             compress=compress, closezip=TRUE)
    add.gdsn(gfile, "snp.position", ranges(SummarizedExperiment::rowRanges(se))@start,
             compress=compress, closezip=TRUE) ## ??
    if (all(c("ALLELE1", "ALLELE2") %in% names(rowData(se))))
        add.gdsn(gfile, "snp.allele",
                 val = paste(## mcols(rowRanges(se))$ALLELE1,
                             ## mcols(rowRanges(se))$ALLELE2,
                             rowData(se)$ALLELE1,
                             rowData(se)$ALLELE2,
                             sep="/"),
                 compress=compress, closezip=TRUE
                 )  ## "," separated.
}
 
.write_sedata_as_gdsnode <- function(data, name, ff, gds_path, chunk_size, nrow, ncol, compress) {
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
            pre <- ifelse(ff=="SEQ_ARRAY", "sample.annotation", "sample.annot")
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
            dirnode <- addfolder.gdsn(dirnode, name=directory, type="directory")
        } else {
            dirnode <- index.gdsn(dirnode, directory)
        }
    }
    ## if a matrix / array, do this in 'chunks', e.g., of 10,000 rows
    idx <- seq(1, nrow(data), by = chunk_size)
    for (start in idx) {
        ridx = (start - 1) +
            seq_len( min(start + chunk_size, nrow(data) + 1) - start )
        ## ridxl <- logical(NROW(data))
        ## ridxl[ridx] <- TRUE
        ## value <- array(data[ridxl], c(sum(ridxl), dim(data)[-1]))
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
                                 val = value, compress = compress, check=T)
        } else {
            datanode <- index.gdsn(gfile, name)
            append.gdsn(datanode, val=value, check=T) ##?? FIXME?
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
            add.gdsn(dirnode, "@data", val=rep(1L, nrow), storage="uint8", compress=compress, visible=FALSE)
        ## n <- .AddVar(storage.option, varGeno, "@data", storage="uint8", visible=FALSE)
    }
}

###
    ## Write gds file from SE with all contents
    ###

.write_se_as_gds <- function(se, fileFormat, gds_path, chunk_size, compress, verbose){
    ### assays
    for (elt in seq_along(assays(se))){
        if (verbose) message("Assay ", elt, "/", length(assays(se)), ": ",
                            names(assays(se))[elt] )
        .write_sedata_as_gdsnode(assays(se)[[elt]],
                                 names(assays(se))[elt],
                                 fileFormat,
                                 gds_path,
                                 chunk_size,
                                 nrow(se),
                                 ncol(se),
                                 compress)
    }
    ### rowRanges
    rnms <- names(rowData(se))
    if (fileFormat == "SEQ_ARRAY") id.al <- match(c("REF", "ALT"), rnms)
    ## generated already with ".initiate_seqgds"
    ## if (fileFormat == "SNP_ARRAY") id.al <- grep("ALLELE", rnms)
    if (fileFormat == "SNP_ARRAY") id.al <- match(c("ID", "ALLELE1","ALLELE2"), rnms)
    ## generated already with ".initiate_snpgds"
    for (elt in rnms[-id.al]){
        if (verbose) message( "rowData: ", elt )
        .write_sedata_as_gdsnode(
            rowData(se)[[elt]],
            elt,
            fileFormat,
            gds_path,
            chunk_size,
            nrow(se),
            ncol(se),
            compress)
    }
    ### colData
    for (elt in names(SummarizedExperiment::colData(se))){
        if (verbose) message( "colData: ", elt )
        .write_sedata_as_gdsnode(
            SummarizedExperiment::colData(se)[[elt]],
            elt,
            fileFormat,
            gds_path,
            chunk_size,
            nrow(se),
            ncol(se),
            compress)
    }
}

.write_se_as_newse <- function(se, gds_path, fileFormat, colDataOnDisk, rowDataOnDisk){
    ### save assay data as GDSArray. if assay is already GDSArray or on-disk, only change the "file" slot to be newly generated gds file path.
    for (i in seq_along(assays(se))){
        if (is(assays(se)[[i]], "DelayedArray")){
            gdsfile(seed(assays(se)[[i]])) <- gds_path
        }else {
        assays(se)[[i]] <- GDSArray(gds_path, name=names(assays(se))[i])
        }
    }
    names(assays(se)) <- names(assays(se))

    ## save GDSArray-based colData if previous in-memory. If colData already GDSArray, only change the "file" slot to be the newly generated gds file path.
    if (colDataOnDisk){
        if (all(vapply(SummarizedExperiment::colData(se),
                      function(x) is(x, "DelayedArray"), logical(1)))){
            for (i in seq_len(ncol(SummarizedExperiment::colData(se)))){
                gdsfile(seed(SummarizedExperiment::colData(se)[[i]])) <- gds_path
            }
        }else {
            coldata <- .colData_gdsdata(
                gds_path,
                fileFormat,
                names(SummarizedExperiment::colData(se)),
                colDataOnDisk)
            SummarizedExperiment::colData(se) <- coldata
        }
    }
    ## save GDSArray-based rowData if previous in-memory. If rowData already GDSArray, only change the "file" slot to be the newly generated gds file path.
    if (rowDataOnDisk){
        if (all(vapply(rowData(se), function(x) is(x, "DelayedArray"), logical(1)))){
            for (i in seq_len(ncol(rowData(se)))){
                gdsfile(seed(rowData(se)[[i]])) <- gds_path
            }
        }else {
            rowDataColumns <- sub("[0-9]", "", names(rowData(se)))
            ## "snp_array", sub "ALLELE1/2" with "ALLELE"
            rowDataColumns <- rowDataColumns[!grepl("info_", rowDataColumns)]
            ## "seq_array", remove all "info_" columns.
            rowRange <- .rowRanges_gdsdata(gds_path, fileFormat, rowDataColumns, rowDataOnDisk)
            ## add "info_" columns for "seq_array".
            if (fileFormat == "SEQ_ARRAY"){
                infoColumns <- names(rowData(se))[grepl("info_", names(rowData(se)))]
                infoColumns <- sub("info_", "", infoColumns)
                infocols <- .info_seqgds(gds_path, infoColumns, rowDataOnDisk)
                mcols(rowRange) <- DataFrame(mcols(rowRange), infocols)
            }
            rowRanges(se) <- rowRange
        }
    }
    se
}

#' saveGDSSummarizedExperiment
#' Save all the assays in GDS format, including in-memory assays. Delayed assays with delayed operations on them are realized while they are written to disk.
#' @param se A SummarizedExperiment object, with the array data being ordinary array structure.
#' @param dir The directory to save the gds format of the array data, and the newly generated SummarizedExperiment object with array data in GDSArray format.
#' @param replace Whether to replace the directory if it already exists. The default is FALSE.
#' @param fileFormat File format for the output gds file. See details.
#' @param compress the compression method for writing the gds file. The default is "LZMA_RA".
#' @param chunk_size The chunk size (number of rows) when reading GDSArray-based assays from input \code{se} into memory and then write into a new gds file. 
#' @param rowDataOnDisk whether to save the \code{rowData} as DelayedArray object. The default is TRUE.
#' @param colDataOnDisk whether to save the \code{colData} as DelayedArray object. The default is TRUE.
#' @param verbose whether to print the process messages. The default is FALSE.
#' @importFrom SummarizedExperiment colData "colData<-" rowRanges "rowRanges<-" rowData "rowData<-" assays assay "assays<-"
#' @export
#' @details If the input \code{SummarizedExperiment} object has GDSArray-based assay data, there is no need to specify the argument \code{fileFomat}. Otherwise, it takes values of \code{SEQ_ARRAY} for sequencing data or \code{SNP_ARRAY} SNP array data. 

saveGDSSummarizedExperiment <- function(se, dir="my_gds_se", replace=FALSE, fileFormat=NULL, compress="LZMA_RA", chunk_size=10000, rowDataOnDisk=TRUE, colDataOnDisk=TRUE, verbose=FALSE){
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
    gds_path <- file.path(dir, "se.gds")
    
    if (is(assay(se, 1), "DelayedArray"))
        fileFormat <- .get_gdsdata_fileFormat(gdsfile(se)[1])
    
    ## initiate gds file.
    if (fileFormat == "SEQ_ARRAY") .initiate_seqgds(se, gds_path, compress)
    if (fileFormat == "SNP_ARRAY") .initiate_snpgds(se, gds_path, compress)
    gds_path <- tools::file_path_as_absolute(gds_path)
    
    ## write se data into gds file.
    .write_se_as_gds(se, fileFormat, gds_path, chunk_size, compress, verbose)
    
    ## save the new SE paired with new GDS file. 
    se <- .write_se_as_newse(se, gds_path, fileFormat, colDataOnDisk, rowDataOnDisk)
    
    ## save new se file in ".rds"
    rds_path <- file.path(dir, "se.rds")
    ## se@assays <- .shorten_gds_paths(se@assays)
    saveRDS(se, file=rds_path)
    invisible(se)
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
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory containing ", .THE_EXPECTED_STUFF))
    gds_path <- tools::file_path_as_absolute(file.path(dir, "se.gds"))
    rds_path <- tools::file_path_as_absolute(file.path(dir, "se.rds"))
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

