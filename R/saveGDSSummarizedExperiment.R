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

.create_seqgds <- function(se, gds_path, compress){
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
    if(all(c("REF", "ALT") %in% names(rowData(se)))){
        ## ref <- mcols(rowRanges(se))$REF
        ## alt <- mcols(rowRanges(se))$ALT
        ref <- rowData(se)$REF
        alt <- rowData(se)$ALT
        if(is(ref, "List"))  ## "DNAStringSet" in-memory way
            ref <- unlist(lapply(ref, function(x) paste(x, collapse=",")))
        if(is(alt, "List"))  ## "DNAStringSetList" in-memory way
            alt <- unlist(lapply(alt, function(x) paste(x, collapse=",")))
        add.gdsn(gfile, "allele", val = paste(ref, alt, sep=","),
                 compress=compress, closezip=TRUE)
    }
}

.create_snpgds <- function(se, gds_path, compress){
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
    if("ID" %in% names(rowData(se)))
        add.gdsn(gfile, "snp.rs.id", as.character(rowData(se)$ID),
                 compress=compress, closezip=TRUE)
    add.gdsn(gfile, "snp.chromosome", as.character(seqnames(SummarizedExperiment::rowRanges(se))),
             compress=compress, closezip=TRUE)
    add.gdsn(gfile, "snp.position", ranges(SummarizedExperiment::rowRanges(se))@start,
             compress=compress, closezip=TRUE) ## ??
    if(all(c("ALLELE1", "ALLELE2") %in% names(rowData(se))))
        add.gdsn(gfile, "snp.allele",
                 val = paste(## mcols(rowRanges(se))$ALLELE1,
                             ## mcols(rowRanges(se))$ALLELE2,
                             rowData(se)$ALLELE1,
                             rowData(se)$ALLELE2,
                             sep="/"),
                 compress=compress, closezip=TRUE
                 )  ## "," separated.
}
 
.write_sedata_as_gds <- function(data, name, ff, gds_path, chunk_size, nrow, ncol, compress) {
    gfile <- openfn.gds(gds_path, readonly=FALSE)
    on.exit(closefn.gds(gfile))

    if(is(data, "GDSArray")){
        name <- seed(data)@name
        perm <- seed(data)@permute
    } else{
        if(length(data) == nrow) {
            if(name %in% c("ID", "QUAL", "FILTER")) name <- tolower(name)
            name <- paste0("annotation/", sub("_", "/", name))
            ## sub for "info_" nodes.
        } else if(length(data) == ncol){
            pre <- ifelse(ff=="SEQ_ARRAY", "sample.annotation", "sample.annot")
            name <- paste(pre, name, sep="/")
        }
        perm <- TRUE  ## for im-memory assay data, write to gds as reverse dimension.
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
    ## if a matrix / array, do this in 'chunks', e.g., of 10,000 rows
    idx <- seq(1, nrow(data), by = chunk_size)
    for (start in idx) {
        ridx = (start - 1) +
            seq_len( min(start + chunk_size, nrow(data) + 1) - start )
        ## ridxl <- logical(NROW(data))
        ## ridxl[ridx] <- TRUE
        ## value <- array(data[ridxl], c(sum(ridxl), dim(data)[-1]))
        if(length(dim(data)) == 1){
            value <- unname(as.array(data[ridx, drop=FALSE]))
        }else if(length(dim(data)) == 2){
            value <- unname(as.array(data[ridx,, drop=FALSE]))
        }else if(length(dim(data)) == 3){
            value <- unname(as.array(data[ridx,,, drop=FALSE]))
        }
        value <- aperm(value)  ## permute because adding rows in chunks.
        if(start == 1){
            datanode <- add.gdsn(dirnode, tail(datapath, 1),
                                 val = value, compress = compress, check=T)
        } else {
            datanode <- index.gdsn(gfile, name)
            append.gdsn(datanode, val=value, check=T) ##?? FIXME?
        }
    }
    ## permute back if permute == FALSE
    if(!perm){
        datanode <- index.gdsn(gfile, name)
        readmode.gdsn(datanode)
        permdim.gdsn(datanode, rev(seq_along(dim(data))))
    }
    ## add attribute for "snp_array" genotype node. 
    if(name == "genotype"){
        ord <- ifelse(perm, "sample.order", "snp.order")
        put.attr.gdsn(datanode, ord)
    }
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
#' @param fileFormat File format for the output gds file. See details.
#' @param rowDataOnDisk whether to save the \code{rowData} as DelayedArray object. The default is TRUE.
#' @param colDataOnDisk whether to save the \code{colData} as DelayedArray object. The default is TRUE.
#' @param verbose whether to print the process messages. The default is FALSE.
#' @importFrom SummarizedExperiment colData rowRanges rowData assays assay "assays<-"
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

    if(is(assay(se, 1), "GDSArray"))
        fileFormat <- .get_gdsdata_fileFormat(gdsfile(se)[1])
    
    ## write and save se data as gds file.
    if(fileFormat == "SEQ_ARRAY") .create_seqgds(se, gds_path, compress)
    if(fileFormat == "SNP_ARRAY") .create_snpgds(se, gds_path, compress)
    ## error: unable to find an inherited method for function ‘rowRanges’ for signature ‘"RangedSummarizedExperiment"

    ###
    ## Write gds file from SE with all contents
    ###

    ### assays
    for (elt in seq_along(assays(se))){
        if(verbose) message("Assay ", elt, "/", length(assays(se)), ": ",
                            names(assays(se))[elt] )
        .write_sedata_as_gds(assays(se)[[elt]], names(assays(se))[elt], fileFormat, gds_path, chunk_size, nrow(se), ncol(se), compress)
    }
    ### rowRanges
    rnms <- names(rowData(se))
    if(fileFormat == "SEQ_ARRAY") id.al <- match(c("REF", "ALT"), rnms)
    ## generated already with ".create_seqgds"
    if(fileFormat == "SNP_ARRAY") id.al <- grep("ALLELE", rnms)
    ## generated already with ".create_snpgds"
    for (elt in rnms[-id.al]){
        if(verbose) message( "rowData: ", elt )
        .write_sedata_as_gds(rowData(se)[[elt]], elt, fileFormat, gds_path, chunk_size, nrow(se), ncol(se), compress)
    }
    ### colData
    for (elt in names(SummarizedExperiment::colData(se))){
        if(verbose) message( "colData: ", elt )
        .write_sedata_as_gds(colData(se)[[elt]], elt, fileFormat, gds_path, chunk_size, nrow(se), ncol(se), compress)
    }

    
    ###
    ## save the new SE with new GDS file. 
    ### 

    ### save assay data as GDSArray, skip when assay is already GDSArray or on-disk.
    ## nassay <- length(assays(se))
    ## namesAssay <- names(assays(se))
    for (i in seq_along(assays(se))){
        if(is.array(assays(se)[[i]]))
            assays(se)[[i]] <- GDSArray(gds_path, name=names(assays(se))[i])
    }
    names(assays(se)) <- names(assays(se))

    ### save GDSArray-based colData if previous in-memory
    if(colDataOnDisk & ! "GDSArray" %in% lapply(colData(se), class)){
        colData <- .colData_gdsdata(gds_path, fileFormat, names(colData(se)), colDataOnDisk)
        colData(se) <- colData
    }
    ### save GDSArray-based rowData if previous in-memory
    if(rowDataOnDisk & ! "GDSArray" %in% lapply(rowData(se), class)){
        rowDataColumns <- sub("[0-9]", "", names(rowData(se)))
        ## "snp_array", sub "ALLELE1/2" with "ALLELE"
        rowDataColumns <- rowDataColumns[!grepl("info_", rowDataColumns)]
        ## "seq_array", remove all "info_" columns.
        rowRange <- .rowRanges_gdsdata(gds_path, fileFormat, rowDataColumns, rowDataOnDisk)
        ## add "info_" columns for "seq_array".
        if(fileFormat == "SEQ_ARRAY"){
            infoColumns <- names(rowData(se))[grepl("info_", names(rowData(se)))]
            infoColumns <- sub("info_", "", infoColumns)
            infocols <- .info_seqgds(gds_path, infoColumns, rowDataOnDisk)
            mcols(rowRange) <- DataFrame(mcols(rowRange), infocols)
        }
        rowRanges(se) <- rowRange
    }
    ## save new se file in ".rds"
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
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory containing ", .THE_EXPECTED_STUFF))
    gds_path <- file.path(dir, "se.gds")
    rds_path <- file.path(dir, "se.rds")
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

