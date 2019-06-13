#' @import GDSArray
#' @importFrom GenomicRanges GRanges ranges seqnames 
#' @importFrom IRanges IRanges
#' @import gdsfmt
#' @import SNPRelate
#' @rawNamespace import(SeqArray, except = c(colData, rowRanges))

.granges_gdsdata <- function(gdsfile, fileFormat, ...){
    if(fileFormat == "SNP_ARRAY"){
        f <- snpgdsOpen(gdsfile)
        on.exit(snpgdsClose(f))
        vid <- read.gdsn(index.gdsn(f, "snp.id"))
        chr <- read.gdsn(index.gdsn(f, "snp.chromosome"))
        pos <- read.gdsn(index.gdsn(f, "snp.position"))
        ## reflen <- nchar(seqGetData(x, "$ref"))  
        ## reflen[reflen < 1L] <- 1L
        reflen <- 1L  ## all snps, length == 1L
        gr <- GenomicRanges::GRanges(seqnames=chr,
                                     ranges=IRanges(start=pos, end=pos+reflen-1L),
                                     ...)
        names(gr) <- as.integer(vid)
    }else if(fileFormat == "SEQ_ARRAY"){
        f <- seqOpen(gdsfile)
        on.exit(seqClose(f))
        gr <- SeqArray::granges(f)
    }
    gr
}

.varnode_snpgds_inmem <- function(snpgdsfile, name){
    f <- openfn.gds(snpgdsfile)
    on.exit(closefn.gds(f))
    if(name %in% "id") node <- "snp.rs.id"
    if(name %in% "allele") node <- "snp.allele"
    res <- read.gdsn(index.gdsn(f, node))
    resDF <- setNames(DataFrame(res), toupper(name))
    if(node == "snp.allele"){
        a <- strsplit(res, split="/")
        ## genotype data, possible for >2 alt? if yes, use sub() here
        ## and use DNAStringSetList class for "allele2"
        a1 <- Biostrings::DNAStringSet(unlist(a)[c(TRUE, FALSE)])
        a2 <- Biostrings::DNAStringSet(unlist(a)[c(FALSE, TRUE)]) 
        resDF <- setNames(DataFrame(a1, a2), paste0("ALLELE", seq_len(2)))
    }
    resDF  ## returns a DataFrame with names.
}

.varnode_seqgds_inmem <- function(seqgdsfile, name){
    f <- seqOpen(seqgdsfile)
    on.exit(seqClose(f))
    if(name %in% "id") res <- seqGetData(f, paste0("annotation/", name))
    if(name %in% "ref") res <- ref(f)
    if(name %in% "alt") res <- alt(f)
    if(name %in% "qual") res <- qual(f)
    if(name %in% "filter") res <- filt(f)
    resDF <- setNames(DataFrame(res), toupper(name))
    resDF  ## returns a DataFrame with names. 
}

#' @importFrom methods new
.varnode_gdsdata_ondisk <- function(gdsfile, fileFormat, name){
    f <- openfn.gds(gdsfile)
    on.exit(closefn.gds(f))
    if(fileFormat == "SNP_ARRAY"){
        varid <- read.gdsn(index.gdsn(f, "snp.id"))
        if(name %in% "id") node <- "snp.rs.id"
        if(name %in% "allele") node <- "snp.allele"
    }else if(fileFormat == "SEQ_ARRAY"){
        varid <- read.gdsn(index.gdsn(f, "variant.id"))
        varnodes <- ls.gdsn(index.gdsn(f, "annotation"))
        if(name %in% varnodes) node <- paste0("annotation/", name)
        if(name %in% c("alt", "ref")) node <- "allele"
    }
    seed <- new("GDSArraySeed",
                file=gdsfile,
                name=node,
                dim=objdesp.gdsn(index.gdsn(f, node))$dim,
                dimnames=list(varid),  ## for snp/variant nodes only.
                permute=FALSE,
                first_val="ANY")
    GDSArray(seed)  ## return a DelayedArray/GDSArray object without names.
}

#' @importMethodsFrom SeqArray info
#' @import DelayedDataFrame
.info_seqgds <- function(seqArrayFile, infoColumns, rowDataOnDisk){
    f <- seqOpen(seqArrayFile)
    on.exit(seqClose(f))
    infonodes <- ls.gdsn(index.gdsn(f, "annotation/info"))
    if (is.null(infoColumns)) {
        infoColumns <- infonodes
    } else {
        idx <- toupper(infoColumns) %in% infonodes
        if(any(!idx)){
            warning("\n", 'The "infoColumns" argument of "',
                    paste(infoColumns[!idx], collapse = ", "),
                    '" does not exist!', "\n",
                    'Please use showAvailable(file, "infoColumns") ',
                    'to get the available columns for "infoColumns."', "\n")
        }
        infoColumns <- toupper(infoColumns[idx])
        if(length(infoColumns) == 0)
            infoColumns <- infonodes
    }
    infonodes <- paste0("annotation/info/", infoColumns)
    if(rowDataOnDisk){
        res <- lapply(infonodes, function(x)
            GDSArray(
                file=new("GDSArraySeed",
                         file=seqArrayFile,
                         name=x,
                         dim=objdesp.gdsn(index.gdsn(f, x))$dim,
                         dimnames=list(seqGetData(f, "variant.id")),
                         permute=FALSE,
                         first_val="ANY")))
        res1 <- DelayedDataFrame(lapply(res, I))
    }else{
        res1 <- SeqArray::info(f, info=infoColumns)
    }
    setNames(res1, paste0("info_", infoColumns))
    ## return a DataFrame with names.
    }

#' @importFrom Biostrings DNAStringSet
#' @import SNPRelate
#' @importMethodsFrom DelayedArray sub
.rowRanges_gdsdata <- function(file, fileFormat, rowDataColumns, rowDataOnDisk){
    rr <- .granges_gdsdata(file, fileFormat)
    ## following code generates the mcols(SummarizedExperiment::rowRanges(se))
    if (is.character(rowDataColumns) && length(rowDataColumns) == 0) { ## character(0)
        f <- openfn.gds(file)
        on.exit(closefn.gds(f))
        stopifnot(inherits(f, "gds.class"))
        sample.id <- read.gdsn(index.gdsn(f, "sample.id"))
        if (rowDataOnDisk) {
            resDF <- DelayedDataFrame(row.names=sample.id)
        } else {
            resDF <- DataFrame(row.names=sample.id)
        }
        ## mcols(rr) <- resDF
    } else {
        if (is.null(rowDataColumns)) {
            rowDataColumns <- tolower(showAvailable(file)$rowDataColumns)
        } else {
            idx.within <- toupper(rowDataColumns) %in%
                showAvailable(file)$rowDataColumns
            if(any(!idx.within)){
                warning('The snp annotation of "',
                        paste(rowDataColumns[!idx.within], collapse = ", "),
                        '" does not exist!', "\n",
                        'Please use showAvailable(file, "rowDataColumns") ',
                        'to get the available columns for "rowData."', "\n")
            }
            rowDataColumns <- tolower(rowDataColumns[idx.within])
            if(length(rowDataColumns)==0)
                rowDataColumns <- tolower(showAvailable(file)$rowDataColumns)
        }
        if(rowDataOnDisk){
            res <- setNames(
                lapply(rowDataColumns, function(x)
                    .varnode_gdsdata_ondisk(file, fileFormat, name=x)),
                toupper(rowDataColumns))
            resDF <- DelayedDataFrame(lapply(res, I))
            if("ALLELE" %in% names(resDF)){
                resDF$ALLELE1 <- sub("/.$", "", resDF$ALLELE)
                resDF$ALLELE2 <- sub("[TCGA]*/", "", resDF$ALLELE)
                resDF[["ALLELE"]] <- NULL 
            }
            if("REF" %in% names(resDF)){
                resDF$REF <- sub(",.*", "", resDF$REF)
            }
            if("ALT" %in% names(resDF)){
                resDF$ALT <- sub("[TCGA]*,", "", resDF$ALT)
            }
        }else{ ## rowDataOnDisk = FALSE...
            if(fileFormat == "SNP_ARRAY"){
                resDF <- DataFrame(lapply(rowDataColumns, function(x)
                    .varnode_snpgds_inmem(file, x)))
            }else if(fileFormat == "SEQ_ARRAY"){
                resDF <- DataFrame(lapply(rowDataColumns, function(x)
                    .varnode_seqgds_inmem(file, x)))
            }
        }
        mcols(rr) <- resDF
    }
    rr
}

.sampnode_gdsdata_ondisk <- function(file, fileFormat, name)
{
    pre <- ifelse(fileFormat == "SNP_ARRAY", "sample.annot",
           ifelse(fileFormat == "SEQ_ARRAY", "sample.annotation", NULL))
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    sampnodes <- ls.gdsn(index.gdsn(f, pre))
    if(name %in% sampnodes){
        node <- paste0(pre, "/", name)
    }else{
        node <- name
    }
    seed <- new(
        "GDSArraySeed",
        file=file,
        name=node,
        dim=objdesp.gdsn(index.gdsn(f, node))$dim,
        dimnames=list(read.gdsn(index.gdsn(f, "sample.id"))),
        ## for sample-related nodes only.
        permute=FALSE,
        first_val="ANY")
    GDSArray(seed)  ## returns Delayed/GDSArray, not DF.
}

###
## colData for samples
###

.colData_gdsdata <- function(file, fileFormat, colDataColumns, colDataOnDisk)
{
    if (is.character(colDataColumns) && length(colDataColumns) == 0) { ## character(0)
        f <- openfn.gds(file)
        on.exit(closefn.gds(f))
        stopifnot(inherits(f, "gds.class"))
        sample.id <- read.gdsn(index.gdsn(f, "sample.id"))
        if (colDataOnDisk) {
            DelayedDataFrame(row.names=sample.id)
        } else {
            DataFrame(row.names=sample.id)
        }
    } else {
        if (is.null(colDataColumns)) {
            colDataColumns <- showAvailable(file)$colDataColumns
        } else {
            idx.within <- colDataColumns %in% showAvailable(file)$colDataColumns
            if (any(!idx.within)) {
                warning("\n", 'The sample annotation of "',
                        paste(colDataColumns[!idx.within], collapse = ", "),
                        '" does not exist!', "\n",
                        'Please use showAvailable(file, "colDataColumns") ',
                        'to get the available columns for "colData."',
                        "\n")
                colDataColumns <- colDataColumns[idx.within]
                if (length(colDataColumns) == 0)
                    colDataColumns <- showAvailable(file)$colDataColumns
            }
        } 
        if (colDataOnDisk) {
            sample.id <- .sampnode_gdsdata_ondisk(
                file, fileFormat, "sample.id")
            annot <- setNames(
                lapply(colDataColumns, function(x)
                    .sampnode_gdsdata_ondisk(file, fileFormat, x)),
                colDataColumns)
            DelayedDataFrame(lapply(annot, I),
                             row.names=as.character(sample.id))
        } else {
            f <- openfn.gds(file)
            on.exit(closefn.gds(f))
            stopifnot(inherits(f, "gds.class"))
            sample.id <- read.gdsn(index.gdsn(f, "sample.id"))
            pre <- ifelse(fileFormat == "SNP_ARRAY", "sample.annot",
                   ifelse(fileFormat == "SEQ_ARRAY", "sample.annotation", NULL))
            node <- paste0(pre, "/", colDataColumns)
            annot <- lapply(node, function(x) read.gdsn(index.gdsn(f, x)))
            names(annot) <- colDataColumns
            DataFrame(annot, row.names=sample.id)
        }
    }
}
#' ShowAvailable
#' 
#' The function to show the available entries for the arguments within
#' \code{makeVariantExperimentFromGDS}
#' @name showAvailable
#' @rdname makeVariantExperimentFromGDS
#' @param file the path to the gds.class file.
#' @param args the arguments in
#'     \code{makeVariantExperimentFromGDS}.
#' @examples
#' ## snp gds file
#' gds <- SNPRelate::snpgdsExampleFileName()
#' showAvailable(gds)
#'
#' ## sequencing gds file
#' gds <- SeqArray::seqExampleFileName("gds")
#' showAvailable(gds)
#'
#' @importFrom IRanges CharacterList
#' @export
#' 
showAvailable <- function(file,
                          args=c("name", "rowDataColumns", "colDataColumns", "infoColumns")){
    ## check if character.
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the gds file where the dataset is located."))
    args <- match.arg(args, several.ok=TRUE)
    ff <- GDSArray:::.get_gdsdata_fileFormat(file)
    res <- CharacterList()
    if("name" %in% args){
        if(ff == "SNP_ARRAY"){
            assaynodes <- "genotype"
        }else if(ff == "SEQ_ARRAY"){
            assaynodes <- GDSArray:::.get_gdsdata_non1D_array(file)
        }
        res$name <- assaynodes
    }
    if("rowDataColumns" %in% args){
        if(ff == "SNP_ARRAY"){
            rdnodes <- c("ID", "ALLELE")
        }else if(ff == "SEQ_ARRAY"){
            rdnodes <- c("ID", "ALT", "REF", "QUAL", "FILTER")
        }
        res$rowDataColumns <- rdnodes
    }
    if (any(c("colDataColumns", "infoColumns") %in% args)) {
        f <- openfn.gds(file)
        on.exit(closefn.gds(f))
    }
    if ("colDataColumns" %in% args) {
        fdnode <- ifelse(ff == "SNP_ARRAY", "sample.annot",
                  ifelse(ff == "SEQ_ARRAY", "sample.annotation", NA))
        cdnodes <- ls.gdsn(index.gdsn(f, fdnode))
        res$colDataColumns <- cdnodes
    }
    if("infoColumns" %in% args && ff == "SEQ_ARRAY"){
        infonodes <- ls.gdsn(index.gdsn(f, "annotation/info"))
        res$infoColumns <- infonodes
    }
    res
}

#' makeVariantExperimentFromGDS
#' 
#' Conversion of gds file into SummarizedExperiment.
#' @param name the components of the gds file that will be represented
#'     as \code{GDSArray} file.
#' @param rowDataColumns which columns of \code{rowData} to
#'     import. The default is NULL to read in all variant annotation
#'     info.
#' @param colDataColumns which columns of \code{colData} to
#'     import. The default is NULL to read in all sample related
#'     annotation info.
#' @param infoColumns which columns of \code{infoColumns} to
#'     import. The default is NULL to read in all info columns.
#' @param rowDataOnDisk whether to save the \code{rowData} as
#'     DelayedArray object. The default is TRUE.
#' @param colDataOnDisk whether to save the \code{colData} as
#'     DelayedArray object. The default is TRUE.
#' @return An \code{VariantExperiment} object.
#' @importFrom tools file_path_as_absolute
## #' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stats setNames
#' @examples
#' file <- SNPRelate::snpgdsExampleFileName()
#' se <- makeVariantExperimentFromGDS(file)
#' rowData(se)
#' colData(se)
#' metadata(se)
#' ## Only read specific columns for feature annotation.
#' showAvailable(file)
#' se1 <- makeVariantExperimentFromGDS(file, rowDataColumns=c("ALLELE"))
#' SummarizedExperiment::rowRanges(se1)

#' file <- SeqArray::seqExampleFileName(type="gds")
#' se <- makeVariantExperimentFromGDS(file)
#' ## all assay data
#' names(assays(se))
#' showAvailable(file)
#'
#' ## only read specific columns for feature / sample annotation. 
#' names <- showAvailable(file, "name")$name
#' rowdatacols <- showAvailable(file, "rowDataColumns")$rowDataColumns
#' coldatacols <- showAvailable(file, "colDataColumns")$colDataColumns
#' infocols <- showAvailable(file, "infoColumns")$infoColumns
#' se1 <- makeVariantExperimentFromGDS(
#' file,
#' name = names[2],
#' rowDataColumns = rowdatacols[1:3],
#' colDataColumns = coldatacols[1],
#' infoColumns = infocols[c(1,3,5,7)],
#' rowDataOnDisk = FALSE,
#' colDataOnDisk = FALSE)
#' assay(se1)
#' 
#' ## the rowData(se1) and colData(se1) are now in DataFrame format 
#' rowData(se1)
#' colData(se1)

#' @export
#' 
makeVariantExperimentFromGDS <- function(file, name=NULL,
                                         rowDataColumns = NULL,
                                         colDataColumns = NULL,
                                         infoColumns = NULL,
                                         rowDataOnDisk = TRUE,
                                         colDataOnDisk = TRUE)
{
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the gds file where the dataset is located."))
    file <- tools::file_path_as_absolute(file)
    stopifnot(is.character(name) | is.null(name))
    ## if (!isSingleStringOrNA(name))
    ##     stop("'name' must be a single string or NA")
    if(!isTRUEorFALSE(colDataOnDisk))
        stop("`colDataOnDisk` must be logical.")
    if(!isTRUEorFALSE(rowDataOnDisk))
        stop("`rowDataOnDisk` must be logical.")
    ## check which extensive gds format? SNPGDSFileClass or seqVarGDSClass? 
    ff <- GDSArray:::.get_gdsdata_fileFormat(file)
    if(is.null(name)){
        if(ff == "SNP_ARRAY"){
            name <- "genotype"
        }else if(ff == "SEQ_ARRAY"){
            name <- showAvailable(file, "name")$name
        }
    }
    assays <- setNames(lapply(name, function(x) GDSArray(file, x)), name)
    colData <- .colData_gdsdata(file, ff, colDataColumns, colDataOnDisk)
    rowRange <- .rowRanges_gdsdata(file, ff, rowDataColumns, rowDataOnDisk)
    if ((is.null(infoColumns) || length(infoColumns) > 0) && ff == "SEQ_ARRAY") {
        infocols <- .info_seqgds(file, infoColumns, rowDataOnDisk)
        mcols(rowRange) <- DelayedDataFrame(mcols(rowRange), infocols)
    }
    se <- VariantExperiment(
        assays = assays,
        colData = colData,
        rowRanges = rowRange)
}
