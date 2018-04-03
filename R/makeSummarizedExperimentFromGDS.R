#' @import GDSArray
#' @importFrom GenomicRanges GRanges ranges seqnames 
#' @importFrom IRanges IRanges
#' @import gdsfmt
#' @import SNPRelate
#' @import SeqArray
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
        resDF <- setNames(DataFrame(a1, a2), paste0("ALLELE", 1:2))
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

### delayedArray in each column in DF, have individual index, which does not save spaces almost......

#' @importMethodsFrom SeqArray info
.info_seqgds <- function(seqArrayFile, infoColumns, rowDataOnDisk){
    f <- seqOpen(seqArrayFile)
    on.exit(seqClose(f))
    infonodes <- ls.gdsn(index.gdsn(f, "annotation/info"))
    if(length(infoColumns) > 0){
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
    }else{
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
#' @import SeqArray
#' @importMethodsFrom DelayedArray sub
.rowRanges_gdsdata <- function(file, fileFormat, rowDataColumns, rowDataOnDisk){
    rr <- .granges_gdsdata(file, fileFormat)
    ## following code generates the mcols(SummarizedExperiment::rowRanges(se))
    if(!is.null(rowDataColumns)){
        if(length(rowDataColumns) > 0){
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
        }else{
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
    if (!is.null(colDataColumns)) {
        if (length(colDataColumns) > 0) {
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
        } else {
            colDataColumns <- showAvailable(file)$colDataColumns
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
            DelayedDataFrame(annot, row.names=sample.id)
        }
    } else {
        f <- openfn.gds(file)
        on.exit(closefn.gds(f))
        stopifnot(inherits(f, "gds.class"))
        sample.id <- read.gdsn(index.gdsn(f, "sample.id"))
        DelayedDataFrame(row.names=sample.id)
    }
}

#' ShowAvailable
#' 
#' The function to show the available entries for the arguments within \code{makeSummarizedExperimentFromGDS}
#' @param file the path to the gds.class file.
#' @param args the arguments in \code{makeSummarizedExperimentFromGDS}.
#' @export
showAvailable <- function(file, args=c("name", "rowDataColumns", "colDataColumns", "infoColumns")){
    ## check if character.
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the gds file where the dataset is located."))
    args <- match.arg(args, several.ok=TRUE)
    ff <- GDSArray:::.get_gdsdata_fileFormat(file)
    res <- list()
    if("name" %in% args){
        if(ff == "SNP_ARRAY"){
            assaynodes <- "genotype"
        }else if(ff == "SEQ_ARRAY"){
            assaynodes <- GDSArray:::.get_gdsdata_arrayNodes(file)
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

#' makeSummarizedExperimentFromGDS
#' 
#' Conversion of gds file into SummarizedExperiment.
#' @param file the path to the gds.class file.
#' @param name the components of the gds file that will be represented
#'     as \code{GDSArray} file.
#' @param rowDataColumns which columns of \code{rowData} to
#'     import. The default is ALL.
#' @param colDataColumns which columns of \code{colData} to
#'     import. The default is ALL.
#' @param infoColumns which columns of \code{infoColumns} to import. The default is
#'     ALL.
#' @param rowDataOnDisk whether to save the \code{rowData} as
#'     DelayedArray object. The default is TRUE.
#' @param colDataOnDisk whether to save the \code{colData} as
#'     DelayedArray object. The default is TRUE.
#' @importFrom tools file_path_as_absolute
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stats setNames
#' 
#' @examples
#' \dontrun{
#' file <- SNPRelate::snpgdsExampleFileName()
#' se <- makeSummarizedExperimentFromGDS(file)
#' rowData(se)
#' SummarizedExperiment::colData(se)
#' metadata(se)
#' showAvailable(file)
#' se1 <- makeSummarizedExperimentFromGDS(file, rowDataColumns=c("ALLELE"))
#' SummarizedExperiment::rowRanges(se1)

#' file <- SeqArray::seqExampleFileName(type="gds")
#' se <- makeSummarizedExperimentFromGDS(file)
#' names(assays(se))
#' showAvailable(file)
#' names <- showAvailable(file, "name")$name
#' rowdatacols <- showAvailable(file, "rowDataColumns")$rowDataColumns
#' coldatacols <- showAvailable(file, "colDataColumns")$colDataColumns
#' infocols <- showAvailable(file, "infoColumns")$infoColumns
#' se1 <- makeSummarizedExperimentFromGDS(
#' file,
#' name = names[2],
#' rowDataColumns = rowdatacols[1:3],
#' colDataColumns = coldatacols[1],
#' infoColumns = infocols[c(1,3,5,7)],
#' rowDataOnDisk = FALSE,
#' colDataOnDisk = FALSE)
#' assay(se1)
#' rowRanges(se1)
#' colData(se1)
#' print(object.size(se[TRUE, TRUE]), units="auto")
#' }
#' @export
#' 
makeSummarizedExperimentFromGDS <- function(file, name=NULL, rowDataColumns=character(), colDataColumns=character(), infoColumns=character(), rowDataOnDisk=TRUE, colDataOnDisk=TRUE){
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
    if(ff == "SEQ_ARRAY"){
        infocols <- .info_seqgds(file, infoColumns, rowDataOnDisk)
        mcols(rowRange) <- DelayedDataFrame(mcols(rowRange), infocols)
    }
    se <- VariantExperiment(
        assays = assays,
        colData = colData,
        rowRanges = rowRange)
}
