## file <- SNPRelate::snpgdsExampleFileName()
## f <- openfn.gds(file) ## "gds.class"
## f <- snpgdsOpen(file) ## "SNPGDSFileClass", "gds.class"
## fileSeqArray <- SeqArray::seqExampleFileName("gds")  ## "SeqVarGDSClass"

#' @rdname GDSArray
#' @exportMethod gdsfile
setMethod("gdsfile", "SummarizedExperiment", function(x) {
    vapply(assays(x), gdsfile, character(1))
})

#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
.granges_gdsdata <- function(gdsfile, fileFormat, ...){
    if(fileFormat == "SNP_ARRAY"){
        f <- SNPRelate::snpgdsOpen(gdsfile)
        on.exit(SNPRelate::snpgdsClose(f))
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
        f <- SeqArray::seqOpen(gdsfile)
        on.exit(SeqArray::seqClose(f))
        gr <- SeqArray::granges(f)
    }
    gr
}

.alleles_snpgds <- function(snpgdsfile){
    ff <- .get_gdsdata_fileFormat(snpgdsfile)
    stopifnot(ff == "SNP_ARRAY")
    f <- SNPRelate::snpgdsOpen(snpgdsfile)
    on.exit(SNPRelate::snpgdsClose(f))
    als <- read.gdsn(index.gdsn(f, "snp.allele"))
    als.d <- data.frame(
        do.call(rbind, strsplit(als, split="/")), stringsAsFactors = FALSE)
    res <- DataFrame(Biostrings::DNAStringSet(als.d[,1]),
                     Biostrings::DNAStringSet(als.d[,2]))
    names(res) <- paste0("allele", 1:2)
    res
}

.varnode_gdsdata_ondisk <- function(gdsfile, fileFormat, name){
    f <- gdsfmt::openfn.gds(gdsfile)
    on.exit(gdsfmt::closefn.gds(f))
    if(fileFormat == "SNP_ARRAY"){
        varid <- read.gdsn(index.gdsn(f, "snp.id"))
        if(name %in% "id") node <- "snp.rs.id"
        if(name %in% "allele") node <- "snp.allele"
    }else if(fileFormat == "SEQ_ARRAY"){
        varid <- read.gdsn(index.gdsn(f, "variant.id"))
        varnodes <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(f, "annotation"))
        if(name %in% varnodes) node <- paste0("annotation/", name)
        if(name %in% c("alt", "ref")) node <- "allele"
    }
    seed <- new("GDSArraySeed",
                file=gdsfile,
                name=node,
                dim=.get_gdsdata_dim(f, node),
                dimnames=list(varid),  ## for snp/variant nodes only.
                permute=FALSE,
                first_val="ANY")
    DelayedArray(seed)  ## return a DelayedArray/GDSArray object without names.
}

### delayedArray in each column in DF, have individual index, which does not save spaces almost......

.info_seqgds <- function(seqArrayFile, infoColumns, rowDataOnDisk){
    f <- SeqArray::seqOpen(seqArrayFile)
    on.exit(SeqArray::seqClose(f))
    infonodes <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(f, "annotation/info"))
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
            DelayedArray(
                seed=new("GDSArraySeed",
                         file=seqArrayFile,
                         name=x,
                         dim=.get_gdsdata_dim(f, x), 
                         dimnames=list(SeqArray::seqGetData(f, "variant.id")),
                         permute=FALSE,
                         first_val="ANY")))
        res1 <- DataFrame(lapply(res, function(x)DataFrame(I(x))))
    }else{
        res1 <- SeqArray::info(f, info=infoColumns)
    }
    setNames(res1, paste0("info_", infoColumns))
    ## return a DataFrame with names.
    }

## put fmt array data into assays(se).  -- todo!
.fmt_seqarray <- function(seqArrayFile, fmt, fmtOnDisk){
    f <- SeqArray::seqOpen(seqArrayFile)
    on.exit(SeqArray::seqClose(f))
    fmtnodes <- gdsfmt::ls.gdsn(index.gdsn(f, "annotation/format"))
    res <- seqGetData(f, "annotation/format/DP")  ## any other format nodes except DP?
    res$length <- as.matrix(res$length)
    res$data <- t(res$data)
    if(fmtOnDisk){
        res$length <- DelayedArray(res$length)
        res$data <- DelayedArray(res$data)
    }
    if(length(fmt)>0){
        idx <- fmt %in% fmtnodes
        if(any(!idx)){
            warning("\n", 'The "fmt" argument of "',
                    paste(fmt[!idx], collapse = ", "),
                    '" does not exist!', "\n",
                    'Please use showAvailable(file, "fmt") ',
                    'to get the available columns for "fmt."', "\n")
        }
    }
    res      
}

## snpgds: snp.rs.id, snp.allele
## seqarray: allele, annotation/id, annotation/qual, annotation/filter.

#' @importFrom Biostrings DNAStringSet
#' @import SNPRelate
#' @import SeqArray
#' @importMethodsFrom DelayedArray sub
.rowRanges_gdsdata <- function(file, fileFormat, rowDataColumns, rowDataOnDisk){
    rr <- .granges_gdsdata(file, fileFormat)
    ## following code generates the mcols(SummarizedExperiment::rowRanges(se))
    if(!is.null(rowDataColumns)){
        if(length(rowDataColumns) > 0){
            ## meta <- DataFrame(
            ##     ID = read.gdsn(index.gdsn(f, "snp.rs.id")),
            ##     ALLELE1 = .alleles_snpgds(file)$allele1,
            ##     ALLELE2 = .alleles_snpgds(file)$allele2  ## DNAStringSet class
            ## )
            idx.within <- toupper(rowDataColumns) %in%
                showAvailable(file)$rowDataColumns
            if(any(!idx.within)){
                warning('The snp annotation of "',
                        paste(rowDataColumns[!idx.within], collapse = ", "),
                        '" does not exist!', "\n",
                        'Please use showAvailable(file, "colDataColumns") ',
                        'to get the available columns for "colData."', "\n")
            }
            rowDataColumns <- tolower(rowDataColumns[idx.within])
            if(length(rowDataColumns)==0)
                rowDataColumns <- tolower(showAvailable(file)$rowDataColumns)
        }else{
            rowDataColumns <- tolower(showAvailable(file)$rowDataColumns)
        }
        ## if(rowDataOnDisk){
        res <- setNames(
            lapply(rowDataColumns, function(x)
                .varnode_gdsdata_ondisk(file, fileFormat, name=x)),
            toupper(rowDataColumns))
        resDF <- DataFrame(lapply(res, I))
        ## }else{ ## rowDataOnDisk = FALSE...
        ## res
        ## resDF 
        ## }
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
        mcols(rr) <- resDF
    }
    rr
}

.sampnode_gdsdata_ondisk <- function(file, fileFormat, name){
    pre <- ifelse(fileFormat == "SNP_ARRAY", "sample.annot",
           ifelse(fileFormat == "SEQ_ARRAY", "sample.annotation", NULL))
    f <- gdsfmt::openfn.gds(file)
    on.exit(gdsfmt::closefn.gds(f))
    sampnodes <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(f, pre))
    if(name %in% sampnodes){
        node <- paste0(pre, "/", name)
    }else{
        node <- name
    }
    seed <- new("GDSArraySeed",
                file=file,
                name=node,
                dim=.get_gdsdata_dim(f, node),
                dimnames=list(gdsfmt::read.gdsn(index.gdsn(f, "sample.id"))),
                ## for sample-related nodes only.
                permute=FALSE,
                first_val="ANY")
    da <- DelayedArray(seed)
    setNames(DataFrame(I(da), check.names=FALSE), name)
}

###
## colData for samples
###

.colData_gdsdata <- function(file, fileFormat, colDataColumns, colDataOnDisk){
    if(!is.null(colDataColumns)){
        if(length(colDataColumns) > 0){
            idx.within <- colDataColumns %in% showAvailable(file)$colDataColumns
            if(any(!idx.within)){
                warning("\n", 'The sample annotation of "',
                        paste(colDataColumns[!idx.within], collapse = ", "),
                        '" does not exist!', "\n",
                        'Please use showAvailable(file, "colDataColumns") ',
                        'to get the available columns for "colData."', "\n")
                colDataColumns <- colDataColumns[idx.within]
                if(length(colDataColumns) == 0)
                    colDataColumns <- showAvailable(file)$colDataColumns
            }
        }else{
            colDataColumns <- showAvailable(file)$colDataColumns
        }
        if(colDataOnDisk){
            sample.id <- .sampnode_gdsdata_ondisk(file, fileFormat, "sample.id")
            annot <- lapply(colDataColumns, function(x)
                .sampnode_gdsdata_ondisk(file, fileFormat, x))
            DataFrame(annot, row.names=as.character(sample.id$sample.id))
        }else{
            f <- gdsfmt::openfn.gds(file)
            on.exit(gdsfmt::closefn.gds(f))
            stopifnot(inherits(f, "gds.class"))
            sample.id <- read.gdsn(index.gdsn(f, "sample.id"))
            pre <- ifelse(fileFormat == "SNP_ARRAY", "sample.annot",
                   ifelse(fileFormat == "SEQ_ARRAY", "sample.annotation", NULL))
            node <- paste0(pre, "/", colDataColumns)
            annot <- lapply(node, function(x) read.gdsn(index.gdsn(f, x)))
            names(annot) <- colDataColumns
            DataFrame(annot, row.names=sample.id)
        }
    }else{
        f <- gdsfmt::openfn.gds(file)
        on.exit(gdsfmt::closefn.gds(f))
        stopifnot(inherits(f, "gds.class"))
        sample.id <- read.gdsn(index.gdsn(f, "sample.id"))
        ## DataFrame(matrix(0, nrow=length(sample.id), ncol=0), row.names=sample.id)
        DataFrame(row.names=sample.id)
    }
}

## for "DNAStringSetList", paste elements with "," for multiple alternative. (SE from "SeqArray", rowData column of "ALT", is "DNAStringSetList" format)
## .inmem_DF_to_DelayedArray_DF <- function(inmemdf){
##     classes <- sapply(inmemdf, class)
##     ind <- which(classes == "DNAStringSetList")
##     if(length(ind) > 0){
##         for(i in ind){
##             a <- inmemdf[[ind]]
##             t <- lapply(as.list(a), function(x)paste(x, collapse=","))
##             inmemdf[[ind]] <- unlist(t)
##         }
##     }
##     l <- lapply(inmemdf, function(x)DataFrame(I(DelayedArray(DataFrame(x)))))
##     dadf <- DataFrame(l)
##     names(dadf) <- names(inmemdf)
##     dadf
## }

#' ShowAvailable
#' 
#' The function to show the available entries for the arguments within \code{makeSummarizedExperimentFromGDS}
#' @param file the path to the gds.class file.
#' @param args the arguments in \code{makeSummarizedExperimentFromGDS}.
#' @export
showAvailable <- function(file, args=c("rowDataColumns", "colDataColumns", "infoColumns", "fmt")){
    ## check if character.
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the gds file where the dataset is located."))
    args <- match.arg(args, several.ok=TRUE)
    ff <- .get_gdsdata_fileFormat(file)
    f <- gdsfmt::openfn.gds(file)
    on.exit(gdsfmt::closefn.gds(f))
    res <- list()
    if("rowDataColumns" %in% args){
        if(ff == "SNP_ARRAY"){
            rdnodes <- c("ID", "ALLELE")
        }else if(ff == "SEQ_ARRAY"){
            rdnodes <- c("ID", "ALT", "REF", "QUAL", "FILTER")
        }
        res$rowDataColumns <- rdnodes
    }
    if("colDataColumns" %in% args){
        if(ff == "SNP_ARRAY"){
            cdnodes <- ls.gdsn(index.gdsn(f, "sample.annot"))
        }else if(ff == "SEQ_ARRAY"){
            cdnodes <- ls.gdsn(index.gdsn(f, "sample.annotation"))
        }
        res$colDataColumns <- cdnodes
    }
    if("infoColumns" %in% args){
        if(ff == "SNP_ARRAY"){
            infonodes <- NULL
        }else if(ff == "SEQ_ARRAY"){
            infonodes <- gdsfmt::ls.gdsn(index.gdsn(f, "annotation/info"))
            res$infoColumns <- infonodes
        }
    }
    if("fmt" %in% args){
        if(ff == "SNP_ARRAY"){
            fmtnodes <- NULL
        }else if(ff == "SEQ_ARRAY"){
            fmtnodes <- gdsfmt::ls.gdsn(index.gdsn(f, "annotation/format"))
            res$fmt <- fmtnodes
        }
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
#' @param fmt which data of \code{fmt} to import. The default is ALL.
#' @param rowDataOnDisk whether to save the \code{rowData} as
#'     DelayedArray object. The default is TRUE.
#' @param colDataOnDisk whether to save the \code{colData} as
#'     DelayedArray object. The default is TRUE.
#' @param fmtDataOnDisk whether to save the \code{fmtData} as
#'     DelayedArray object. The default is TRUE.
#' @importFrom tools file_path_as_absolute
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stats setNames
#' 
#' @examples
#' \dontrun{
#' file <- system.file(package="SNPRelate", "extdata", "hapmap_geno.gds")
#' se <- makeSummarizedExperimentFromGDS(file)
#' rowData(se)
#' SummarizedExperiment::colData(se)
#' metadata(se)
#' showAvailable(file)
#' se1 <- makeSummarizedExperimentFromGDS(file, rowDataColumns=c("ALLELE"))
#' SummarizedExperiment::rowRanges(se1)

#' file <- SeqArray::seqExampleFileName(type="gds")
#' se <- makeSummarizedExperimentFromGDS(file)
#' showAvailable(file)
#' rowdatacols <- showAvailable(file, "rowDataColumns")$rowDataColumns
#' coldatacols <- showAvailable(file, "colDataColumns")$colDataColumns
#' infocols <- showAvailable(file, "infoColumns")$infoColumns
#' fmt <- showAvailable(file, "fmt")$fmt
#' se1 <- makeSummarizedExperimentFromGDS(
#' file,
#' rowDataColumns=rowdatacols[1:3],
#' colDataColumns = coldatacols[1],
#' infoColumns = infocols[c(3,5,7)],
#' fmt = fmt
#' )
#' metadata(se1)
#' SummarizedExperiment::rowRanges(se1)
#' }
#' @export
#' 
makeSummarizedExperimentFromGDS <- function(file, name=NA, rowDataColumns=character(), colDataColumns=character(), infoColumns=character(), fmt=character(), rowDataOnDisk=TRUE, colDataOnDisk=TRUE, fmtDataOnDisk=TRUE){
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the gds file where the dataset is located."))
    file <- tools::file_path_as_absolute(file)
    if (!isSingleStringOrNA(name))
        stop("'name' must be a single string or NA")
    if(!isTRUEorFALSE(colDataOnDisk))
        stop("`colDataOnDisk` must be logical.")
    if(!isTRUEorFALSE(rowDataOnDisk))
        stop("`rowDataOnDisk` must be logical.")
    ## check which extensive gds format? SNPGDSFileClass or seqVarGDSClass? 
    ff <- .get_gdsdata_fileFormat(file)
    if(is.na(name)){
        name <- ifelse(
            ff == "SNP_ARRAY", "genotype",
                ifelse(ff == "SEQ_ARRAY", "genotype/data", NA))
    }
    assay <- setNames(list(GDSArray(file, name)), name)
    colData <- .colData_gdsdata(file, ff, colDataColumns, colDataOnDisk)
    rowRange <- .rowRanges_gdsdata(file, ff, rowDataColumns, rowDataOnDisk)
    if(ff == "SEQ_ARRAY"){
        infocols <- .info_seqgds(file, infoColumns, rowDataOnDisk)
        mcols(rowRange) <- DataFrame(mcols(rowRange), infocols)
        }
    ## on-disk representation of row/colData in DelayedArray back-end of DataFrame.
    ## if(rowDataOnDisk){
    ##     mcols(rowRange) <- .inmem_DF_to_DelayedArray_DF(mcols(rowRange))
    ## }
    ## if(colDataOnDisk){
    ##     colData <- .inmem_DF_to_DelayedArray_DF(colData)
    ## }
    se <- SummarizedExperiment(assays = assay, colData=colData, rowRanges = rowRange)
    ## save INFO as metadata for SeqVarGDSClass, on-disk by default. 
    ## if(ff == "SEQ_ARRAY"){
    ##     infodata <- .info_seqarray_ondisk(file, info=info)
    ##     metadata(se)$info <- infodata
    ##     fmtdata <- .fmt_seqarray(file, fmt=fmt, fmtOnDisk=fmtDataOnDisk)
    ##     metadata(se)$depthLength <- fmtdata$length
    ##     metadata(se)$depth <- fmtdata$data
    ## }
    se
}
