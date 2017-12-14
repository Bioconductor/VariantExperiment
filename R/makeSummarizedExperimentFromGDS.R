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

.varnode_snpgds_ondisk <- function(snpgdsFile, name){
    f <- SNPRelate::snpgdsOpen(snpgdsFile)
    on.exit(SNPRelate::snpgdsClose(f))
    ## varnodes <- gdsfmt::ls.gdsn(f)
    if(name %in% "id") node <- "snp.rs.id"
    if(name %in% "allele") node <- "snp.allele"
    seed <- new("GDSArraySeed",
                file=snpgdsFile,
                name=node,
                dim=.get_gdsdata_dim(f, node),
                dimnames=list(read.gdsn(index.gdsn(f, "snp.id"))),
                ## for snp nodes only.
                permute=FALSE,
                first_val="ANY")
    da <- DelayedArray(seed)
    setNames(DataFrame(I(da)), toupper(name))
    ##### todo: if "allele", then output a DataFrame with 2 columns of GDSArray.
}

.varnode_seqgds_ondisk <- function(seqArrayFile, name){
    f <- SeqArray::seqOpen(seqArrayFile)
    on.exit(SeqArray::seqClose(f))
    varnodes <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(f, "annotation"))
    if(name %in% varnodes) node <- paste0("annotation/", name)
    if(name %in% c("alt", "ref")) node <- "allele"
    seed <- new("GDSArraySeed",
                file=seqArrayFile,
                name=node,  ## only good for 1-dimension variant related nodes.
                            ## e.g., "annotation/id", "annotation/qual",
                            ## "annotation/filter", "ref", "alt".
                dim=.get_gdsdata_dim(f, node),
                dimnames=list(SeqArray::seqGetData(f, "variant.id")),
                ## for var nodes only.
                permute=FALSE,
                first_val="ANY")
    da <- DelayedArray(seed)
    if(name %in% "ref"){
        da <- sub(",.*", "", da)
        ## da <- as(sub(",.*", "", da), "GDSAarray")
    }else if(name %in% "alt"){
        da <- sub("[TCGA]*,", "", da)
        ## da <- as(sub("[TCGA]*,", "", da), "GDSArray")
    }
    setNames(DataFrame(I(da)), toupper(name)) ## maybe remove here?
}

.info_seqgds_ondisk <- function(seqArrayFile, info){
    f <- SeqArray::seqOpen(seqArrayFile)
    on.exit(SeqArray::seqClose(f))
    infonodes <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(f, "annotation/info"))
    if(length(info)>0){
        idx <- toupper(info) %in% infonodes
        if(any(!idx)){
            warning("\n", 'The "info" argument of "',
                    paste(info[!idx], collapse = ", "),
                    '" does not exist!', "\n",
                    'Please use showAvailable(file, "info") ',
                    'to get the available columns for "info."', "\n")
        }
        infonodes <- toupper(info[idx])
    }
    info <- paste0("annotation/info/", infonodes)
    res <- lapply(info, function(x)
        DelayedArray(
            seed=new("GDSArraySeed",
                     file=seqArrayFile,
                     name=x,
                     dim=.get_gdsdata_dim(f, x), 
                     ## objdesp.gdsn(index.gdsn(f, x))$dim,
                     dimnames=list(SeqArray::seqGetData(f, "variant.id")),
                     permute=FALSE,
                     first_val="ANY")))
    res1 <- lapply(res, function(x)DataFrame(I(x)))
    setNames(DataFrame(res1), paste0("info_", infonodes))
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
.rowRanges_gdsdata <- function(file, fileFormat, rowDataColumns, infoColumns
                               ## , rowDataOnDisk
                               ){
    rr <- .granges_gdsdata(file, fileFormat)
    ## following code generates the mcols(SummarizedExperiment::rowRanges(se))
    if(fileFormat == "SNP_ARRAY"){
        if(!is.null(rowDataColumns)){
            rowDataColumns <- tolower(rowDataColumns)
            if(length(rowDataColumns) > 0){
                ## meta <- DataFrame(
                ##     ID = read.gdsn(index.gdsn(f, "snp.rs.id")),
                ##     ALLELE1 = .alleles_snpgds(file)$allele1,
                ##     ALLELE2 = .alleles_snpgds(file)$allele2  ## DNAStringSet class
                ## )
                idx.within <- toupper(rowDataColumns) %in% c("ID", "ALLELE")
                if(any(!idx.within)){
                    warning('The snp annotation of "',
                            paste(rowDataColumns[!idx.within], collapse = ", "),
                            '" does not exist!',
                            'Please use showAvailable(file, "colDataColumns") ',
                            'to get the available columns for "colData."', "\n")
                }
            }else{
                rowDataColumns <- tolower(c("ID", "ALLELE"))
                res <- lapply(rowDataColumns, function(x)
                       .varnode_snpgds_ondisk(file, name=x))
                mcols(rr) <- res
        }
    }else if(fileFormat == "SEQ_ARRAY"){
        if(!is.null(rowDataColumns)){
            ## rowDataColumns <- tolower(rowDataColumns)
            if(length(rowDataColumns) > 0){
                idx.within <- toupper(rowDataColumns) %in%
                    c("ID", "REF", "ALT", "QUAL", "FILTER")
                if(any(!idx.within)){
                    warning('The variant annotation of "',
                            paste(rowDataColumns[!idx.within], collapse = ", "),
                            '" does not exist!', "\n",
                            'Please use showAvailable(file, "rowDataColumns") ',
                            'to get the available columns for "rowData."', "\n")
                    rowDataColumns <- tolower(rowDataColumns[idx.within])
                }
            }else{
                rowDataColumns <- tolower(c("ID", "REF", "ALT", "QUAL", "FILTER"))
            }
            ## if(rowDataOnDisk){
            res <- DataFrame(
                lapply(rowDataColumns, function(x)
                    .varnode_seqgds_ondisk(file, name=x)))
            mcols(rr) <- res
            if(!is.null(infoColumns)){
                info <- .info_seqgds_ondisk(file, infoColumns)
                mcols(rr) <- DataFrame(res, info)
                ## }
                ## }else{ ## rowDataOnDisk = FALSE...
            }
        }else{
            ## if(!is.null(infoColumns)){
            if(rowDataOnDisk){
                mcols(rr) <- .info_seqgds_ondisk(file, infoColumns)
                ##  }else{## rowDataOnDisk == FALSE...
            }
        }
    }
    }
    rr
}
    
###
## colData for samples
###

.colData_gdsdata <- function(file, fileFormat, colDataColumns){
    if(fileFormat == "SNP_ARRAY"){
        f <- snpgdsOpen(file)
        on.exit(snpgdsClose(f))
        node <- index.gdsn(f, "sample.annot", silent=TRUE)
    }else if(fileFormat == "SEQ_ARRAY"){
        f <- seqOpen(file)
        on.exit(seqClose(f))
        node <- index.gdsn(f, "sample.annotation", silent=TRUE)
    }
    stopifnot(inherits(f, "gds.class"))
    sample.id <- read.gdsn(index.gdsn(f, "sample.id"))

    if (!is.null(node)){
        vars <- ls.gdsn(node)
        if(length(colDataColumns) > 0){
            vars <- vars[vars %in% colDataColumns]
            idx <- colDataColumns %in% vars
            if(any(!idx)){
                warning("\n", 'The sample annotation of "',
                        paste(colDataColumns[!idx], collapse = ", "),
                        '" does not exist!', "\n",
                        'Please use showAvailable(file, "colDataColumns") ',
                        'to get the available columns for "colData."', "\n")
            }
        }
    }else{
        vars <- character()
    }
    if (length(vars) > 0){
        annot <- lapply(vars, function(v) {
            ## read.gdsn(index.gdsn(f, paste0("sample.annot/", v))
            read.gdsn(index.gdsn(node, v))
        })
        names(annot) <- vars
        DataFrame(Samples=seq_along(sample.id), annot, row.names=sample.id)
    } else {
        DataFrame(Samples=seq_along(sample.id), row.names=sample.id)
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
showAvailable <- function(file, args=c("rowDataColumns", "colDataColumns", "info", "fmt")){
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
    if("info" %in% args){
        if(ff == "SNP_ARRAY"){
            infonodes <- NULL
        }else if(ff == "SEQ_ARRAY"){
            infonodes <- gdsfmt::ls.gdsn(index.gdsn(f, "annotation/info"))
            res$info <- infonodes
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
#' @param infoColumns which columns of \code{info} to import. The default is
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
#' infocols <- showAvailable(file, "info")$info
#' fmt <- showAvailable(file, "fmt")$fmt
#' se1 <- makeSummarizedExperimentFromGDS(
#' file,
#' rowDataColumns=rowdatacols[1:3],
#' colDataColumns = coldatacols[1],
#' info = infocols[c(3,5,7)],
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
    colData <- .colData_gdsdata(file, fileFormat=ff, colDataColumns=colDataColumns)
    rowRange <- .rowRanges_gdsdata(file, fileFormat=ff, rowDataColumns = rowDataColumns, infoColumns = infoColumns
                                   ##, rowDataOnDisk = rowDataOnDisk
                                   )
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
