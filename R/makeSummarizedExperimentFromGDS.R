## file <- SNPRelate::snpgdsExampleFileName()
## f <- openfn.gds(file) ## "gds.class"
## f <- snpgdsOpen(file) ## "SNPGDSFileClass", "gds.class"
## fileSeqArray <- SeqArray::seqExampleFileName("gds")  ## "SeqVarGDSClass"

setMethod("gdsfile", "SummarizedExperiment", function(x) {
    vapply(assays(gse1), gdsfile, character(1))
})

#' @importFrom GenomicRanges GRanges
.granges_snpgds <- function(gdsfile, ...){
    stopifnot(inherits(gdsfile, "gds.class"))
    stopifnot(inherits(gdsfile, "SNPGDSFileClass"))
    vid <- read.gdsn(index.gdsn(gdsfile, "snp.id"))
    chr <- read.gdsn(index.gdsn(gdsfile, "snp.chromosome"))
    pos <- read.gdsn(index.gdsn(gdsfile, "snp.position"))
    ## reflen <- nchar(seqGetData(x, "$ref"))  
    ## reflen[reflen < 1L] <- 1L
    reflen <- 1L  ## all snps, length == 1L
    gr <- GenomicRanges::GRanges(seqnames=chr,
                                 ranges=IRanges(start=pos, end=pos+reflen-1L),
                  ...)
    names(gr) <- as.integer(vid)
    gr
}

.alleles_snpgds <- function(gdsfile){
    stopifnot(inherits(gdsfile, "gds.class"))
    stopifnot(inherits(gdsfile, "SNPGDSFileClass"))
    als <- read.gdsn(index.gdsn(gdsfile, "snp.allele"))
    als.d <- data.frame(
        do.call(rbind, strsplit(als, split="/")), stringsAsFactors = FALSE)
    names(als.d) <- paste0("allele", 1:2)
    als.d
}

## snpgds: snp.rs.id, snp.allele
## seqarray: allele, annotation/id, annotation/qual, annotation/filter.

#' @importFrom Biostrings DNAStringSet
#' @import SNPRelate
#' @import SeqArray
.rowRanges_gdsdata <- function(file, fileFormat, rowDataColumns){
    if(fileFormat == "SNP_ARRAY"){
        f <- SNPRelate::snpgdsOpen(file)
        on.exit(SNPRelate::snpgdsClose(f))
    }else if(fileFormat == "SEQ_ARRAY"){
        f <- SeqArray::seqOpen(file)
        on.exit(SeqArray::seqClose(f))
    }
    rowDataColumns <- toupper(rowDataColumns)
    if(length(rowDataColumns) == 0){
        idx.var <- setNames(rep(TRUE, 6), c("id", "allele", "ref", "alt", "qual", "filter")) 
        ## idx.id <- idx.allele <- idx.ref <- idx.alt <- idx.qual <- idx.filter <- TRUE
    }else{
        idx.var <- c(
            id = any(grepl("ID", rowDataColumns)),
            allele = any(grepl("ALLELE", rowDataColumns)),
            ref = any(grepl("REF", rowDataColumns)),
            alt = any(grepl("ALT", rowDataColumns)),
            qual = any(grepl("QUAL", rowDataColumns)),
            filter = any(grepl("FILTER", rowDataColumns))
        )
    }
    if(inherits(f, "SNPGDSFileClass")){
        meta <- DataFrame(ID = read.gdsn(index.gdsn(f, "snp.rs.id")),
                          ALLELE1 = Biostrings::DNAStringSet(.alleles_snpgds(f)$allele1),
                          ALLELE2 = Biostrings::DNAStringSet(.alleles_snpgds(f)$allele2)
                          ## DNAStringSetList class
                          )
        idx.var <- idx.var[1:2]
        if(any(!idx.var)){
            rm.var <- names(idx.var[!idx.var])
            meta <- meta[-match(toupper(rm.var), names(meta))]
        }
        idx.other <- rowDataColumns %in% c("ID", "ALLELE")
        if(any(!idx.other)){
            warning("The snp annotation of '",
                    tolower(paste(rowDataColumns[!idx.other], collapse = ", ")),
                    "' does not exist!")
        }
        rr <- .granges_snpgds(f)
        mcols(rr) <- meta
    }else if(inherits(f, "SeqVarGDSClass")){
        meta <- DataFrame(ID = SeqArray::seqGetData(f, "annotation/id"),
                          REF = SeqArray::ref(f),
                          ALT = SeqArray::alt(f),
                          QUAL = SeqArray::qual(f),
                          FILTER = SeqArray::filt(f))
        idx.var <- idx.var[-2]
        if(any(!idx.var)){
            rm.var <- names(idx.var[!idx.var])
            meta <- meta[-match(toupper(rm.var), names(meta))]
        }
        idx.other <- rowDataColumns %in% c("ID", "REF", "ALT", "QUAL", "FILTER")
        if(any(!idx.other)){
            warning("The variant annotation of '",
                    tolower(paste(rowDataColumns[!idx.other], collapse = ", ")),
                    "' does not exist!")
        }
        rr <- SeqArray::granges(f)
        mcols(rr) <- meta
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
                warning("\n", "The sample annotation of '",
                        paste(colDataColumns[!idx], collapse = " , "),
                        "' does not exist!", "\n")
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

.info_seqarray <- function(seqArrayFile, info){
    f <- SeqArray::seqOpen(seqArrayFile)
    on.exit(SeqArray::seqClose(f))
    mdata <- SeqArray::info(f)
    infonodes <- gdsfmt::ls.gdsn(index.gdsn(f, "annotation/info"))
    if(length(info)>0){
        if(any(! info %in% infonodes)){
            warning("the `info` argument should be NULL ",
                    "OR include only the following info: ",
                    paste(infonodes, collapse=", "))
        }
        mdata <- mdata[names(mdata) %in% info]
    }
    mdata
}


## for "DNAStringSetList", paste elements with "," for multiple alternative. (SE from "SeqArray", rowData column of "ALT", is "DNAStringSetList" format)
.inmem_DF_to_DelayedArray_DF <- function(inmemdf){
    classes <- sapply(inmemdf, class)
    ind <- which(classes == "DNAStringSetList")
    if(length(ind) > 0){
        for(i in ind){
            a <- inmemdf[[ind]]
            t <- lapply(as.list(a), function(x)paste(x, collapse=","))
            inmemdf[[ind]] <- unlist(t)
        }
    }
    l <- lapply(inmemdf, function(x)DataFrame(I(DelayedArray(DataFrame(x)))))
    dadf <- DataFrame(l)
    names(dadf) <- names(inmemdf)
    dadf
}

#' makeSummarizedExperimentFromGDS
#' 
#' Conversion of gds file into SummarizedExperiment
#' @param file the path to the gds.class file.
#' @param name the components of the gds file that will be represented as \code{GDSArray} file.
#' @param rowDataColumns which columns of \code{rowData} to import.
#' @param colDataColumns which columns of \code{colData} to import.
#' @param rowDataOnDisk whether to save the \code{rowData} as DelayedArray object. The default is TRUE.
#' @param colDataOnDisk whether to save the \code{colData} as DelayedArray object. The default is TRUE. 
#' @importFrom tools file_path_as_absolute
#' @export
#' 
makeSummarizedExperimentFromGDS <- function(file, name=NA, rowDataColumns=character(), colDataColumns=character(), info=character(), rowDataOnDisk=TRUE, colDataOnDisk=TRUE){
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
    rowRange <- .rowRanges_gdsdata(file, fileFormat=ff, rowDataColumns = rowDataColumns)
    ## on-disk representation of row/colData in DelayedArray back-end of DataFrame.
    if(rowDataOnDisk){
        mcols(rowRange) <- .inmem_DF_to_DelayedArray_DF(mcols(rowRange))
    }
    if(colDataOnDisk){
        colData <- .inmem_DF_to_DelayedArray_DF(colData)
    }

    ## save INFO as metadata for SeqVarGDSClass
    if(ff == "SEQ_ARRAY"){
        metadata <- .info_seqarray(file, info=info)
    }else{
        metadata <- list()
    }
    SummarizedExperiment(assay = assay, colData=colData, rowRanges = rowRange, metadata = metadata)
}





