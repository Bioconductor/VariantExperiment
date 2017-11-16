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

#' @import VariantAnnotation
.rowRanges_gdsdata <- function(file, fileFormat, rowDataColumns){
    if(fileFormat == "SNP_ARRAY"){
        f <- snpgdsOpen(file)
        on.exit(snpgdsClose(f))
    }else if(fileFormat == "SEQ_ARRAY"){
        f <- seqOpen(file)
        on.exit(seqClose(f))
    }
    rowDataColumns <- toupper(rowDataColumns)
    if(length(rowDataColumns) == 0){
        idx.id <- idx.allele <- idx.ref <- idx.alt <- idx.qual <- idx.filter <- TRUE
    }else{
    idx.id <- any(grepl("ID", rowDataColumns))
    idx.allele <- any(grepl("ALLELE", rowDataColumns))
    idx.ref <- any(grepl("REF", rowDataColumns))
    idx.alt <- any(grepl("ALT", rowDataColumns))
    idx.qual <- any(grepl("QUAL", rowDataColumns))
    idx.filter <- any(grepl("FILTER", rowDataColumns))
    }
    if(inherits(f, "SNPGDSFileClass")){
        meta <- DataFrame(ID = read.gdsn(index.gdsn(f, "snp.rs.id")),
                          ALLELE1 = DNAStringSet(.alleles_snpgds(f)$allele1),
                          ALLELE2 = DNAStringSet(.alleles_snpgds(f)$allele2)
                          ## DNAStringSetList class
                          )
        if(!idx.allele) meta <- meta[-grep("ALLELE", names(meta))]
        if(!idx.id) meta <- meta[-grep("ID", names(meta))]
        idx.other <- rowDataColumns %in% c("ID", "ALLELE")
        if(any(!idx.other)){
            warning("The snp annotation of '",
                    tolower(paste(rowDataColumns[!idx.other], collapse = " , ")),
                    "' does not exist!")
        }
        rr <- .granges_snpgds(f)
        mcols(rr) <- meta
    }else if(inherits(f, "SeqVarGDSClass")){
        meta <- DataFrame(ID = seqGetData(f, "annotation/id"),
                          REF = SeqArray::ref(f),
                          ALT = SeqArray::alt(f),
                          QUAL = SeqArray::qual(f),
                          FILTER = SeqArray::filt(f))
        if(!idx.ref) meta <- meta[-grep("REF", names(meta))]
        if(!idx.alt) meta <- meta[-grep("ALT", names(meta))]
        if(!idx.qual) meta <- meta[-grep("QUAL", names(meta))]
        if(!idx.filter) meta <- meta[-grep("FILTER", names(meta))]
        if(!idx.id) meta <- meta[-grep("ID", names(meta))]
        idx.other <- rowDataColumns %in% c("ID", "REF", "ALT", "QUAL", "FILTER")
        if(any(!idx.other)){
            warning("The variant annotation of '",
                    tolower(paste(rowDataColumns[!idx.other], collapse = " , ")),
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
    ## stopifnot(inherits(f, "SNPGDSFileClass"))
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

#' makeSummarizedExperimentFromGdsfmt
#' 
#' Conversion of gds file into SummarizedExperiment
#' @param file the path to the gds.class file.
#' @param name the components of the gds file that will be represented as \code{GDSArray} file.
#' @param frompkg From which package the gds file is generated.
#' @export
#' 
makeSummarizedExperimentFromGdsfmt <- function(file, name=NA, rowDataColumns=character(), colDataColumns=character()){
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the gds file where the dataset is located."))
    file <- file_path_as_absolute(file)
    if (!isSingleStringOrNA(name))
        stop("'name' must be a single string or NA")

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
    ## if(ff == "SNP_ARRAY"){
    ##     f <- snpgdsOpen(file)
    ##     on.exit(snpgdsClose(f))
    ##     colData <- .colData_snpgds(f)
    ##     rowRange <- .rowRanges_snpgds(f)
    ## }else if(ff == "SEQ_ARRAY"){
    ##     f <- seqOpen(file)
    ##     on.exit(seqClose(f))
    ##     colData <- SeqArray::colData(f)
    ##     rowRange <- SeqArray::rowRanges(f)
    ## }
    SummarizedExperiment(assay = assay, colData=colData, rowRanges = rowRange)
}


## file <- SNPRelate::sngdsExampleFileName()
## gdsa <- GDSArray(file)
## f <- snpgdsOpen(file)
## colData <- .colData_snpgds(f)
## rowRange <- .rowRanges_snpgds(f)
## SummarizedExperiment(assay = gdsa, colData = colData, rowRanges = rowRange)

## file1 <- SeqArray::seqExampleFileName("gds")
## gdsa1 <- GDSArray(file1)
## f1 <- seqOpen(file1)
## colData1 <- .colData_snpgds(f1)
## rowRange1 <- .rowRanges_snpgds(f1)
## SummarizedExperiment(assay = gdsa1, colData = colData1, rowRanges = rowRange1)



