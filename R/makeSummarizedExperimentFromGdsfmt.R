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

#' @import VariantAnnotation
.rowRanges_snpgds <- function(gdsfile){
    .granges_snpgds(
        gdsfile,
        ID = read.gdsn(index.gdsn(gdsfile, "snp.rs.id")),
        allele1 = DNAStringSet(.alleles_snpgds(gdsfile)$allele1), ## DNAStringSetList class
        ## usually minor? check snp array data. GWAS::.bim? 
        allele2 = DNAStringSet(.alleles_snpgds(gdsfile)$allele2))
    ## usually major?
}

###
## colData for samples
###

.colData_snpgds <- function(gdsfile){
    stopifnot(inherits(gdsfile, "gds.class"))
    stopifnot(inherits(gdsfile, "SNPGDSFileClass"))
    
    sample.id <- read.gdsn(index.gdsn(gdsfile, "sample.id"))
    node <- index.gdsn(gdsfile, "sample.annot", silent=TRUE)
    if (!is.null(node)){
        vars <- ls.gdsn(node)
    }else{
        vars <- character()
    }
    if (length(vars) > 0){
        annot <- lapply(vars, function(v) {
            read.gdsn(index.gdsn(gdsfile, paste0("sample.annot/", v)))
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
makeSummarizedExperimentFromGdsfmt <- function(file, name=NA, frompkg = c("SNPRelate", "SeqArray")){
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the gds file where the dataset is located."))
    file <- file_path_as_absolute(file)
    if (!isSingleStringOrNA(name))
        stop("'type' must be a single string or NA")

    if (is.na(name)){
        name <- "genotype"
        ## the default value for "name"? genotype for snpgds. ? for SeqArray? 
    } 
    assay <- setNames(list(GDSArray(file, name)), name)
    
    frompkg <- match.arg(frompkg)
    if(frompkg == "SNPRelate"){  ## snpgdsFileClass
        f <- snpgdsOpen(file)
        on.exit(snpgdsClose(f))
        colData <- .colData_snpgds(f)
        rowRange <- .rowRanges_snpgds(f)
        ## Irange with meta data (ID, allele1, allele2)
    }else if(frompkg == "SeqArray"){  ## SeqVarGDSClass
        f <- seqOpen(file)
        on.exit(seqClose(f))
        colData <- SeqArray::colData(f)
        rowRange <- SeqArray::rowRanges(f)
    }
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



