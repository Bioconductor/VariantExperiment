## In "makeVariantExperimentFromSNPGDS()", we have customized the "ID"
## (from "snp.rs.id" node if available, "ALLELE1" and "ALLELE2" (from
## "snp.allele" node if available) columns as part of the
## mcols(rowRanges(ve)). "snp.id" was used as default sample id, and
## "sample.id" was used as default feature id. These 2 nodes are
## required by "SNPGDSFileClass" in SNPRelate.

#' @import GDSArray
#' @importFrom GenomicRanges GRanges ranges seqnames 
#' @importFrom IRanges IRanges
#' @import gdsfmt
#' @import SNPRelate

## .granges_gdsdata
.granges_snpgds <- function(snpgdsfile, ...){
    f <- snpgdsOpen(snpgdsfile)
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
    gr
}

## return a DataFrame
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

## #' @importFrom methods new

.varnode_snpgds <- function(snpgdsfile, name) {
    f <- openfn.gds(snpgdsfile)
    on.exit(closefn.gds(f))
    varid <- read.gdsn(index.gdsn(f, "snp.id"))
    if(name %in% "id") node <- "snp.rs.id"
    if(name %in% "allele") node <- "snp.allele"
    node
}
## .varnode_gdsdata_ondisk
.varnode_snpgds_ondisk <- function(snpgdsfile, name){
    node <- .varnode_snpgds(snpgdsfile, name)
    GDSArray(snpgdsfile, node)
}

.rowDataColumns_check <- function(gdsfile, rowDataColumns) {
    if (is.null(rowDataColumns)) {
        rowDataColumns <- tolower(showAvailable(gdsfile)$rowDataColumns)
    } else {
        idx.within <- toupper(rowDataColumns) %in%
            showAvailable(gdsfile)$rowDataColumns
        if(any(!idx.within)) {
            misrdcs <- paste(rowDataColumns[!idx.within], collapse = ", ")
            warning(.rowDataColumns_check_msg(misrdcs))
            rowDataColumns <- rowDataColumns[idx.within]
            ## if none, return character(0)
        }
    }
    tolower(rowDataColumns)
}

.rowDataColumns_check_msg <- function(missingRDCs) {
    paste('The feature (snp/variant) annotation of "',
          missingRDCs, 
          '" does not exist!', "\n",
          'Please use showAvailable(file, "rowDataColumns") ',
          'to get the available columns for "rowData()."', "\n")
}

## "snp.id", "variant.id"
.empty_rowData_DF <- function(gdsfile, ftnode, rowDataOnDisk) {
    f <- openfn.gds(gdsfile)
    on.exit(closefn.gds(f))
    ft.id <- read.gdsn(index.gdsn(f, ftnode))
    if (rowDataOnDisk) {
        DelayedDataFrame(row.names=ft.id)
    } else {
        DataFrame(row.names=ft.id)
    }
}

#' @importFrom Biostrings DNAStringSet
#' @import SNPRelate
#' @importMethodsFrom DelayedArray sub

## .rowRanges_gdsdata
.rowRanges_snpgds <- function(snpgdsfile, rowDataColumns, rowDataOnDisk){
    
    rr <- .granges_snpgds(snpgdsfile)
    ## following code generates the mcols(SummarizedExperiment::rowRanges(se))
    rowDataColumns <- .rowDataColumns_check(snpgdsfile, rowDataColumns)

    ## if no available rowDataColumns are selected, i.e.,
    ## rowDataColumns = character(0), return an empty (Delayed)DataFrame
    ## for mcols()
    if (is.character(rowDataColumns) && length(rowDataColumns) == 0) {
        resDF <- .empty_rowData_DF(snpgdsfile, "snp.id", rowDataOnDisk)
        mcols(rr) <- resDF
        return(rr)
    }
    ## if there are valid rowDataColumns
    if(rowDataOnDisk){
        res <- setNames(
            lapply(rowDataColumns, function(x)
                .varnode_snpgds_ondisk(snpgdsfile, name=x)),
            toupper(rowDataColumns))
        resDF <- DelayedDataFrame(lapply(res, I))
        if("ALLELE" %in% names(resDF)){ ## for snpgds
            resDF$ALLELE1 <- sub("/.$", "", resDF$ALLELE)
            resDF$ALLELE2 <- sub("[TCGA]*/", "", resDF$ALLELE)
            resDF[["ALLELE"]] <- NULL 
        }
    }else{ ## rowDataOnDisk = FALSE...
        resDF <- DataFrame(lapply(rowDataColumns, function(x)
            .varnode_snpgds_inmem(snpgdsfile, x)))
    }
    mcols(rr) <- resDF
    rr
}

.sampnodes_snpgds <- function(snpgdsfile) {
    f <- openfn.gds(snpgdsfile)
    on.exit(closefn.gds(f))
    ls.gdsn(index.gdsn(f, "sample.annot"))
    
}

## returns Delayed/GDSArray, not DF.
.sampnode_snpgds_ondisk <- function(snpgdsfile, name) {
    sampnodes <- .sampnodes_snpgds(snpgdsfile)
    if(name %in% sampnodes){
        node <- paste0("sample.annot/", name)
    }else{
        node <- name
    }
    GDSArray(snpgdsfile, node)
}

###
## colData for samples
###

.colDataColumns_check <- function(gdsfile, colDataColumns) {
    if (is.null(colDataColumns)) {
        colDataColumns <- showAvailable(gdsfile)$colDataColumns
    } else {
        idx.within <- colDataColumns %in%
            showAvailable(gdsfile)$colDataColumns
        if(any(!idx.within)) {
            misrdcs <- paste(colDataColumns[!idx.within], collapse = ", ")
            warning(.colDataColumns_check_msg(misrdcs))
            colDataColumns <- colDataColumns[idx.within]
            ## if none, return character(0)
        }
    }
    colDataColumns
}

.colDataColumns_check_msg <- function(missingCDCs) {
    paste("\n", 'The sample annotation of "',
          missingCDCs,
          '" does not exist!', "\n",
          'Please use showAvailable(file, "colDataColumns") ',
          'to get the available columns for "colData()."', "\n")
}

.empty_colData_DF <- function(gdsfile, smpnode, colDataOnDisk) {
    f <- openfn.gds(gdsfile)
    on.exit(closefn.gds(f))
    sample.id <- read.gdsn(index.gdsn(f, smpnode))
    if (colDataOnDisk) {
        DelayedDataFrame(row.names=sample.id)
    } else {
        DataFrame(row.names=sample.id)
    }
}

.colData_snpgds <- function(snpgdsfile, colDataColumns, colDataOnDisk) {
    colDataColumns <- .colDataColumns_check(snpgdsfile, colDataColumns)

    ## if no available colDataColumns are selected, i.e.,
    ## colDataColumns = character(0), return an empty (Delayed)DataFrame    
    if (is.character(colDataColumns) && length(colDataColumns) == 0) { ## character(0)
        .empty_colData_DF(snpgdsfile, "sample.id", colDataOnDisk)
    } else { ## if there are valid rowDataColumns
        if (colDataOnDisk) {
            sample.id <- .sampnode_snpgds_ondisk(snpgdsfile, "sample.id")
            annot <- setNames(
                lapply(colDataColumns, function(x)
                    .sampnode_snpgds_ondisk(snpgdsfile, x)),
                colDataColumns)
            DelayedDataFrame(lapply(annot, I),
                             row.names=as.character(sample.id))
        } else {
            f <- openfn.gds(snpgdsfile)
            on.exit(closefn.gds(f))
            sample.id <- read.gdsn(index.gdsn(f, "sample.id"))
            node <- paste0("sample.annot/", colDataColumns)
            annot <- lapply(node, function(x) read.gdsn(index.gdsn(f, x)))
            names(annot) <- colDataColumns
            DataFrame(annot, row.names=sample.id)
        }
    }
}

## http://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
## SNPRelate requires 5 variables: sample.id, snp.id, snp.position, snp.chromosome, genotype
## 2 optional variables: snp.rs.id, snp.allele (when merging gentypes from diff platforms)
.showAvailable_snparray <- function(file,
                                    args = c("assayNames",
                                             "rowDataColumns",
                                             "colDataColumns")) {
    res <- CharacterList()
    if("assayNames" %in% args)
        res$assayNames <- "genotype"
    if("rowDataColumns" %in% args) {
        res$rowDataColumns <- c("ID", "ALLELE")
    }
    if ("colDataColumns" %in% args) {
        f <- openfn.gds(file)
        on.exit(closefn.gds(f))
        res$colDataColumns <- ls.gdsn(index.gdsn(f, "sample.annot"))
    }
    res
}

#' makeVariantExperimentFromGDS
#' 
#' Conversion of snp gds file into SummarizedExperiment.
#' @param assayNames the gds node name that will be read into the
#'     \code{assays} slot and be represented as \code{DelayedArray}
#'     object.
#' @param rowDataColumns which columns of \code{rowData} to
#'     import. The default is NULL to read in all variant annotation
#'     info.
#' @param colDataColumns which columns of \code{colData} to
#'     import. The default is NULL to read in all sample related
#'     annotation info.
#' @param rowDataOnDisk whether to save the \code{rowData} as
#'     DelayedArray object. The default is TRUE.
#' @param colDataOnDisk whether to save the \code{colData} as
#'     DelayedArray object. The default is TRUE.
#' @return An \code{VariantExperiment} object.
#' @importFrom tools file_path_as_absolute
## #' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stats setNames
#' @export
#' @examples
#' file <- SNPRelate::snpgdsExampleFileName()
#' se <- makeVariantExperimentFromGDS(file)
#' rowData(se)
#' colData(se)
#' metadata(se)
#' ## Only read specific columns for feature annotation.
#' showAvailable(file)
#' se1 <- makeVariantExperimentFromGDS(file, rowDataColumns=c("ALLELE"))
#' rowRanges(se1)

makeVariantExperimentFromSNPGDS <- function(file, assayNames=NULL,
                                         rowDataColumns = NULL,
                                         colDataColumns = NULL,
                                         ## infoColumns = NULL,
                                         rowDataOnDisk = TRUE,
                                         colDataOnDisk = TRUE)
{
    ## checkings
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the gds file where the dataset is located."))
    file <- tools::file_path_as_absolute(file)
    stopifnot(is.character(assayNames) | is.null(assayNames))
    ## if (!isSingleStringOrNA(assayNames))
    ##     stop("'assayNames' must be a single string or NA")
    if(!isTRUEorFALSE(colDataOnDisk))
        stop("`colDataOnDisk` must be logical.")
    if(!isTRUEorFALSE(rowDataOnDisk))
        stop("`rowDataOnDisk` must be logical.")

    ## check which extensive gds format? SNPGDSFileClass or seqVarGDSClass? 
    ## ff <- .get_gds_fileFormat(file)

    allnames <- showAvailable(file)$assayNames
    if (is.null(assayNames)) {
        assayNames <- allnames
    } else {
        assayNames <- match.arg(assayNames, assayNames)
    }

    ## colData 
    colData <- .colData_snpgds(file, colDataColumns, colDataOnDisk)

    ## rowRange with info data if available for seqVarGDSClass
    rowRange <- .rowRanges_snpgds(file, rowDataColumns, rowDataOnDisk)

    ## if ((is.null(infoColumns) || length(infoColumns) > 0) && ff == "SEQ_ARRAY") {
    ##     infocols <- .info_seqgds(file, infoColumns, rowDataOnDisk)
    ##     mcols(rowRange) <- cbind(mcols(rowRange), infocols)
    ## }

    ## assay data
    ## Make sure all assays are in correct dimensions (feature*sample*else) 
    assays <- setNames(lapply(assayNames, function(x) GDSArray(file, x)), assayNames)
    ans_nrow <- length(rowRange)
    ans_ncol <- nrow(colData)
    permFun <- function(x, dim1, dim2) {
        pos <- match(c(dim1, dim2), dim(x))
        if (length(dim(x)) > 2) {
            aperm(x, perm = c(pos, setdiff(seq_along(dim(x)), pos)))
        } else {
            aperm(x, perm = pos)
        }
    }
    assays <- lapply(assays, permFun, dim1 = ans_nrow, dim2 = ans_ncol)
    
    se <- VariantExperiment(
        assays = assays,
        colData = colData,
        rowRanges = rowRange)
}
