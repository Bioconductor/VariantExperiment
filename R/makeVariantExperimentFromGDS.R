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

.varnodes <- function(gdsfile, fileFormat, name){
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
    node
}
.varnode_gdsdata_ondisk <- function(gdsfile, fileFormat, name){
    node <- .varnodes(gdsfile, fileFormat, name)
    GDSArray(gdsfile, node)
}

#' @importMethodsFrom SeqArray info
#' @import DelayedDataFrame
.infonodes <- function(seqArrayFile){
    f <- seqOpen(seqArrayFile)
    on.exit(seqClose(f))
    ls.gdsn(index.gdsn(f, "annotation/info"))
}
.infonodes_val <- function(seqArrayFile, name) {
    f <- seqOpen(seqArrayFile)
    on.exit(seqClose(f))
    SeqArray::info(f, info = name)
}

.info_seqgds <- function(seqArrayFile, infoColumns, rowDataOnDisk){
    infonodes <- .infonodes(seqArrayFile)
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
    ans_dim <- .get_gdsnode_desp(seqArrayFile, "variant.id", "dim")
    if(rowDataOnDisk){
        res <- lapply(infonodes, function(x) GDSArray(seqArrayFile, x))
        ## for "SEQ_ARRAY", dimensions don't have restrictions, so check and adjust. 
        dim_short <- lengths(res) != ans_dim
        if (any(dim_short)) {
            node <- infoColumns[dim_short]
            node_val <- .infonodes_val(seqArrayFile, node)[,1]
            res[[which(dim_short)]] <- DelayedArray(array(unlist(node_val)))
        }
        res1 <- DelayedDataFrame(lapply(res, I))
    }else{
        res1 <- .infonodes_val(seqArrayFile, infoColumns)
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

.sampnodes <- function(file, fileFormat){
    pre <- ifelse(fileFormat == "SNP_ARRAY", "sample.annot",
           ifelse(fileFormat == "SEQ_ARRAY", "sample.annotation", NULL))
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    ls.gdsn(index.gdsn(f, pre))
    
}
.sampnode_gdsdata_ondisk <- function(file, fileFormat, name)
{
    pre <- ifelse(fileFormat == "SNP_ARRAY", "sample.annot",
           ifelse(fileFormat == "SEQ_ARRAY", "sample.annotation", NULL))
    sampnodes <- .sampnodes(file, fileFormat)
    if(name %in% sampnodes){
        node <- paste0(pre, "/", name)
    }else{
        node <- name
    }
    GDSArray(file, node)  ## returns Delayed/GDSArray, not DF.
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

.get_gds_fileFormat <- function(file)
{
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    ff <- get.attr.gdsn(f$root)$FileFormat
    ff
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
    ff <- .get_gds_fileFormat(file)
    switch(ff, "SEQ_ARRAY" = .showAvailable_seqarray(file, args),
           "SNP_ARRAY" = .showAvailable_snparray(file, args))
}

.get_gdsnode_desp <- function(file, node, desp)
{
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    objdesp <- objdesp.gdsn(index.gdsn(f, node))
    desp <- match.arg(desp, names(objdesp))
    objdesp[[desp]]
}

## http://bioconductor.org/packages/release/bioc/vignettes/SeqArray/inst/doc/SeqArrayTutorial.html
## SeqArray requires 6 variables: sample.id, variant.id, position,
## chromosome, allele, genotype(folder).
.showAvailable_seqarray <- function(file,
                                    args = c("name",
                                             "rowDataColumns",
                                             "colDataColumns",
                                             "infoColumns")) {
    res <- CharacterList()
    if ("name" %in% args) {
        names <- gdsnodes(file)
        isarray <- vapply(names,
                          function(x) .get_gdsnode_desp(file, x, "is.array"),
                          logical(1))
        dims <- lapply(names, function(x) .get_gdsnode_desp(file, x, "dim"))
        res$name <- names[isarray & lengths(dims) > 1 & 
              ! vapply(dims, function(x) any(x == 0L), logical(1)) &
              !grepl("~", names)]
    }
    if ("rowDataColumns" %in% args) {
        res$rowDataColumns <- c("ID", "ALT", "REF", "QUAL", "FILTER")
    }
        if (any(c("colDataColumns", "infoColumns") %in% args)) {
        f <- openfn.gds(file)
        on.exit(closefn.gds(f))
    }
    if ("colDataColumns" %in% args)
        res$colDataColumns <- ls.gdsn(index.gdsn(f, "sample.annotation"))
    if("infoColumns" %in% args)
        res$infoColumns <- ls.gdsn(index.gdsn(f, "annotation/info"))
    res
}

## http://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
## SNPRelate requires 5 variables: sample.id, snp.id, snp.position, snp.chromosome, genotype
## 2 optional variables: snp.rs.id, snp.allele (when merging gentypes from diff platforms)
.showAvailable_snparray <- function(file,
                          args = c("name",
                                   "rowDataColumns",
                                   "colDataColumns")) {
    res <- CharacterList()
    if("name" %in% args)
        res$name <- "genotype"
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
#' ## se <- makeVariantExperimentFromGDS(file)
#' ## rowData(se)
#' ## colData(se)
#' ## metadata(se)
#' ## Only read specific columns for feature annotation.
#' showAvailable(file)
#' ## se1 <- makeVariantExperimentFromGDS(file, rowDataColumns=c("ALLELE"))
#' ## SummarizedExperiment::rowRanges(se1)

#' file <- SeqArray::seqExampleFileName(type="gds")
#' ## se <- makeVariantExperimentFromGDS(file)
#' ## all assay data
#' ## names(assays(se))
#' ## showAvailable(file)
#'
#' ## only read specific columns for feature / sample annotation. 
#' names <- showAvailable(file, "name")$name
#' rowdatacols <- showAvailable(file, "rowDataColumns")$rowDataColumns
#' coldatacols <- showAvailable(file, "colDataColumns")$colDataColumns
#' infocols <- showAvailable(file, "infoColumns")$infoColumns
#' ## se1 <- makeVariantExperimentFromGDS(
#' ## file,
#' ## name = names[2],
#' ## rowDataColumns = rowdatacols[1:3],
#' ## colDataColumns = coldatacols[1],
#' ## infoColumns = infocols[c(1,3,5,7)],
#' ## rowDataOnDisk = FALSE,
#' ## colDataOnDisk = FALSE)
#' ## assay(se1)
#' 
#' ## the rowData(se1) and colData(se1) are now in DataFrame format 
#' ## rowData(se1)
#' ## colData(se1)

#' @export
#' 
makeVariantExperimentFromGDS <- function(file, name=NULL,
                                         rowDataColumns = NULL,
                                         colDataColumns = NULL,
                                         infoColumns = NULL,
                                         rowDataOnDisk = TRUE,
                                         colDataOnDisk = TRUE)
{
    ## checkings
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
    ff <- .get_gds_fileFormat(file)
    allnames <- showAvailable(file)$name
    if (is.null(name)) {
        name <- allnames
    } else {
        name <- match.arg(name, names)
    }

    ## colData 
    colData <- .colData_gdsdata(file, ff, colDataColumns, colDataOnDisk)

    ## rowRange with info data if available for seqVarGDSClass
    rowRange <- .rowRanges_gdsdata(file, ff, rowDataColumns, rowDataOnDisk)
    if ((is.null(infoColumns) || length(infoColumns) > 0) && ff == "SEQ_ARRAY") {
        infocols <- .info_seqgds(file, infoColumns, rowDataOnDisk)
        mcols(rowRange) <- cbind(mcols(rowRange), infocols)
    }

    ## assay data
    ## Make sure all assays are in correct dimensions (feature*sample*else) 
    assays <- setNames(lapply(name, function(x) GDSArray(file, x)), name)
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
