## In "makeVariantExperimentFromSEQGDS()", we have customized the info
## columns as part of the mcols(rowRanges(ve)). "variant.id" was used
## as default sample id, and "sample.id" was used as default feature
## id. These 2 nodes are required by "SeqVarGDSClass" in SeqArray. 


#' @rawNamespace import(SeqArray, except = c(colData, rowRanges))

.granges_seqgds <- function(seqgdsfile, ...){
    f <- seqOpen(seqgdsfile)
    on.exit(seqClose(f))
    SeqArray::granges(f)
}

## return a DataFrame
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

.varnode_seqgds <- function(seqgdsfile, name){
    f <- openfn.gds(seqgdsfile)
    on.exit(closefn.gds(f))
    varid <- read.gdsn(index.gdsn(f, "variant.id"))
    varnodes <- ls.gdsn(index.gdsn(f, "annotation"))
    if(name %in% varnodes) node <- paste0("annotation/", name)
    if(name %in% c("alt", "ref")) node <- "allele"
    node
}

.varnode_seqgds_ondisk <- function(seqgdsfile, name){
    node <- .varnode_seqgds(seqgdsfile, name)
    GDSArray(seqgdsfile, node)
}

#' @importMethodsFrom SeqArray info
#' @import DelayedDataFrame
.infonodes <- function(seqgdsfile){
    f <- seqOpen(seqgdsfile)
    on.exit(seqClose(f))
    ls.gdsn(index.gdsn(f, "annotation/info"))
}
.infonodes_val <- function(seqgdsfile, name) {
    f <- seqOpen(seqgdsfile)
    on.exit(seqClose(f))
    SeqArray::info(f, info = name)
}

.infoColumns_check <- function(seqgdsfile, infoColumns) {
    infonodes <- .infonodes(seqgdsfile)
    if (is.null(infoColumns)) {
        infoColumns <- infonodes
    } else {
        idx <- toupper(infoColumns) %in% infonodes
        if(any(!idx)){
            misics <- paste(infoColumns[!idx], collapse = ", ")
            warning(.infoColumns_check_msg(misics))
            infoColumns <- toupper(infoColumns[idx])
        }       
    }
    infoColumns
}

.infoColumns_check_msg <- function(missingICs) {
    paste("\n", 'The "infoColumns" argument of "',
          missingICs, 
          '" does not exist!', "\n",
          'Please use showAvailable(file, "infoColumns") ',
          'to get the available columns for "infoColumns."', "\n")
}

## return a DataFrame with names.
.infoColumns_seqgds <- function(seqgdsfile, infoColumns, rowDataOnDisk){
    infoColumns <- .infoColumns_check(seqgdsfile, infoColumns)
    infonodes <- paste0("annotation/info/", infoColumns)

    ## for "SEQ_ARRAY", dimensions don't have restrictions, so check and adjust.
    ans_dim <- .get_gdsnode_desp(seqgdsfile, "variant.id", "dim")
    if(rowDataOnDisk){
        res <- lapply(infonodes, function(x) GDSArray(seqgdsfile, x))
        dim_short <- lengths(res) != ans_dim
        if (any(dim_short)) {
            node <- infoColumns[dim_short]
            node_val <- .infonodes_val(seqgdsfile, node)[,1]
            res[[which(dim_short)]] <- DelayedArray(array(unlist(node_val)))
        }
        res1 <- DelayedDataFrame(lapply(res, I))
    }else{
        res1 <- .infonodes_val(seqgdsfile, infoColumns)
    }
    setNames(res1, paste0("info_", infoColumns))

    }

#' @importFrom Biostrings DNAStringSet
#' @import SNPRelate
#' @importMethodsFrom DelayedArray sub
.rowRanges_seqgds <- function(seqgdsfile, rowDataColumns, rowDataOnDisk){
    rr <- .granges_seqgds(seqgdsfile)
    ## following code generates the mcols(SummarizedExperiment::rowRanges(se))
    rowDataColumns <- .rowDataColumns_check(seqgdsfile, rowDataColumns)

    ## if no available rowDataColumns are selected, i.e.,
    ## rowDataColumns = character(0), return an empty (Delayed)DataFrame
    ## for mcols()
    if (is.character(rowDataColumns) && length(rowDataColumns) == 0) { ## character(0)
        resDF <- .empty_rowData_DF(seqgdsfile, "variant.id", rowDataOnDisk)
        mcol(rr) <- resDF
        return(rr)
    }
    ## if there are valid rowDataColumns
    if(rowDataOnDisk){
        res <- setNames(
            lapply(rowDataColumns, function(x)
                .varnode_seqgds_ondisk(seqgdsfile, name=x)),
            toupper(rowDataColumns))
        resDF <- DelayedDataFrame(lapply(res, I))
        if("REF" %in% names(resDF)){ ## for seqgds
            resDF$REF <- sub(",.*", "", resDF$REF)
        }
        if("ALT" %in% names(resDF)){ ## for seqgds
            resDF$ALT <- sub("[TCGA]*,", "", resDF$ALT)
        }
    }else{ ## rowDataOnDisk = FALSE...
        resDF <- DataFrame(lapply(rowDataColumns, function(x)
            .varnode_seqgds_inmem(file, x)))
    }
    mcols(rr) <- resDF
    rr
}

.sampnodes_seqgds <- function(seqgdsfile) {
    f <- openfn.gds(seqgdsfile)
    on.exit(closefn.gds(f))
    ls.gdsn(index.gdsn(f, "sample.annotation"))
}

## returns Delayed/GDSArray, not DF.
.sampnode_seqgds_ondisk <- function(seqgdsfile, name) {
    sampnodes <- .sampnodes_seqgds(seqgdsfile)
    if(name %in% sampnodes){
        node <- paste0("sample.annotation/", name)
    }else{
        node <- name
    }
    GDSArray(seqgdsfile, node)  
}

###
## colData for samples
###

.colData_seqgds <- function(seqgdsfile, colDataColumns, colDataOnDisk) {
    colDataColumns <- .colDataColumns_check(seqgdsfile, colDataColumns)

    ## if no available colDataColumns are selected, i.e.,
    ## colDataColumns = character(0), return an empty (Delayed)DataFrame    
    if (is.character(colDataColumns) && length(colDataColumns) == 0) { ## character(0)
        .empty_colData_DF(seqgdsfile, colDataOnDisk)
    } else {
        ## if there are valid rowDataColumns
        if (colDataOnDisk) {
            sample.id <- .sampnode_seqgds_ondisk(seqgdsfile, "sample.id")
            annot <- setNames(
                lapply(colDataColumns, function(x)
                    .sampnode_seqgds_ondisk(seqgdsfile, x)),
                colDataColumns)
            DelayedDataFrame(lapply(annot, I),
                             row.names=as.character(sample.id))
        } else {
            f <- openfn.gds(seqgdsfile)
            on.exit(closefn.gds(f))
            sample.id <- read.gdsn(index.gdsn(f, "sample.id"))
            node <- paste0("sample.annotation/", colDataColumns)
            annot <- lapply(node, function(x) read.gdsn(index.gdsn(f, x)))
            names(annot) <- colDataColumns
            DataFrame(annot, row.names=sample.id)
        }
    }
}

## http://bioconductor.org/packages/release/bioc/vignettes/SeqArray/inst/doc/SeqArrayTutorial.html
## SeqArray requires 6 variables: sample.id, variant.id, position,
## chromosome, allele, genotype(folder).
.showAvailable_seqarray <- function(file,
                                    args = c("assayNames",
                                             "rowDataColumns",
                                             "colDataColumns",
                                             "infoColumns")) {
    res <- CharacterList()
    if ("assayNames" %in% args) {
        res$assayNames <- .get_gds_arraynodes(file)
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

#' makeVariantExperimentFromGDS
#' 
#' Conversion of seq gds file into SummarizedExperiment.
#' @inheritParams makeVariantExperimentFromSNPGDS 
#' @param infoColumns which columns of \code{infoColumns} to
#' @return An \code{VariantExperiment} object.
#' @export
#' @examples
#' file <- SeqArray::seqExampleFileName(type="gds")
#' se <- makeVariantExperimentFromGDS(file)
#' ## all assay data
#' names(assays(se))
#' showAvailable(file)
#'
#' ## only read specific columns for feature / sample annotation. 
#' assayNamess <- showAvailable(file)$assayNames
#' rowdatacols <- showAvailable(file)$rowDataColumns
#' coldatacols <- showAvailable(file)$colDataColumns
#' infocols <- showAvailable(file)$infoColumns
#' se1 <- makeVariantExperimentFromGDS(
#' file,
#' assayNames = assayNamess[2],
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

makeVariantExperimentFromSEQGDS <- function(file, assayNames=NULL,
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
    colData <- .colData_seqgds(file, colDataColumns, colDataOnDisk)

    ## rowRange with info data if available for seqVarGDSClass
    rowRange <- .rowRanges_seqgds(file, rowDataColumns, rowDataOnDisk)
    ## if ((is.null(infoColumns) || length(infoColumns) > 0) && ff == "SEQ_ARRAY") {
    if ((is.null(infoColumns) || length(infoColumns) > 0)) {
        infocols <- .infoColumns_seqgds(file, infoColumns, rowDataOnDisk)
        mcols(rowRange) <- cbind(mcols(rowRange), infocols)
    }

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


## todo: makeVariantExperimentFromSEQGDS(), makeVariantExperimentFromSNPGDS()
