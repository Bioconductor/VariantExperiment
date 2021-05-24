## In "maAkeVariantExperimentFromSEQGDS()", we have customized:

## 1) the "ID", "ALT", "REF", "QUAL", "FILTER" columns (from
## "annotation/id", "allele", "annotation/qual", "annotaiton/filter")
## nodes;
## 2) the info columns (from "annotation/info/*" nodes) as part of the
## mcols(rowRanges(ve)) with column names with "info_" prefix.
## "variant.id" was used as default feature id, and "sample.id" was
## used as default sample id. These 2 nodes are required by
## "SeqVarGDSClass" in SeqArdray.

#' @rawNamespace import(SeqArray, except = c(colData, rowRanges))

.granges_seqgds <- function(seqgdsfile, ...){
    f <- seqOpen(seqgdsfile)
    on.exit(seqClose(f))
    SeqArray::granges(f)
}

## return a DataFrame
.varnode_seqgds_inmem <- function(seqgdsfile, node){
    f <- seqOpen(seqgdsfile)
    on.exit(seqClose(f))
    res <- switch(node,
                  'annotation/id' = list(seqGetData(f, "annotation/id")),
                  'annotation/qual' = list(qual(f)),
                  'annotation/filter' = list(filt(f)),
                  allele = list(ref(f), alt(f)))
    if (node == "allele")
        node <- c("REF", "ALT")
    res <- setNames(res, node)
    resDF <- DataFrame(res)
    resDF  ## returns a DataFrame with names. 
}

#' @importFrom methods new
#' @importMethodsFrom SeqArray info
#' @import DelayedDataFrame

.infonodes_val <- function(seqgdsfile, name) {
    f <- seqOpen(seqgdsfile)
    on.exit(seqClose(f))
    SeqArray::info(f, info = name)
}

.infoColumns_check <- function(seqgdsfile, infoColumns) {
    infonodes <- showAvailable(seqgdsfile, "infoColumns")[[1]]
    if (is.null(infoColumns)) {
        infoColumns <- infonodes
    } else {
        idx <- infoColumns %in% infonodes
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
.infoColumns_seqgds <- function(seqgdsfile, ftnode, infoColumns, rowDataOnDisk){
    infoColumns <- .infoColumns_check(seqgdsfile, infoColumns)

    if(is.character(infoColumns) && length(infoColumns) == 0) {
        resDF <- .empty_rowData_DF(seqgdsfile, ftnode, rowDataOnDisk)
        return(resDF)
    }
    
    infonodes <- paste0("annotation/info/", infoColumns)
    ans_dim <- .get_gdsnode_desp(seqgdsfile, ftnode, "dim")
    if(rowDataOnDisk){
        res <- lapply(infonodes, function(x) GDSArray(seqgdsfile, x))
        ## only keep nodes that are same length as "variant.id"!!  for
        ## SeqArray::seqExampleFileName("gds"), the
        ## "annotation/info/AA" doesn't match (1328 vs 1348).
        idx.keep <- lengths(res) == ans_dim
        res <- res[idx.keep] 
        infoColumns <- infoColumns[idx.keep]
            res1 <- DelayedDataFrame(lapply(res, I))
    }else{
        res1 <- .infonodes_val(seqgdsfile, infoColumns)
    }
    setNames(res1, paste0("info.", infoColumns))    
}

#' @importFrom Biostrings DNAStringSet
#' @import SNPRelate
#' @importMethodsFrom DelayedArray sub
#' @importFrom S4Vectors mcols<- 
.rowRanges_seqgds <- function(seqgdsfile, ftnode = "variant.id", rowDataColumns, rowDataOnDisk){
    rr <- .granges_seqgds(seqgdsfile)
    ## following code generates the mcols(SummarizedExperiment::rowRanges(se))
    rowDataColumns <- .rowDataColumns_check(seqgdsfile, ftnode, rowDataColumns)

    ## if no available rowDataColumns are selected, i.e.,
    ## rowDataColumns = character(0), return an empty (Delayed)DataFrame
    ## for mcols()
    if (is.character(rowDataColumns) && length(rowDataColumns) == 0) { ## character(0)
        resDF <- .empty_rowData_DF(seqgdsfile, ftnode, rowDataOnDisk)
        mcol(rr) <- resDF
        return(rr)
    }
    ## if there are valid rowDataColumns
    if(rowDataOnDisk){
        res <- setNames(
            lapply(rowDataColumns, function(x) GDSArray(seqgdsfile, x)), 
            rowDataColumns)
        resDF <- DelayedDataFrame(lapply(res, I))
        if("allele" %in% rowDataColumns){ ## for seqgds
            resDF$REF <- sub(",.*", "", resDF$allele)
            resDF$ALT <- sub("[TCGA]*,", "", resDF$allele)
            resDF$allele <- NULL
            ## resDF <- resDF[, c("REF", "ALT", names(resDF)[!names(resDF) %in% c("REF", "ALT", "allele")])]
        }
    }else{ ## rowDataOnDisk = FALSE...
        resDF <- DataFrame(lapply(rowDataColumns, function(x)
            .varnode_seqgds_inmem(seqgdsfile, x)))
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

.colData_seqgds <- function(seqgdsfile, smpnode, colDataColumns, colDataOnDisk) {
    colDataColumns <- .colDataColumns_check(seqgdsfile, colDataColumns, smpnode)

    ## if no available colDataColumns are selected, i.e.,
    ## colDataColumns = character(0), return an empty (Delayed)DataFrame    
    if (is.character(colDataColumns) && length(colDataColumns) == 0) { ## character(0)
        .empty_colData_DF(seqgdsfile, smpnode, colDataOnDisk)
    } else {
        ## if there are valid rowDataColumns
        if (colDataOnDisk) {
            sample.id <- .sampnode_seqgds_ondisk(seqgdsfile, smpnode)
            annot <- setNames(
                lapply(colDataColumns, function(x)
                    .sampnode_seqgds_ondisk(seqgdsfile, x)),
                colDataColumns)
            DelayedDataFrame(lapply(annot, I),
                             row.names=as.character(sample.id))
        } else {
            f <- openfn.gds(seqgdsfile)
            on.exit(closefn.gds(f))
            sample.id <- read.gdsn(index.gdsn(f, smpnode))
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
                                             "infoColumns"),
                                    ftnode = "variant.id",
                                    smpnode = "sample.id") {
    res <- CharacterList()
    allnodes <- gdsnodes(file)
    if ("assayNames" %in% args) {
        res$assayNames <- .get_gds_arraynodes(file)
    }
    if ("rowDataColumns" %in% args) {
        ## "allele" is a required node in "SEQ_ARRAY";
        ## "annotation/id", "annotation/qual", "annotation/filter"are
        ## optional nodes in "SEQ_ARRAY".
        rowcols <- c("allele",
                     allnodes[allnodes %in% paste("annotation", c("id", "qual", "filter"), sep="/")])
        rowlens <- vapply(rowcols, function(x) .get_gdsnode_desp(file, x, "dim"), integer(1))
        rowcols <- rowcols[rowlens == .get_gdsnode_desp(file, ftnode, "dim")]
        res$rowDataColumns <- rowcols
        ## res$rowDataColumns <- c("ID", "ALT", "REF", "QUAL", "FILTER")
    }
    if ("colDataColumns" %in% args) {
        res$colDataColumns <- .showAvailable_check(file, "colDataColumns", ftnode, smpnode)
    }
    if("infoColumns" %in% args) {
        res$infoColumns <- .showAvailable_check(file, "infoColumns", ftnode, smpnode)
    }
    res
}

#' @rdname makeVariantExperimentFromGDS
#' @export

makeVariantExperimentFromSEQGDS <- function(file, ftnode = "variant.id",
                                            smpnode = "sample.id",
                                            assayNames=NULL,
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
    if(!isTRUEorFALSE(colDataOnDisk))
        stop("`colDataOnDisk` must be logical.")
    if(!isTRUEorFALSE(rowDataOnDisk))
        stop("`rowDataOnDisk` must be logical.")
    
    ## assay data
    allnames <- showAvailable(file, "assayNames")[[1]]
    if (is.null(assayNames)) {
        assayNames <- allnames
    } else {
        assayNames <- match.arg(assayNames, assayNames)
    }
    assays <- setNames(lapply(assayNames, function(x) GDSArray(file, x)), assayNames)
    
    ## colData 
    colData <- .colData_seqgds(file, smpnode, colDataColumns, colDataOnDisk)
    
    ## rowRange with info data if available for seqVarGDSClass
    rowRange <- .rowRanges_seqgds(file, ftnode, rowDataColumns, rowDataOnDisk)
    if ((is.null(infoColumns) || length(infoColumns) > 0)) {
        infocols <- .infoColumns_seqgds(file, ftnode, infoColumns, rowDataOnDisk)
        mcols(rowRange) <- cbind(mcols(rowRange), infocols)
    }

    ## Make sure all assays are in correct dimensions (feature*sample*else) 
    ans_nrow <- length(rowRange)
    ans_ncol <- nrow(colData)
    assays <- lapply(assays, .permdim, dim1 = ans_nrow, dim2 = ans_ncol)
    
    se <- VariantExperiment(
        assays = assays,
        colData = colData,
        rowRanges = rowRange)
}
