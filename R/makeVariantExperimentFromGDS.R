## makeVariantExperimentFromGDS requires that the input file has nodes for:
## 1. sample id (needed in generating colData)
## 2. feature id (needed in rowRanges, currently is matched by "id" for all feature related nodes)
## 3. chromosome, position (needed in rowRanges, currently are matched using "chrom" and "pos")
## 4. stop position? (reflen=1L used in .granges_generalgds)

.granges_generalgds <- function(gdsfile, len.ftr, reflen = 1L, ...){
    ftnodes <- .get_gds_annonodes(gdsfile, len.ftr)
    f <- openfn.gds(gdsfile)
    on.exit(closefn.gds(f))
    ## Here trying to match for "chromosome" and "position" related
    ## nodes. If no match, return error.
    vidnd <- index.gdsn(f, ftnodes[grep("id", ftnodes)[1]], silent = TRUE)
    chrnd <- index.gdsn(f, ftnodes[grep("chrom", ftnodes)[1]], silent = TRUE)
    posnd <- index.gdsn(f, ftnodes[grep("pos", ftnodes)[1]], silent = TRUE)
    if (!is.null(chrnd) & !is.null(posnd)) {
        chr <- read.gdsn(chrnd)
        pos <- read.gdsn(posnd)
        gr <- GenomicRanges::GRanges(seqnames=chr,
                                     ranges=IRanges(start=pos, end=pos+reflen-1L),
                                     ...)
        if (!is.null(vidnd)){
            vid <- read.gdsn(vidnd)
            names(gr) <- as.integer(vid)
        }
        gr
    } else {
        stop("Can not find the related gds nodes containing chromosome and position information!")
    }
}

.rowRanges_generalgds <- function(gdsfile, ftnode, rowDataColumns, rowDataOnDisk) {
    len.ftr <- .get_gdsnode_desp(gdsfile, ftnode, "dim")
    if (length(len.ft) > 1 | any(len.ft == 0L))
        stop("Wrong feature node name is provided!")
    rr <- .granges_generalgds(gdsfile, len.ftr)
    rowDataColumns <- .rowDataColumns_check(gdsfile, rowDataColumns)
    ## if no available rowDataColumns are selected, i.e.,
    ## rowDataColumns = character(0), return an empty (Delayed)DataFrame
    ## for mcols()
    if (is.character(rowDataColumns) && length(rowDataColumns) == 0) {
        resDF <- .empty_rowData_DF(gdsfile, ftnode, rowDataOnDisk)
        mcols(rr) <- resDF
        return(rr)
    }
    ## if there are valid rowDataColumns
    rowDataColumns <- rowDataColumns[rowDataColumns != ftnode]    
    if(rowDataOnDisk){
        res <- setNames(
            lapply(rowDataColumns, function(x) GDSArray(gdsfile, x)), 
            rowDataColumns)
        resDF <- DelayedDataFrame(lapply(res, I))
    }else{ ## rowDataOnDisk = FALSE...
        f <- openfn.gds(gdsfile)
        on.exit(closefn.gds(f))
        resDF <- setNames(
            lapply(rowDataColumns, function(x) read.gdsn(index.gdsn(f, rowDataColumns))),
            rowDataColumns)
    }
    mcols(rr) <- resDF
    rr    
}

.colData_generalgds <- function(gdsfile, smpnode, colDataColumns, colDataOnDisk) {
    colDataColumns <- .colDataColumns_check(gdsfile, colDataColumns)

    ## if no available colDataColumns are selected, i.e.,
    ## colDataColumns = character(0), return an empty
    ## (Delayed)DataFrame with "sample.id" as rownames.
    ## FIXME: here assuming "sample.id" exists!!! 

    if (is.character(colDataColumns) && length(colDataColumns) == 0) { ## character(0)
        .empty_colData_DF(gdsfile, smpnode, colDataOnDisk)
    } else {
        colDataColumns <- colDataColumns[colDataColumns != smpnode]
        ## if there are valid rowDataColumns
        if (colDataOnDisk) {
            sample.id <- GDSArray(gdsfile, smpnode)
            annot <- setNames(
                lapply(colDataColumns, function(x) GDSArray(gdsfile, x)), 
                colDataColumns)
            DelayedDataFrame(lapply(annot, I),
                             row.names=as.character(sample.id))
        } else {
            f <- openfn.gds(gdsfile)
            on.exit(closefn.gds(f))
            sample.id <- read.gdsn(index.gdsn(f, smpnode))
            annot <- lapply(colDataColumns, function(x) read.gdsn(index.gdsn(f, x)))
            names(annot) <- colDataColumns
            DataFrame(annot, row.names=sample.id)
        }
    }
}

makeVariantExperimentFromGDS <- function(file, ftnode, smpnode,
                                         assayNames=NULL,
                                         rowDataColumns = NULL,
                                         colDataColumns = NULL,
                                         rowDataOnDisk = TRUE,
                                         colDataOnDisk = TRUE,
                                         infoColumns = NULL ## only used when "SEQ_ARRAY"
                                         ){ 
    ## check which extensive gds format? SNPGDSFileClass or seqVarGDSClass? 
    ff <- .get_gds_fileFormat(file)
    if (ff == "SEQ_ARRAY") {
        return(makeVariantExperimentFromSEQGDS(file, assayNames,
                                               rowDataColumns,
                                               colDataColumns,
                                               infoColumns, 
                                               rowDataOnDisk,
                                               colDataOnDisk))
    } else if (ff == "SNP_ARRAY") {
        return(makeVariantExperimentFromSNPGDS(file, assayNames,
                                               rowDataColumns,
                                               colDataColumns,
                                               rowDataOnDisk,
                                               colDataOnDisk))
    }  

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

    ## assays
    all_assays <- showAvailable(file)$assayNames
    if (is.null(assayNames)) {
        assayNames <- all_assays
    } else {
        assayNames <- match.arg(assayNames, assayNames)
    }

    ## colData 
    colData <- .colData_generalgds(file, smpnode, colDataColumns, colDataOnDisk)
    
    ## rowRange with info data if available for seqVarGDSClass
    rowRange <- .rowRanges_generalgds(file, ftnode, rowDataColumns, rowDataOnDisk)

    ## if ((is.null(infoColumns) || length(infoColumns) > 0) && ff == "SEQ_ARRAY") {
    ##     infocols <- .info_seqgds(file, infoColumns, rowDataOnDisk)
    ##     mcols(rowRange) <- cbind(mcols(rowRange), infocols)
    ## }

    ## assay data adjust dimensions into: feature*sample*else 
    assays <- setNames(lapply(assayNames, function(x) GDSArray(file, x)), assayNames)
    ans_nrow <- length(rowRange) ## can also be len.ft
    ans_ncol <- nrow(colData)  ## can also be len.smp
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


