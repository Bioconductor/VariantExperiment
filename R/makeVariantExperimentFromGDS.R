## makeVariantExperimentFromGDS requires that the input file has nodes for:
## 1. sample id (needed in generating colData)
## 2. feature id (needed in rowRanges, currently is matched by "id" for all feature related nodes)
## 3. chromosome, position (needed in rowRanges, currently are matched using "chrom" and "pos")
## 4. stop position? (reflen=1L used in .granges_generalgds)

## grep("chromosome", "position", "id") as
## GenomicRanges::Granges(seqnames=chromosome, ranges = Iranges(start
## = position, end = position + reflen-1L)) where define reflen = 1L
## for snps, but can be other values for other data type)

.granges_generalgds <- function(gdsfile, feature.num, reflen = 1L, ...){
    ftnodes <- .get_gds_annonodes(gdsfile, feature.num)
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
    feature.num <- .get_gdsnode_desp(gdsfile, ftnode, "dim")
    if (length(feature.num) > 1 | any(feature.num == 0L))
        stop("Wrong feature node name is provided!")
    rr <- .granges_generalgds(gdsfile, feature.num)
    rowDataColumns <- .rowDataColumns_check(gdsfile, ftnode, rowDataColumns)
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
            lapply(rowDataColumns, function(x) read.gdsn(index.gdsn(f, x))),
            rowDataColumns)
    }
    mcols(rr) <- resDF
    rr    
}

.colData_generalgds <- function(gdsfile, smpnode, colDataColumns, colDataOnDisk) {
    colDataColumns <- .colDataColumns_check(gdsfile, colDataColumns, smpnode)

    ## if no available colDataColumns are selected, i.e.,
    ## colDataColumns = character(0), return an empty
    ## (Delayed)DataFrame with sample id as rownames.

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

#' makeVariantExperimentFromGDS
#' 
#' Conversion of gds files into SummarizedExperiment object.
#' @rdname makeVariantExperimentFromGDS
#' @param file the GDS file name to be converted.
#' @param ftnode the node name for feature id (e.g., "variant.id",
#'     "snp.id", etc.).
#' @param smpnode the node name for sample id (e.g., "sample.id").
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
#' @param infoColumns which columns of \code{infoColumns} to import
#'     for "SEQ_ARRAY" ("SeqVarGDSClass" gds class). The default is
#'     NULL to read in all available info columns.
#' @return An \code{VariantExperiment} object.
#' @importFrom tools file_path_as_absolute
#' @importFrom stats setNames
#' @export

#' @examples
#'
#' ## gds file from DNA-seq data
#' 
#' seqfile <- SeqArray::seqExampleFileName(type="gds")
#' ve <- makeVariantExperimentFromGDS(seqfile)
#' ## all assay data
#' names(assays(ve))
#' showAvailable(seqfile)
#'
#' ## only read specific columns for feature / sample annotation. 
#'
#' assayNamess <- showAvailable(seqfile)$assayNames
#' rowdatacols <- showAvailable(seqfile)$rowDataColumns
#' coldatacols <- showAvailable(seqfile)$colDataColumns
#' infocols <- showAvailable(seqfile)$infoColumns
#' ve1 <- makeVariantExperimentFromGDS(
#' seqfile,
#' assayNames = assayNamess[2],
#' rowDataColumns = rowdatacols[1:3],
#' colDataColumns = coldatacols[1],
#' infoColumns = infocols[c(1,3,5,7)],
#' rowDataOnDisk = FALSE,
#' colDataOnDisk = FALSE)
#' assay(ve1)
#' 
#' ## the rowData(ve1) and colData(ve1) are now in DataFrame format 
#'
#' rowData(ve1)
#' colData(ve1)
#'
#' ## gds file from genotyping data
#' 
#' snpfile <- SNPRelate::snpgdsExampleFileName()
#' ve <- makeVariantExperimentFromGDS(snpfile)
#' rowData(ve)
#' colData(ve)
#' metadata(ve)
#'
#' ## Only read specific columns for feature annotation.
#'
#' showAvailable(snpfile)
#' ve1 <- makeVariantExperimentFromGDS(snpfile, rowDataColumns=c("snp.allele"))
#' rowRanges(ve1)
#'
#' ## use specific conversion functions for certain gds types
#'
#' veseq <- makeVariantExperimentFromSEQGDS(seqfile)
#' vesnp <- makeVariantExperimentFromSNPGDS(snpfile)

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
    if (!is.null(ff) && ff == "SEQ_ARRAY") {
        return(makeVariantExperimentFromSEQGDS(file,
                                               ftnode = "variant.id",
                                               smpnode = "sample.id",
                                               assayNames,
                                               rowDataColumns,
                                               colDataColumns,
                                               infoColumns, 
                                               rowDataOnDisk,
                                               colDataOnDisk))
    } else if (!is.null(ff) && ff == "SNP_ARRAY") {
        return(makeVariantExperimentFromSNPGDS(file,
                                               ftnode = "snp.id",
                                               smpnode = "sample.id",
                                               assayNames,
                                               rowDataColumns,
                                               colDataColumns,
                                               rowDataOnDisk,
                                               colDataOnDisk))
    }  

    ## ELSE: FOR GENERAL GDS FILES
    
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

    ## ans_nrow <- .get_gdsnode_desp(file, ftnode, "dim")
    ## ans_ncol <- .get_gdsnode_desp(file, smpnode, "dim")
    
    ## assays
    all_assays <- showAvailable(file, ftnode = ftnode, smpnode = smpnode)$assayNames
    if (is.null(assayNames)) {
        assayNames <- all_assays
    } else {
        assayNames <- match.arg(assayNames, assayNames)
    }
    assays <- setNames(lapply(assayNames, function(x) GDSArray(file, x)), assayNames)
    
    ## colData 
    colData <- .colData_generalgds(file, smpnode, colDataColumns, colDataOnDisk)
    
    ## rowRange with info data if available for seqVarGDSClass
    rowRange <- .rowRanges_generalgds(file, ftnode, rowDataColumns, rowDataOnDisk)

    ## assay data adjust dimensions into: feature*sample*else 
    ans_nrow <- length(rowRange) 
    ans_ncol <- nrow(colData) 
    assays <- lapply(assays, .permdim, dim1 = ans_nrow, dim2 = ans_ncol)
    
    se <- VariantExperiment(
        assays = assays,
        colData = colData,
        rowRanges = rowRange)
}


