## input is VCF file, we need to read it into VCFArray for assay data
## (>1d), and read the other annotation columns as 1-d VCFArray and
## put into DelayedDataFrame, and then construct an VariantExperiment
## object.

## We also need to make sure that the downstream analysis works on the
## VCFArray. Currently, the VariantExperiment-methods are inherited
## from the SeqArray which take the GDS as input. But we can make the
## VariantAnnotation methods works for the VariantExperiment object
## with "VCFArray + DelayedDataFrame" combination.

#' VCFtoVE
#' @description \code{VCFtoVE} is the function to convert a vcf file
#'     into \code{VariantExperiment} object. The genotype data will be
#'     written as \code{VCFArray} format, which is saved in the
#'     ‘assays’ slot. The annotation info for variants or samples will
#'     be written as \code{DelayedDataFrame} object, and saved in the
#'     \code{rowData} or \code{colData} slot.
#' @param vcf.fn vcf file name. Could be character string or
#'     \code{VcfFile} object or \code{RangedVcfStack} object.
#' @param out.dir The directory to save the vcf data, and the newly
#'     generated \code{VariantExperiment} object with array data in
#'     \code{VCFArray} format and annotation data in
#'     \code{elayedDataFrame} format. The default is a temporary
#'     folder.
#' @import VCFArray
#' @import VariantAnnotation GenomicFiles
#' @export
#' @examples
#' ## the vcf file
#' vcf <- SeqArray::seqExampleFileName("vcf")
#' ve <- VCFtoVE(vcf)
#' ve1 <- VCFtoVE(VcfFile(vcf))
#' all.equal(ve, ve1)
#' ve0 <- makeSummarizedExperimentFromVCF(vcf)
#' ## RangedVcfStack
#' extdata <- system.file(package = "GenomicFiles", "extdata")
#' files <- dir(extdata, pattern="^CEUtrio.*bgz$", full=TRUE)[1:2]
#' names(files) <- sub(".*_([0-9XY]+).*", "\\1", basename(files))
#' seqinfo <- as(readRDS(file.path(extdata, "seqinfo.rds")), "Seqinfo")
#' stack <- GenomicFiles::VcfStack(files, seqinfo)
#' gr <- as(GenomicFiles::seqinfo(stack)[rownames(stack)], "GRanges")
#' rgstack <- GenomicFiles::RangedVcfStack(stack, rowRanges = gr)
#' ve2 <- VCFtoVE(rgstack)
#'
#' 
VCFtoVE <- function(vcf.fn, out.dir = tempfile(), replace = FALSE,
                    ## header = NULL,
                    info.import = NULL, fmt.import = NULL,
                    sample.info = NULL, ignore.chr.prefix = "chr",
                    reference = NULL, start = 1L, count = -1L,
                    parallel = FALSE, verbose = FALSE)
{
    browser()
    ## "vcf.fn" could be: 1) character string (file name), 2) VcfFile
    ## object, 3) RangedVcfStack object?? 
    ### stopifnot(is.character(vcf.fn), length(vcf.fn)==1L)
    stopifnot(isSingleString(vcf.fn) | is(vcf.fn, "VcfFile") | is(vcf.fn, "RangedVcfStack"))
    if (!isSingleString(out.dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory where to save the \"VariantExperiment\" ",
                  "object (the directory will be created)"))
    ## read VCF file into VCFArray objects
    avail <- vcfFields(vcf.fn)

    ## 2d VCFArrays to be assay data
    assays <- SimpleList(vector("list", length(avail$geno)))
    for (i in seq_along(avail$geno)) {
        assays[[i]] <- VCFArray(vcf.fn, name = avail$geno[i], pfix = "geno")
    }
    names(assays) <- avail$geno  ## FIXME: add prefix here? 
    ## 1d VCFArrays to be rowData/colData in DelayedDataFrame

    varannots <- SimpleList(vector("list", length(c(avail$fixed, avail$info))))
    for (j in seq_along(avail$fixed)) {
        va <- VCFArray(vcf.fn, name = avail$fixed[j], pfix = "fixed")
        varannots[[j]] <- va
    }
    for (k in seq_along(avail$info)) {
        va <- VCFArray(vcf.fn, name = avail$info[k], pfix = "info")
        varannots[[k+length(avail$fixed)]] <- va
    }
    names(varannots) <- c(avail$fixed, paste0("info_", avail$info))
    rowData <- DelayedDataFrame(lapply(varannots, I))
    ## setNames(rowData, paste0("info_", infoColumns))
    
    ## rowRange
    param <- ScanVcfParam(fixed = NA, info = NA, geno = NA)
    if (is(vcf.fn, "RangedVcfStack")) {
        vcfWhich(param) <- rowRanges(vcf.fn)
        readvcf <- readVcfStack(vcf.fn, param = param)
    } else {
        readvcf <- readVcf(vcf.fn, genome = "hg19", param = param)
    }
    gr <- rowRanges(readvcf)

    mcols(gr) <- rowData

    ## colData
    sample.id <- samples(scanVcfHeader(vcf.fn))
    colData <- DelayedDataFrame(row.names = sample.id)
    
    VariantExperiment(assays = assays, rowRanges = gr, colData = colData)
}



