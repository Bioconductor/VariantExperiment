## VCF2VE function
## VCF2VE <- function(vcf.fn, out.fn, header=NULL,
##                    storage.option="LZMA_RA", info.import=NULL, fmt.import=NULL,
##                    genotype.var.name="GT", ignore.chr.prefix="chr",
##                    reference=NULL, start=1L, count=-1L, optimize=TRUE, raise.error=TRUE,
##                    digest=TRUE, parallel=FALSE, verbose=TRUE)

#' @name VCF2VE
#' @rdname VariantExperiment-class.Rd
#' @description \code{VCF2VE} is the function to convert a vcf file
#'     into \code{VariantExperiment} object. By default, the genotype
#'     data will be written as \code{GDSArray} format, which is saved
#'     in the \code{assays} slot. The annotation info for variants
#'     will be written as \code{DelayedDataFrame} object (by default)
#'     or \code{DataFrame} object, and saved in the \code{rowData}
#'     slot. The annotation info for samples will be written as
#'     \code{DelayedDataFrame} object (by default) or \code{DataFrame}
#'     object, and saved in the \code{colData} slot. 
#' @param replace Whether to replace the directory if it already
#'     exists. The default is FALSE.
#' @param out.dir The directory to save the gds format of the vcf
#'     data, and the newly generated VariantExperiment object with
#'     array data in \code{GDSArray} format and annotation data in
#'     \code{DelayedDataFrame} format.
#' @param replace Whether to replace the directory if it already
#'     exists. The default is FALSE.
#' @param head if NULL, ‘header’ is set to be ‘seqVCF_Header(vcf.fn)’,
#'     which is a list (with a class name "SeqVCFHeaderClass", S3
#'     object).
#' @param compress the compression method for writing the gds
#'     file. The default is "LZMA_RA". See ‘?SeqArray::seqVCF2GDS’ for
#'     more details of this argument.
#' @param annotationOnDisk whether to save the annotation info for
#'     samples and variants as Delayed object. The default is TRUE.
#' @param ignore.chr.prefix a vector of character, indicating the
#'     prefix of chromosome which should be ignored, like "chr"; it is
#'     not case-sensitive.
#' @param reference genome reference, like "hg19", "GRCh37"; if the
#'     genome reference is not available in VCF files, users could
#'     specify the reference here.
#' @param start the starting variant if importing part of VCF files.
#' @param count the maximum count of variant if importing part of VCF
#'     files, -1 indicates importing to the end.
#' @param verbose whether to print the process messages. The default
#'     is FALSE.
#' @export
VCF2VE <- function(vcf.fn, out.dir="my_gds_se", replace=FALSE, header=NULL,
                   info.import=NULL, fmt.import=NULL,
                   annotationOnDisk = TRUE,
                   ignore.chr.prefix="chr",
                   reference=NULL, start=1L, count=-1L,
                   parallel=FALSE, verbose=TRUE){
    
    ## check
    stopifnot(is.character(vcf.fn), length(vcf.fn)==1L)
    if (!isSingleString(out.dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory where to save the ", "VariantExperiment",
                  " object (the directory will be created)"))
    if (!isTRUEorFALSE(replace))
        stop("'replace' must be TRUE or FALSE")
    
    ## stopifnot(is.character(out.dir), length(out.dir)==1L)

    ## run VCF to GDS
    .create_dir(out.dir, replace)
    out.gds.fn <- file.path(out.dir, "se.gds")
    seqVCF2GDS(vcf.fn, out.gds.fn, header = header,
               info.import = info.import, fmt.import = fmt.import,
               ignore.chr.prefix = ignore.chr.prefix,
               reference = reference, optimize=TRUE, raise.error=TRUE,
               verbose=verbose)

    ## run GDS to VE
    makeSummarizedExperimentFromGDS(
        file=out.gds.fn, name=NULL,
        ## rowDataColumns = rowDataColumns,
        ## colDataColumns = colDataColumns,
        infoColumns = info.import,  ## ??
        rowDataOnDisk = annotationOnDisk,
        colDataOnDisk = annotationOnDisk)
}

.saveGDSMaybe <- function(gdsfile) {
    file <- gdsfile(gdsfile)[1]
    gdsdim <- c(seqSummary(file, verbose=FALSE)$num.variant,
                seqSummary(file, verbose=FALSE)$num.sample)
    sedim <- dim(gdsfile) ## variants * samples
    if(!identical(sedim, gdsdim))
        stop("use 'saveGDSSummarizedExperiment()' to synchronize on-disk and in-memory representations")
        ## tmpdir <- tempfile()
        ## dir.create(tmpdir)
        ## gdsfile <- saveGDSSummarizedExperiment(gdsfile, dir=tmpdir, replace=TRUE)
        ## warning("updated 'gdsfile' to temporary on-disk location")
    gdsfile
}

.doCompatibleFunction <- function(gdsfile, ..., FUN) {
    gdsfile <- .saveGDSMaybe(gdsfile)
    file <- gdsfile(gdsfile)[1]
    f <- SeqArray::seqOpen(file)
    on.exit(SeqArray::seqClose(f))
    FUN(f, ...)
}

## SeqArray stats: seqAlleleFreq seqAlleleCount seqMissing seqNumAllele 

## #' @importFrom SeqArray seqAlleleFreq seqAlleleCount seqMissing seqNumAllele
#' @import SeqArray 
.seqAlleleFreq <- function(gdsfile, ref.allele=0L, .progress=FALSE, parallel = seqGetParallel()){
    alleleFreq <- .doCompatibleFunction(
        gdsfile, ref.allele=ref.allele, .progress=.progress, parallel=parallel,
        FUN=SeqArray::seqAlleleFreq
    )
    alleleFreq
    ## id <- paste(ref.allele, collapse="_")
    ## rowData(gdsfile)[[paste0("seqAlleleFreq_", id)]] <- alleleFreq
    ## gdsfile
}
setMethod("seqAlleleFreq", "VariantExperiment", .seqAlleleFreq)        

.seqAlleleCount <- function(gdsfile, ref.allele=0L, .progress=FALSE, parallel = seqGetParallel()){
    alleleCount <- .doCompatibleFunction(
        gdsfile, ref.allele=ref.allele, .progress=.progress, parallel=parallel,
        FUN=SeqArray::seqAlleleCount
    )
    alleleCount
    ## id <- paste(ref.allele, collapse="_")
    ## rowData(gdsfile)[[paste0("seqAlleleCount_", id)]] <- alleleCount
    ## gdsfile
}

setMethod("seqAlleleCount", "VariantExperiment", .seqAlleleCount)        

.seqMissing <- function(gdsfile, per.variant=TRUE, .progress=FALSE, parallel = seqGetParallel()){
    seqMissing <- .doCompatibleFunction(
        gdsfile, per.variant=per.variant, .progress =.progress, parallel=parallel,
        FUN=SeqArray::seqMissing
    )
    seqMissing
    ## if (per.variant)
    ##     rowData(gdsfile)[["variant_missing"]] <- seqMissing
    ## else
    ##     VariantExperiment::colData(gdsfile)[["sample_missing"]] <- seqMissing
    ## gdsfile
}
setMethod("seqMissing", "VariantExperiment", .seqMissing)        

.seqNumAllele <- function(gdsfile){
    numAllele <- .doCompatibleFunction(gdsfile, FUN=SeqArray::seqNumAllele)
    numAllele
    ## rowData(gdsfile)[["numAllele"]] <- numAllele
    ## gdsfile
}
setMethod("seqNumAllele", "VariantExperiment", .seqNumAllele)        

## SeqVarTools stats functions:

## already works: alleleFrequency(SeqArray::seqAlleleFreq), duplicateDiscordance, granges(SummarizedExperiment::granges), nAlleles(SeqArray::seqNumAllele), nAllele(SeqArray::seqAlleleCount)

## implemented: hwe, inbreedCoeff, pca, titv, refDosage, altDosage, countSingletons, homozygosity, heterozygosity, meanBySample, missingGenotypeRate, isSNV, isVariant

## remove: getGenotype/getGenotypeAlleles/expandedAltDosage/alleleDosage(SeqVarGDSClass,numeric)/alleleDosage(SeqVarGDSClass,list)

## Iterators(Extends ‘SeqVarData’ to provide iterators over variants.?)
## variantRanges(x): Get the variant ranges.
## variantFilter(x): Get the list of variant indices. 
## currentRanges(x): Get the ranges selected in the current iteration.
## duplicateDiscordance (SeqVarData),
## granges(already works for SE),

## Methods for class "SeqVarData" in package "SeqVarTools":
## SeqVarData(gds, sample.data) ## "sample.data" in class "AnnotatedDataFrame". 
## 1. alternateAlleleDetection
## 2. alleleFrequency
## 3. duplicateDiscordance
## 4. regression ?
## 5. sampleData, sampleData<-
## 6. validateSex
## 7. variantData, variantData<-

## #' @importFrom SeqVarTools hwe inbreedCoeff pca refDosage altDosage countSingletons homozygosity heterozygosity meanBySample missingGenotypeRate
## titv, isSNV, isVariant
.hwe <- function(gdsobj, permute=FALSE){
    hwe <- .doCompatibleFunction(gdsobj, permute=permute, FUN = SeqVarTools::hwe)
    hwe   ## returns a data.frame, will output as is, do not paste into rowData(se)
}
setMethod("hwe", "VariantExperiment", .hwe)

.inbreedCoeff <- function(gdsobj, margin=c("by.variant", "by.sample"),
                          use.names=FALSE){
    inbCoef <- .doCompatibleFunction(
        gdsobj, margin=margin, use.names=use.names, FUN = SeqVarTools::inbreedCoeff)
    inbCoef   ## returns a named (if use.names=TRUE) vector
}
setMethod("inbreedCoeff", "VariantExperiment", .inbreedCoeff)

.pca <- function(gdsobj, eigen.cnt=32){
    pca <- .doCompatibleFunction(gdsobj, eigen.cnt=eigen.cnt, FUN = SeqVarTools::pca)
    pca   ## returns a list, $eigenval (vector), $eigenvect (matrix)
}
setMethod("pca", "VariantExperiment", .pca)

.titv <- function(gdsobj, by.sample=FALSE, use.names=FALSE){
    titv <- .doCompatibleFunction(gdsobj, by.sample=by.sample,
                                  use.names=use.names, FUN = SeqVarTools::titv)
    titv   ## returns a scalar / vector (if by.sample=TRUE)
}
setMethod("titv", "VariantExperiment", .titv)

.refDosage <- function(gdsobj, use.names=TRUE){
    dos <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                  FUN = SeqVarTools::refDosage)
    dos   ## returns a matrix, nsamp * nvar
}
setMethod("refDosage", "VariantExperiment", .refDosage)

.altDosage <- function(gdsobj, use.names=TRUE){
    dos <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                  FUN = SeqVarTools::altDosage)
    dos   ## returns a matrix, nsamp * nvar
}
setMethod("altDosage", "VariantExperiment", .altDosage)

.ctSingleton <- function(gdsobj, use.names=FALSE){
    ct <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                  FUN = SeqVarTools::countSingletons)
    ct   ## returns a vector of the number of singleton variants per sample.
}
setMethod("countSingletons", "VariantExperiment", .ctSingleton)

.heterozygosity <- function(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE){
    margin <- match.arg(margin)
    hetero <- .doCompatibleFunction(gdsobj, margin=margin, use.names=use.names,
                                  FUN = SeqVarTools::heterozygosity)
    hetero   ## returns a vector of the number of singleton variants per sample.
}
setMethod("heterozygosity", "VariantExperiment", .heterozygosity)

.homozygosity <- function(gdsobj, allele=c("any", "ref", "alt"), margin=c("by.variant", "by.sample"), use.names=FALSE){
    allele <- match.arg(allele)
    margin <- match.arg(margin)
    homo <- .doCompatibleFunction(gdsobj, allele=allele, margin=margin,
                                  use.names=use.names,
                                  FUN = SeqVarTools::homozygosity)
    homo   ## returns a vector of the number of singleton variants per sample.
}
setMethod("homozygosity", "VariantExperiment", .homozygosity)

.meanBySample <- function(gdsobj, var.name, use.names=FALSE){
    mean <- .doCompatibleFunction(gdsobj, var.name=var.name, use.names=use.names,
                                  FUN = SeqVarTools::meanBySample)
    mean
} 
setMethod("meanBySample", "VariantExperiment", .meanBySample)

.missingGenotypeRate <- function(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE){
    margin <- match.arg(margin)
    misRate <- .doCompatibleFunction(gdsobj, margin=margin, use.names=use.names,
                                     FUN = SeqVarTools::missingGenotypeRate)
    misRate
}
setMethod("missingGenotypeRate", "VariantExperiment", .missingGenotypeRate)

.isSNV <- function(gdsobj, biallelic=TRUE){
    issnv <- .doCompatibleFunction(gdsobj, biallelic=biallelic,
                                   FUN = SeqVarTools::isSNV)
    issnv
}
setMethod("isSNV", "VariantExperiment", .isSNV)

.isVariant <- function(gdsobj, use.names=FALSE){
    isvar <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                   FUN = SeqVarTools::isVariant)
    isvar
}
setMethod("isVariant", "VariantExperiment", .isVariant)
