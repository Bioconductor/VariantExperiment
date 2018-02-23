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

    id <- paste(ref.allele, collapse="_")
    rowData(gdsfile)[[paste0("seqAlleleFreq_", id)]] <- alleleFreq
    gdsfile
}
setMethod("seqAlleleFreq", "SummarizedExperiment", .seqAlleleFreq)        

.seqAlleleCount <- function(gdsfile, ref.allele=0L, .progress=FALSE, parallel = seqGetParallel()){
    alleleCount <- .doCompatibleFunction(
        gdsfile, ref.allele=ref.allele, .progress=.progress, parallel=parallel,
        FUN=SeqArray::seqAlleleCount
    )

    id <- paste(ref.allele, collapse="_")
    rowData(gdsfile)[[paste0("seqAlleleCount_", id)]] <- alleleCount
    gdsfile
}
setMethod("seqAlleleCount", "SummarizedExperiment", .seqAlleleCount)        

.seqMissing <- function(gdsfile, per.variant=TRUE, .progress=FALSE, parallel = seqGetParallel()){
    seqMissing <- .doCompatibleFunction(
        gdsfile, per.variant=per.variant, .progress =.progress, parallel=parallel,
        FUN=SeqArray::seqMissing
    )
    if (per.variant)
        rowData(gdsfile)[["variant_missing"]] <- seqMissing
    else
        SummarizedExperiment::colData(gdsfile)[["sample_missing"]] <- seqMissing
    gdsfile
}
setMethod("seqMissing", "SummarizedExperiment", .seqMissing)        

.seqNumAllele <- function(gdsfile){
    numAllele <- .doCompatibleFunction(gdsfile, FUN=SeqArray::seqNumAllele)

    rowData(gdsfile)[["numAllele"]] <- numAllele
    gdsfile
}
setMethod("seqNumAllele", "SummarizedExperiment", .seqNumAllele)        

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
setMethod("hwe", "SummarizedExperiment", .hwe)

.inbreedCoeff <- function(gdsobj, margin=c("by.variant", "by.sample"),
                          use.names=FALSE){
    inbCoef <- .doCompatibleFunction(
        gdsobj, margin=margin, use.names=use.names, FUN = SeqVarTools::inbreedCoeff)
    inbCoef   ## returns a named (if use.names=TRUE) vector
}
setMethod("inbreedCoeff", "SummarizedExperiment", .inbreedCoeff)

.pca <- function(gdsobj, eigen.cnt=32){
    pca <- .doCompatibleFunction(gdsobj, eigen.cnt=eigen.cnt, FUN = SeqVarTools::pca)
    pca   ## returns a list, $eigenval (vector), $eigenvect (matrix)
}
setMethod("pca", "SummarizedExperiment", .pca)

.titv <- function(gdsobj, by.sample=FALSE, use.names=FALSE){
    titv <- .doCompatibleFunction(gdsobj, by.sample=by.sample,
                                  use.names=use.names, FUN = SeqVarTools::titv)
    titv   ## returns a scalar / vector (if by.sample=TRUE)
}
setMethod("titv", "SummarizedExperiment", .titv)

.refDosage <- function(gdsobj, use.names=TRUE){
    dos <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                  FUN = SeqVarTools::refDosage)
    dos   ## returns a matrix, nsamp * nvar
}
setMethod("refDosage", "SummarizedExperiment", .refDosage)

.altDosage <- function(gdsobj, use.names=TRUE){
    dos <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                  FUN = SeqVarTools::altDosage)
    dos   ## returns a matrix, nsamp * nvar
}
setMethod("altDosage", "SummarizedExperiment", .altDosage)

.ctSingleton <- function(gdsobj, use.names=FALSE){
    ct <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                  FUN = SeqVarTools::countSingletons)
    ct   ## returns a vector of the number of singleton variants per sample.
}
setMethod("countSingletons", "SummarizedExperiment", .ctSingleton)

.heterozygosity <- function(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE){
    margin <- match.arg(margin)
    hetero <- .doCompatibleFunction(gdsobj, margin=margin, use.names=use.names,
                                  FUN = SeqVarTools::heterozygosity)
    hetero   ## returns a vector of the number of singleton variants per sample.
}
setMethod("heterozygosity", "SummarizedExperiment", .heterozygosity)

.homozygosity <- function(gdsobj, allele=c("any", "ref", "alt"), margin=c("by.variant", "by.sample"), use.names=FALSE){
    allele <- match.arg(allele)
    margin <- match.arg(margin)
    homo <- .doCompatibleFunction(gdsobj, allele=allele, margin=margin,
                                  use.names=use.names,
                                  FUN = SeqVarTools::homozygosity)
    homo   ## returns a vector of the number of singleton variants per sample.
}
setMethod("homozygosity", "SummarizedExperiment", .homozygosity)

.meanBySample <- function(gdsobj, var.name, use.names=FALSE){
    mean <- .doCompatibleFunction(gdsobj, var.name=var.name, use.names=use.names,
                                  FUN = SeqVarTools::meanBySample)
    mean
} 
setMethod("meanBySample", "SummarizedExperiment", .meanBySample)

.missingGenotypeRate <- function(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE){
    margin <- match.arg(margin)
    misRate <- .doCompatibleFunction(gdsobj, margin=margin, use.names=use.names,
                                     FUN = SeqVarTools::missingGenotypeRate)
    misRate
}
setMethod("missingGenotypeRate", "SummarizedExperiment", .missingGenotypeRate)

.isSNV <- function(gdsobj, biallelic=TRUE){
    issnv <- .doCompatibleFunction(gdsobj, biallelic=biallelic,
                                   FUN = SeqVarTools::isSNV)
    issnv
}
setMethod("isSNV", "SummarizedExperiment", .isSNV)

.isVariant <- function(gdsobj, use.names=FALSE){
    isvar <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                   FUN = SeqVarTools::isVariant)
    isvar
}
setMethod("isVariant", "SummarizedExperiment", .isVariant)
