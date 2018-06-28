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

.permDim <- function(object, ve) {
    if (identical(rev(dim(object)), dim(ve))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

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

#' Statistical functions for \code{VariantExperiment} objects.
#' @name seqAlleleFreq
#' @rdname VariantExperiment-stats
#' @aliases seqAlleleFreq,VariantExperiment-method
#' @param gdsfile xx
#' @param ref.allele xx
#' @param .progress xx 
#' @param parallel xx
#' @export
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
#' @name seqAlleleCount
#' @rdname VariantExperiment-stats
#' @aliases seqAlleleCount,VariantExperiment-method
#' @export
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
#' @name seqMissing
#' @rdname VariantExperiment-stats
#' @aliases seqMissing,VariantExperiment-method
#' @param per.variant xx
#' @export
setMethod("seqMissing", "VariantExperiment", .seqMissing)        

.seqNumAllele <- function(gdsfile){
    numAllele <- .doCompatibleFunction(gdsfile, FUN=SeqArray::seqNumAllele)
    numAllele
    ## rowData(gdsfile)[["numAllele"]] <- numAllele
    ## gdsfile
}
#' @name seqNumAllele
#' @rdname VariantExperiment-stats
#' @aliases seqNumAllele,VariantExperiment-method
#' @export
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

#' @importMethodsFrom SeqVarTools hwe inbreedCoeff pca refDosage
#'     altDosage countSingletons homozygosity heterozygosity
#'     meanBySample titv isSNV isVariant

.hwe <- function(gdsobj, permute=FALSE){
    hwe <- .doCompatibleFunction(gdsobj, permute=permute,
                                 FUN = SeqVarTools::hwe)
    hwe   ## returns a data.frame, will output as is, do not paste
          ## into rowData(se)
}

#' @name hwe
#' @rdname VariantExperiment-stats
#' @aliases hwe,VariantExperiment-method
#' @param gdsobj xx
#' @param permute xx
#' @export
setMethod("hwe", "VariantExperiment", .hwe)

.inbreedCoeff <- function(gdsobj, margin=c("by.variant", "by.sample"),
                          use.names=FALSE){
    inbCoef <- .doCompatibleFunction(
        gdsobj, margin=margin, use.names=use.names, FUN = SeqVarTools::inbreedCoeff)
    inbCoef   ## returns a named (if use.names=TRUE) vector
}

#' @name inbreedCoeff
#' @rdname VariantExperiment-stats
#' @aliases inbreedCoeff,VariantExperiment-method
#' @param margin xx
#' @param use.names xx
#' @export
setMethod("inbreedCoeff", "VariantExperiment", .inbreedCoeff)

.pca <- function(gdsobj, eigen.cnt=32){
    pca <- .doCompatibleFunction(gdsobj, eigen.cnt=eigen.cnt, FUN = SeqVarTools::pca)
    pca   ## returns a list, $eigenval (vector), $eigenvect (matrix)
}

#' @name pca
#' @rdname VariantExperiment-stats
#' @aliases pca,VariantExperiment-method
#' @param eigen.cnt xx
#' @export
setMethod("pca", "VariantExperiment", .pca)

.titv <- function(gdsobj, by.sample=FALSE, use.names=FALSE){
    titv <- .doCompatibleFunction(gdsobj, by.sample=by.sample,
                                  use.names=use.names, FUN = SeqVarTools::titv)
    titv   ## returns a scalar / vector (if by.sample=TRUE)
}

#' @name titv
#' @rdname VariantExperiment-stats
#' @aliases titv,variantExperiment-method
#' @param by.sample xx
#' @export
setMethod("titv", "VariantExperiment", .titv)

.refDosage <- function(gdsobj, use.names=TRUE){
    dos <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                  FUN = SeqVarTools::refDosage)
    if (.permDim(dos, gdsobj))
        return(t(dos))   ## always return a matrix with same dimension of ve
    dos
}

#' @name refDosage
#' @rdname VariantExperiment-stats
#' @aliases refDosage,VariantExperiment-method
#' @export
setMethod("refDosage", "VariantExperiment", .refDosage)

.altDosage <- function(gdsobj, use.names=TRUE, sparse=FALSE){
    dos <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                 sparse=sparse,
                                 FUN = SeqVarTools::altDosage)
    if (.permDim(dos, gdsobj))
        return(t(dos))   ## always return a matrix with same dimension of ve
    dos  
}

#' @name altDosage
#' @rdname VariantExperiment-stats
#' @aliases altDosage,VariantExperiment-method
#' @param sparse xx
#' @export
setMethod("altDosage", "VariantExperiment", .altDosage)

.ctSingleton <- function(gdsobj, use.names=FALSE){
    ct <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                  FUN = SeqVarTools::countSingletons)
    ct   ## returns a vector of the number of singleton variants per sample.
}

#' @name countSingletons
#' @rdname VariantExperiment-stats
#' @aliases countSingletons,VariantExperiment-method
#' @export
setMethod("countSingletons", "VariantExperiment", .ctSingleton)

.heterozygosity <- function(gdsobj,
                            margin=c("by.variant", "by.sample"),
                            use.names=FALSE){
    margin <- match.arg(margin)
    hetero <- .doCompatibleFunction(gdsobj, margin=margin, use.names=use.names,
                                  FUN = SeqVarTools::heterozygosity)
    hetero   ## returns a vector of the number of singleton variants per sample.
}

#' @name heterozygosity
#' @rdname VariantExperiment-stats
#' @aliases heterozygosity,VariantExperiment-method
#' @export
setMethod("heterozygosity", "VariantExperiment", .heterozygosity)

.homozygosity <- function(gdsobj, allele=c("any", "ref", "alt"),
                          margin=c("by.variant", "by.sample"),
                          use.names=FALSE){
    allele <- match.arg(allele)
    margin <- match.arg(margin)
    homo <- .doCompatibleFunction(gdsobj, allele=allele, margin=margin,
                                  use.names=use.names,
                                  FUN = SeqVarTools::homozygosity)
    homo   ## returns a vector of the number of singleton variants per sample.
}

#' @name homozygosity
#' @rdname VariantExperiment-stats
#' @aliases homozygosity,VariantExperiment-method
#' @param allele xx
#' @export
setMethod("homozygosity", "VariantExperiment", .homozygosity)

.meanBySample <- function(gdsobj, var.name, use.names=FALSE){
    if (! var.name %in% names(assays(gdsobj))) {
        stop(paste0('the "var.names" argument must take values from ',
                    paste(names(assays(gdsobj)), collapse="")))
    }
    var.name <- sub("/data", "", var.name)
    mean <- .doCompatibleFunction(gdsobj, var.name=var.name,
                                  use.names=use.names,
                                  FUN = SeqVarTools::meanBySample)
    mean
} 

#' @name meanBySample
#' @rdname VariantExperiment-stats
#' @aliases meanBySample,VariantExperiment-method
#' @param var.name xx
#' @export
setMethod("meanBySample", "VariantExperiment", .meanBySample)

.isSNV <- function(gdsobj, biallelic=TRUE){
    issnv <- .doCompatibleFunction(gdsobj, biallelic=biallelic,
                                   FUN = SeqVarTools::isSNV)
    issnv
}

#' @name isSNV
#' @rdname VariantExperiment-stats
#' @aliases isSNV,VariantExperiment-method
#' @param biallelic xx
#' @export
setMethod("isSNV", "VariantExperiment", .isSNV)

.isVariant <- function(gdsobj, use.names=FALSE){
    isvar <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                   FUN = SeqVarTools::isVariant)
    if (.permDim(isvar, gdsobj))
        return(t(isvar))   ## always return a matrix with same dimension of ve
    isvar
}

#' @name isVariant
#' @rdname VariantExperiment-stats
#' @aliases isVariant,VariantExperiment-method
#' @export
setMethod("isVariant", "VariantExperiment", .isVariant)
