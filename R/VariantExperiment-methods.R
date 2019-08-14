.saveGDSMaybe <- function(gdsfile) {
    file <- gdsfile(gdsfile)
    gdsdim <- c(seqSummary(file, verbose=FALSE)$num.variant,
                seqSummary(file, verbose=FALSE)$num.sample)
    sedim <- dim(gdsfile) ## variants * samples
    if(!identical(sedim, gdsdim))
        stop("use 'saveVariantExperiment()'",
             " to synchronize on-disk and in-memory representations")
    gdsfile
}
        
.doCompatibleFunction <- function(gdsfile, ..., FUN) {
    gdsfile <- .saveGDSMaybe(gdsfile)
    file <- gdsfile(gdsfile)
    f <- SeqArray::seqOpen(file)
    on.exit(SeqArray::seqClose(f))
    FUN(f, ...)
}

.permDim <- function(object, ve) {
    identical(rev(dim(object)), dim(ve))
}

.seqAlleleFreq <- function(gdsfile, ref.allele=0L, .progress=FALSE,
                               parallel = seqGetParallel()){
    alleleFreq <- .doCompatibleFunction(
        gdsfile, ref.allele=ref.allele, .progress=.progress, parallel=parallel,
        FUN=SeqArray::seqAlleleFreq
    )
    alleleFreq
}

#' Statistical functions for \code{VariantExperiment} objects.
#' @name VariantExperiment-methods
#' @rdname VariantExperiment-methods
#' @aliases seqAlleleFreq seqAlleleFreq,VariantExperiment-method
#' @param gdsfile an \code{VariantExperiment} object that with
#'     synchronized gds file.
#' @param ref.allele a single numeric value, a numeric vector or a
#'     character vector; see \code{?SeqArray::seqAlleleFreq} for more
#'     details.
#' @param .progress Logical, show process information if \code{TRUE}.
#' @param parallel A logical value to indicate serial processing
#'     (\code{FALSE}) or multicore processing (\code{TRUE}). Takes
#'     numeric value or other value; see \code{?SeqArray::seqParallel}
#'     for more details.
#' @return Statistical results in \code{vector} or \code{data.frame} format. 
#' @export
#' @examples
#' gds <- SeqArray::seqExampleFileName("gds")
#' ve <- makeVariantExperimentFromGDS(gds)
#' ve
#' 
#' ## sample missing rate
#' mr.samp <- seqMissing(ve, per.variant = FALSE)
#' head(mr.samp)
#'
#' ## hwe
#' hwe <- hwe(ve)
#' head(hwe)
#'
#' ## titv ratio by sample / overall
#' titv <- titv(ve, by.sample=TRUE)
#' head(titv)
#' titv(ve, by.sample=FALSE)
#'
#' ## countSingletons
#' countSingletons(ve)

setMethod("seqAlleleFreq", "VariantExperiment", .seqAlleleFreq)        

.seqAlleleCount <- function(gdsfile, ref.allele=0L, .progress=FALSE, parallel = seqGetParallel()){
    alleleCount <- .doCompatibleFunction(
        gdsfile, ref.allele=ref.allele, .progress=.progress, parallel=parallel,
        FUN=SeqArray::seqAlleleCount
    )
    alleleCount
}
#' @rdname VariantExperiment-methods
#' @aliases seqAlleleCount seqAlleleCount,VariantExperiment-method
#' @export
setMethod("seqAlleleCount", "VariantExperiment", .seqAlleleCount)        

.seqMissing <- function(gdsfile, per.variant=TRUE, .progress=FALSE,
                        parallel = seqGetParallel())
{
    seqMissing <- .doCompatibleFunction(gdsfile,
                                        per.variant=per.variant,
                                        .progress =.progress,
                                        parallel=parallel,
                                        FUN=SeqArray::seqMissing )
    seqMissing
}

#' @rdname VariantExperiment-methods
#' @aliases seqMissing seqMissing,VariantExperiment-method
#' @param per.variant A logical value to indicate whether to calculate
#'     missing rate for variant (\code{TRUE}), or samples
#'     (\code{FALSE}).
#' @export
setMethod("seqMissing", "VariantExperiment", .seqMissing)        

.seqNumAllele <- function(gdsfile){
    numAllele <- .doCompatibleFunction(gdsfile,
                                       FUN=SeqArray::seqNumAllele)
    numAllele
}

#' @rdname VariantExperiment-methods
#' @aliases seqNumAllele,VariantExperiment-method
#' @export
setMethod("seqNumAllele", "VariantExperiment", .seqNumAllele)        

## SeqVarTools stats functions:

## already works: alleleFrequency(SeqArray::seqAlleleFreq),
## duplicateDiscordance, granges(SummarizedExperiment::granges),
## nAlleles(SeqArray::seqNumAllele), nAllele(SeqArray::seqAlleleCount)

## implemented: hwe, inbreedCoeff, pca, titv, refDosage, altDosage,
## countSingletons, homozygosity, heterozygosity, meanBySample, isSNV,
## isVariant

## removed:
## getGenotype/getGenotypeAlleles/expandedAltDosage/alleleDosage
## (SeqVarGDSClass,numeric)/alleleDosage(SeqVarGDSClass,list)

## Iterators(Extends ‘SeqVarData’ to provide iterators over variants.?)
## variantRanges(x): Get the variant ranges.
## variantFilter(x): Get the list of variant indices. 
## currentRanges(x): Get the ranges selected in the current iteration.
## duplicateDiscordance (SeqVarData),
## granges(already works for SE),

## Methods for class "SeqVarData" in package "SeqVarTools":
## SeqVarData(gds, sample.data)
## "sample.data" in class "AnnotatedDataFrame". 
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

#' @rdname VariantExperiment-methods
#' @aliases hwe,VariantExperiment-method
#' @param gdsobj same as above \code{gdsfile} argument.
#' @param permute A logical value indicating whether to permute the
#'     genotypes. See \code{?SeqVarTools::hwe} for more details.
#' @export
setMethod("hwe", "VariantExperiment", .hwe)

.inbreedCoeff <- function(gdsobj, margin=c("by.variant", "by.sample"),
                          use.names=FALSE){
    inbCoef <- .doCompatibleFunction(gdsobj, margin=margin,
                                     use.names=use.names,
                                     FUN = SeqVarTools::inbreedCoeff)
    inbCoef   ## returns a named (if use.names=TRUE) vector
}

#' @rdname VariantExperiment-methods
#' @aliases inbreedCoeff,VariantExperiment-method
#' @param margin "by.variant" OR "by.sample" to indicate
#'     whether the calculation should be done over all samples for
#'     each variant, or over all variants for each sample. See
#'     \code{?SeqVarTools::inbreedCoeff} for more details.
#' @param use.names A logical value indicating whether to assign
#'     variant or sample IDs as names of the output vector.
#' @export
setMethod("inbreedCoeff", "VariantExperiment", .inbreedCoeff)

.pca <- function(gdsobj, eigen.cnt=32){
    pca <- .doCompatibleFunction(gdsobj, eigen.cnt=eigen.cnt,
                                 FUN = SeqVarTools::pca)
    pca   ## returns a list, $eigenval (vector), $eigenvect (matrix)
}

#' @rdname VariantExperiment-methods
#' @aliases pca,VariantExperiment-method
#' @param eigen.cnt An integer value indicating how many eigenvalues
#'     and eignvectors to return. The default is 32.
#' @export
setMethod("pca", "VariantExperiment", .pca)

.titv <- function(gdsobj, by.sample=FALSE, use.names=FALSE){
    titv <- .doCompatibleFunction(gdsobj, by.sample=by.sample,
                                  use.names=use.names,
                                  FUN = SeqVarTools::titv)
    titv   ## returns a scalar / vector (if by.sample=TRUE)
}

#' @rdname VariantExperiment-methods
#' @aliases titv,variantExperiment-method
#' @param by.sample A logical value indicating whether TiTv should be
#'     calculated by sample or overall for the entire
#'     \code{VariantExperiment} object. See \code{?SeqVarTools::titv}
#'     for more details.
#' @export
setMethod("titv", "VariantExperiment", .titv)

.refDosage <- function(gdsobj, use.names=TRUE){
    dos <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                  FUN = SeqVarTools::refDosage)
    if (.permDim(dos, gdsobj))
        return(t(dos))   ## always return a matrix with same dimension
                         ## of ve
    dos
}

#' @rdname VariantExperiment-methods
#' @aliases refDosage,VariantExperiment-method
#' @export
setMethod("refDosage", "VariantExperiment", .refDosage)

.altDosage <- function(gdsobj, use.names=TRUE, sparse=FALSE){
    dos <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                 sparse=sparse,
                                 FUN = SeqVarTools::altDosage)
    if (.permDim(dos, gdsobj))
        return(t(dos))   ## always return a matrix with same dimension
                         ## of ve
    dos  
}

#' @rdname VariantExperiment-methods
#' @aliases altDosage,VariantExperiment-method
#' @param sparse A Logical value indicating whether or not to return
#'     the alterate allele dosage as a sparse matrix. In most cases,
#'     it will dramatically reduce the size of the returned
#'     object. See \code{?SeqVarTools::altDosage} for more details.
#' @export
setMethod("altDosage", "VariantExperiment", .altDosage)

.ctSingleton <- function(gdsobj, use.names=FALSE){
    ct <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                  FUN = SeqVarTools::countSingletons)
    ct   ## returns a vector of the number of singleton variants per
         ## sample.
}

#' @rdname VariantExperiment-methods
#' @aliases countSingletons,VariantExperiment-method
#' @export
setMethod("countSingletons", "VariantExperiment", .ctSingleton)

.heterozygosity <- function(gdsobj,
                            margin=c("by.variant", "by.sample"),
                            use.names=FALSE){
    margin <- match.arg(margin)
    hetero <- .doCompatibleFunction(gdsobj, margin=margin,
                                    use.names=use.names,
                                    FUN = SeqVarTools::heterozygosity)
    hetero   ## returns a vector of the number of singleton variants per sample.
}

#' @rdname VariantExperiment-methods
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
    homo   ## returns a vector of the number of singleton variants per
           ## sample.
}

#' @rdname VariantExperiment-methods
#' @aliases homozygosity,VariantExperiment-method
#' @param allele Choose from "any", "ref," or "alt," to indicate which
#'     alleles to consider when calculating homozygosity. See
#'     \code{?SeqVarTools::homozygosity} for more details.

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

#' @rdname VariantExperiment-methods
#' @aliases meanBySample,VariantExperiment-method
#' @param var.name Character string with name of the variable. Choose
#'     from \code{names(assays(VE_Object))}. See
#'     \code{?SeqVarTools::meanBySample} for more details.
#' @export
setMethod("meanBySample", "VariantExperiment", .meanBySample)

.isSNV <- function(gdsobj, biallelic=TRUE){
    issnv <- .doCompatibleFunction(gdsobj, biallelic=biallelic,
                                   FUN = SeqVarTools::isSNV)
    issnv
}

#' @rdname VariantExperiment-methods
#' @aliases isSNV,VariantExperiment-method
#' @param biallelic A logical indicating whether to only consider
#'     biallelic SNVs. See \code{?SeqVarTools::isSNV} for more
#'     details.
#' @export
setMethod("isSNV", "VariantExperiment", .isSNV)

.isVariant <- function(gdsobj, use.names=FALSE){
    isvar <- .doCompatibleFunction(gdsobj, use.names=use.names,
                                   FUN = SeqVarTools::isVariant)
    if (.permDim(isvar, gdsobj))
        return(t(isvar))   ## always return a matrix with same
                           ## dimension of ve
    isvar
}

#' @rdname VariantExperiment-methods
#' @aliases isVariant,VariantExperiment-method
#' @export
setMethod("isVariant", "VariantExperiment", .isVariant)
